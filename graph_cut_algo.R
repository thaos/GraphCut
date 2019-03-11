get_adjacent_gridpoints <- function(idx, ni = 8, nj = 32){
  i <- idx %% ni
  i <- ifelse(i == 0, ni, i)
  j <- (idx - 1) %/% ni + 1
  get_adjacent_gridpoints_arrayInd <- function(i, j, ni = 8, nj = 32){
    ans <- data.frame(
      i = c(i - 1, i + 1, i , i),
      j = c(j, j, j - 1, j + 1)
    )
    subset(ans, (i >= 1 & i <= ni & j >= 1 & j <= nj))
  }
  adj_arrayInd <- get_adjacent_gridpoints_arrayInd(i, j, ni = ni, nj = nj)
  adj <- adj_arrayInd$i + ni * (adj_arrayInd$j - 1)
  return(adj)
}
get_adjacent_gridpoints(1, ni = 8, nj = 5)
get_adjacent_gridpoints(8, ni = 8, nj = 5)
get_adjacent_gridpoints(8*5, ni = 8, nj = 5)
get_adjacent_gridpoints(8*5 - 7, ni = 8, nj = 5)
get_adjacent_gridpoints(20, ni = 8, nj = 5)

create_arcs <- function(diff_overlap, overlap_id, canvas){
  for(from in seq_along(diff_overlap)){
    adj <- get_adjacent_gridpoints(
      overlap_id[from], ni = nrow(canvas), nj = ncol(canvas)
    )
    adj <- adj[adj %in% overlap_id]
    to <- vapply(adj, function(x) which(overlap_id == x), FUN.VALUE = 1)
    weight <- abs(diff_overlap[from]) + abs(diff_overlap[to])
    if(from == 1){
      arcs <- data.frame(from = from, to = to, weight = weight)
    } else {
      arcs <- rbind(arcs, data.frame(from = from, to = to, weight = weight))
    }
  }
  subset(arcs, from < to)
}

# Add Source and Sink#
add_sink_and_source <- function(arcs, overlap_id, source_nodes, sink_nodes){
  rbind(
    arcs,
    data.frame(from = length(overlap_id) + 1, to = source_nodes, weight = 10E6),
    data.frame(from = length(overlap_id) + 2, to = sink_nodes, weight = 10E6)
  )
}

reorder_nodes <- function(arcs){
  arcs <- arcs[order(arcs[, 1], arcs[, 2]), ]
}

rearrange_node <- function(cutset){
  cutset <- apply(
    cutset, 1, 
    function(x){
      x <- as.vector(x)
      if(x[1] > x[2]){
        ans <- x[c(2, 1, 3)]
      } else {
        ans <- x
      }
      return(ans)
    }
  ) %>%
    t()
  colnames(cutset) <- colnames(cutset)
  return(cutset)
}

local_to_global_arcs <- function(arcs_local, overlap_id){
  arcs_global <- arcs_local
  arcs_global[, 1:2] <- overlap_id[c(arcs_local[, 1:2])]
  return(arcs_global)
}

add_patch_to_list <- function(patch_list, number, patch, patch_id){
  rbind(
    patch_list,
    data.frame(number = number, patch = I(list(patch)), patch_id = I(list(patch_id)))
  )
}

assemble_overlap_origin <- function(mincut, label,  overlap_origin, overlap_id){
  overlaped_origin <- overlap_origin
  tcut <- mincut$t.cut[mincut$t.cut <= length(overlap_id)]
  overlaped_origin[tcut] <- label
  return(overlaped_origin)
}


assemble_overlap <- function(mincut, overlap_A, overlap_B, overlap_id){
  overlaped <- overlap_A
  tcut <- mincut$t.cut[mincut$t.cut <= length(overlap_id)]
  overlaped[tcut] <- overlap_B[tcut]
  return(overlaped)
}


# adding overlapes to the canvas
update_canvas <- function(canvas, patch, patch_id, overlaped = NULL, overlaped_id = NULL){
  canvas[patch_id] <- patch
  if(!is.null(overlaped) & !is.null(overlaped_id)){
    canvas[overlaped_id] <- overlaped
  }
  return(canvas)
}

# adding overlapes to the canvas
update_canvas_origin <- function(canvas_origin, label, patch_id, overlaped_origin = NULL, overlaped_id = NULL){
  canvas_origin[patch_id] <- label
  if(!is.null(overlaped_origin) & !is.null(overlaped_id)){
    canvas_origin[overlaped_id] <- overlaped_origin
  }
  return(canvas_origin)
}
lines_seams <- function(cutset_global, canvas){
  n1_coord <- arrayInd(cutset_global[, 1], .dim = dim(canvas))
  n2_coord <- arrayInd(cutset_global[, 2], .dim = dim(canvas))
  ncolors = 10
  # wcuts <- cut(cutset_global[, 3], breaks = quantile(cutset_global[,3], probs = seq(0, 1, length.out = ncolors + 1))) 
  if(length(unique(cutset_global[, 3])) < ncolors){
    ncolors <- length(unique(cutset_global[, 3]))
  }
  if(length(unique(cutset_global[, 3])) == 1){
    wcuts <- rep(1, nrow(cutset_global))
  } else {
    wcuts <- cut(cutset_global[, 3], breaks = ncolors) 
  }
  pts <- matrix(NA, ncol = 2, nrow = 2)
  for(i in seq.int(nrow(cutset_global))){
    if(n1_coord[i, 1]  == n2_coord[i, 1]){
      pts[, 2] = mean( c(n1_coord[i, 2], n2_coord[i, 2]))
      pts[, 1] = n1_coord[i, 1] + c(-0.5, 0.5)
    } else {
      pts[, 1] = mean( c(n1_coord[i, 1], n2_coord[i, 1]))
      pts[, 2] = n1_coord[i, 2] + c(-0.5, 0.5)
    }
    lines(pts, col = rev(heat.colors(ncolors))[wcuts[i]], lwd = 2)
  }
}

check_oldseam_in_overlap <- function(cut, overlap_id){
  (cut[1] %in% overlap_id & cut[2] %in% overlap_id)
}

compute_oldseam_weight1 <- function(cut, overlap, overlap_id, canvas, canvas_origin, patch_list){
  val_cut1 <- canvas[cut[1]]
  overlap1 <- canvas_origin[cut[1]]
  ioverlap1 <- which(patch_list$number == overlap1)
  val_cut2 <- patch_list$patch[[ioverlap1]][patch_list$patch_id[[ioverlap1]] == cut[2]]
  abs(val_cut1 - overlap[overlap_id == cut[1]]) +  abs(val_cut2 - overlap[overlap_id == cut[2]])
}
# debug(compute_oldseam_weight1)
compute_oldseam_weight2 <- function(cut, overlap, overlap_id, canvas, canvas_origin, patch_list){
  val_cut2 <- canvas[cut[2]]
  overlap2 <- canvas_origin[cut[2]]
  ioverlap2 <- which(patch_list$number == overlap2)
  val_cut1 <- patch_list$patch[[ioverlap2]][patch_list$patch_id[[ioverlap2]] == cut[1]]
  abs(val_cut2 - overlap[overlap_id == cut[2]]) +  abs(val_cut1 - overlap[overlap_id == cut[1]])
}
# debug(compute_oldseam_weight1)

add_oldcut <- function(cut, overlap, overlap_id, canvas, canvas_origin, patch_list, arcs){
  overlap_id_local <- seq.int(length(overlap_id))
  node1 <- overlap_id_local[overlap_id == cut[1]] 
  node2 <- overlap_id_local[overlap_id == cut[2]]
  iseam <- which(
    node1 == arcs[, 1] & node2 == arcs[, 2] |
      node1 == arcs[, 2] & node2 == arcs[, 1]
  )
  nodeseam_id <- max(arcs[, 1:2]) + 1
  arcs[iseam, 1:2] <- c(length(overlap_id) + 2, nodeseam_id)
  arcs <- rbind(
    arcs,
    c(node1, nodeseam_id, compute_oldseam_weight1(cut = cut, overlap = overlap, overlap_id = overlap_id, canvas = canvas, canvas_origin = canvas_origin, patch_list = patch_list)),
    c(node2, nodeseam_id, compute_oldseam_weight2(cut = cut, overlap = overlap, overlap_id = overlap_id, canvas = canvas, canvas_origin = canvas_origin, patch_list = patch_list))
  )
  return(arcs)
}

update_graph_with_cuts <- function(cutset_global, overlap, overlap_id, canvas, canvas_origin, patch_list, arcs){
  for(i in seq.int(nrow(cutset_global))){
    if(check_oldseam_in_overlap(cut = cutset_global[i, ], overlap_id = overlap_id)){
      arcs <- add_oldcut(
        cut = cutset_global[i, ],
        overlap = overlap, overlap_id = overlap_id,
        canvas = canvas, canvas_origin = canvas_origin,
        patch_list = patch_list, arcs = arcs
      )
    }
  }
  return(arcs)
}

update_cutset_oneseamnode <- function(seam_node, arcs, cutset_local, cutset_global, overlap_id){
  cutset_local <- rearrange_node(cutset_local)
  print(dim(cutset_global))
  cutset_global <- rearrange_node(cutset_global)
  print(dim(cutset_global))
  arcs_with_seamnode <- subset(arcs, from == seam_node | to == seam_node)
  stopifnot(nrow(arcs_with_seamnode) == 3)
  arcs_with_seamnode <- reorder_nodes(arcs_with_seamnode)                                   
  edge_local_coord <- arcs_with_seamnode[1:2, 1]
  # Case with a cut between the sink and  a seam node
  icutset_local <- which(
    apply(cutset_local, 1, function(x) all(x == arcs_with_seamnode[3, 1:2]))
  )
  if(length(icutset_local) == 1){
    return(cutset_global)
  } 
  edge_global_coord <- overlap_id[edge_local_coord]
  icutset_global <- which(
    apply(cutset_global[, 1:2], 1, function(x) all(x == edge_global_coord))
  )
  # Cases with a cut between a regular node and a seam node
  for(irow in 1:2){
    icutset_local <- which(
      apply(cutset_local, 1, function(x) all(x == arcs_with_seamnode[irow, 1:2]))
    )
    if(length(icutset_local) == 1){
      cutset_global[icutset_global, 3] <- cutset_local[icutset_local, 3]
      return(cutset_global)
    }
  }
  return(cutset_global[-icutset_global, ])
}  


remove_extra_nodes <- function(cutset_local, overlap_id){
  irm <- apply(cutset_local[, 1:2], 1, function(x) any(x > length(overlap_id)))
  return(cutset_local[!irm, ])
}

update_cutset_global <- function(arcs, cutset_local, cutset_global, overlap_id) {
  seam_nodes <- sort(unique(unlist(arcs[, 1:2])))
  seam_nodes <- seam_nodes[seam_nodes > (length(overlap_id) + 2)] 
  for (node in seam_nodes){
    cutset_global <- update_cutset_oneseamnode(
      seam_node = node,
      arcs = arcs,
      cutset_local = cutset_local,
      cutset_global = cutset_global,
      overlap_id = overlap_id
    )
  }
  cutset_global <- rbind(
    cutset_global,
    local_to_global_arcs(remove_extra_nodes(cutset_local, overlap_id), overlap_id)
  )
  return(cutset_global)
}

prepare_for_new_patch <- function(xstart, ystart, xlength = 32, ylength = 32 , canvas, canvas_id){
  patch <- canvas[seq(xstart, xstart + xlength - 1),  seq(ystart, ystart + ylength - 1)]
  patch_id <- canvas_id[seq(xstart, xstart + xlength - 1),  seq(ystart, ystart + ylength - 1)]
  overlap <- patch[!is.na(patch)]
  overlap_id <- patch_id[!is.na(patch)]
  overlap_arrInd <- arrayInd(overlap_id, .dim = dim(canvas))
  overlap_box <- canvas[min(overlap_arrInd[, 1]):max(overlap_arrInd[, 1]),
                        min(overlap_arrInd[, 2]):max(overlap_arrInd[, 2])]
  overlap_box_id <- canvas_id[min(overlap_arrInd[, 1]):max(overlap_arrInd[, 1]),
                              min(overlap_arrInd[, 2]):max(overlap_arrInd[, 2])]
  #find source node
  get_frame <- function(matrix){
    frame <- matrix(NA, ncol = ncol(matrix), nrow = nrow(matrix))
    frame[1,] <- 1
    frame[nrow(frame),] <- 1
    frame[, 1] <- 1
    frame[, ncol(frame)] <- 1
    return(frame * matrix)
  }
  adjacent_to_NA <- sapply(
    overlap_id,
    function(id){
      any(is.na(
        canvas[get_adjacent_gridpoints(id, ni = nrow(canvas), nj = ncol(canvas))]
      ))
    }
  )
  sink_id <- overlap_id[adjacent_to_NA]
  overlap_frame_id <- get_frame(overlap_box_id)
  overlap_frame_id <- overlap_frame_id[!is.na(overlap_frame_id)]
  patch_frame <- get_frame(patch)
  source_id <- patch_id[!is.na(patch_frame)]
  source_id <- source_id[source_id %in% overlap_frame_id]
  source_id <- source_id[!(source_id %in% sink_id)]
  stopifnot(length(source_id) > 1)
  if(length(sink_id) == 0){
    xcenter <- nrow(overlap_box_id)/2 + 0.5
    xcenter <- unique(c(floor(xcenter), ceiling(xcenter)))
    print(xcenter)
    ycenter <- ncol(overlap_box_id)/2 + 0.5
    ycenter <- unique(c(floor(ycenter), ceiling(ycenter)))
    print(xcenter)
    centers <- as.matrix(expand.grid(x = xcenter, y = ycenter))
    print(centers)
    sink_id <- apply(
      centers, 1, 
      function(center) overlap_box_id[center[1], center[2]]
    )
  }
  par(mfrow = c(2, 1))
  image(matrix(patch_id %in% source_id, ncol = ncol(patch_id)))
  image(matrix(patch_id %in% sink_id, ncol = ncol(patch_id)))
  source_id_local <- which(overlap_id %in% source_id)
  sink_id_local <- which(overlap_id %in% sink_id)
  list(
    patch = patch, patch_id = patch_id,
    overlap = overlap, overlap_id = overlap_id,
    source_id = source_id, sink_id = sink_id,
    source_id_local = source_id_local, sink_id_local = sink_id_local
  ) 
}


find_matching_patch <- function(patch_id, overlap, overlap_id, training_img){
  dim_patch <- dim(patch_id)
  dim_train_img <- dim(training_img)
  
  ij_scan <- expand.grid(
    i = 1:(dim_train_img[1] - dim_patch[1] + 1),
    j = 1:(dim_train_img[2] - dim_patch[2] + 1)
  )
  ioverlap <- patch_id %in% overlap_id
  compute_dist <- function(i, j, ioverlap){
    candidate_patch <- img[i:(i + dim_patch[1] - 1), j:(j + dim_patch[2] - 1)]
    candidate_overlap <- candidate_patch[ioverlap]
    sum(abs(overlap - candidate_overlap))
  }
  dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j, MoreArgs = list(ioverlap = ioverlap)) 
  dmat[dmat == 0] <- Inf
  imin <- which.min(dmat)
  patch_new <- with(ij_scan[imin, ], img[i:(i + dim_patch[1] - 1), j:(j + dim_patch[2] - 1)])
  overlap_new <- patch_new[ioverlap]
  list(patch_new = patch_new, overlap_new = overlap_new)
}


# main function -----------------------------------------------------------

add_newpatch <- function(
  xstart, ystart, xlength = 32, ylength = 32 ,
  canvas, canvas_origin, canvas_id,
  training_img,
  patch_list = NULL,  cutset_global = NULL
){
  preparations <- prepare_for_new_patch(
    xstart = xstart, ystart = ystart,
    xlength = xlength, ylength = ylength,
    canvas = canvas, canvas_id = canvas_id
  )
  matching_patch <- find_matching_patch(
    patch_id = preparations$patch_id,
    overlap = preparations$overlap,
    overlap_id = preparations$overlap_id,
    training_img = training_img
  )
  patch_id <- preparations$patch_id
  overlap_old <- preparations$overlap
  overlap_id <- preparations$overlap_id
  patch_new <- matching_patch$patch_new
  overlap_new <- matching_patch$overlap_new
  if(is.null(patch_list)){
    ina_canvas <- is.na(canvas)
    patch_old <- canvas[!ina_canvas]
    patch_old_id <- canvas_id[!ina_canvas]
    patch_list <- data.frame(
      number = 1,
      patch = I(list(patch_old)),
      patch_id = I(list(patch_old_id))
    )
  }
  patch_list <- add_patch_to_list(
    patch_list = patch_list, number = max(patch_list$number) + 1,
    patch = patch_new, patch_id = patch_id
  )
  
  overlap_diff <- abs(overlap_old -  overlap_new)
  arcs <- create_arcs(diff_overlap = overlap_diff, overlap_id = overlap_id, canvas = canvas)
  arcs <- reorder_nodes(arcs)
  graph <- graph_from_data_frame(arcs, directed = FALSE)
  plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))
  
  arcs <- add_sink_and_source(
    arcs = arcs, overlap_id = overlap_id,
    source_nodes = preparations$source_id_local,
    sink_nodes = preparations$sink_id_local
  )
  arcs <- reorder_nodes(arcs)
  graph <- graph_from_data_frame(arcs, directed = FALSE)
  plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))
  
  # Check if old seam in overlap
  if(!is.null(cutset_global)){
    cutset_global <- reorder_nodes(cutset_global)
    arcs <- update_graph_with_cuts(
      cutset_global = cutset_global,
      overlap = overlap_new, overlap_id = overlap_id,
      canvas = canvas, canvas_origin =  canvas_origin,
      patch_list = patch_list, arcs = arcs
    )
    arcs <- reorder_nodes(arcs)
  }
  
  # Solve MaxFlow / ST-MinCut
  mincut <- findMinCut(
    nodes = sort(unique(unlist(arcs[, 1:2]))) ,
    arcs = as.matrix(arcs),
    algorithm = "Ford-Fulkerson",
    source.node = length(overlap_old) + 1,
    sink.node = length(overlap_old) + 2,
    directed = FALSE
  )
  
  cutset_local <- mincut$cut.set
  cutset_local <- rearrange_node(cutset_local)
  cutset_local <- reorder_nodes(cutset_local)
  # debug(update_cutset_oneseamnode)
  if(is.null(cutset_global)){
    cutset_global <- local_to_global_arcs(cutset_local, overlap_id)
    cutset_global <- rearrange_node(cutset_global)
  } else{
    cutset_global <- update_cutset_global(
      arcs = arcs, cutset_local = cutset_local,
      cutset_global =  cutset_global, overlap_id = overlap_id
    )
  }
  
  overlaped <- assemble_overlap(
    mincut = mincut, overlap_A = overlap_old,
    overlap_B = overlap_new, overlap_id = overlap_id
  )
  overlaped_origin <- overlap_old
  overlaped_origin[] <- canvas_origin[overlap_id]
  overlaped_origin <- assemble_overlap_origin(
    mincut = mincut, label = max(patch_list$number) + 1,
    overlap_origin = overlaped_origin, overlap_id = overlap_id
  )  
  canvas <- update_canvas(
    canvas = canvas,
    patch = patch_new, patch_id = patch_id,
    overlaped = overlaped, overlaped_id = overlap_id
  ) 
  canvas_origin <- update_canvas_origin(
    canvas_origin = canvas_origin,
    label = max(patch_list$number) + 1,
    patch_id = patch_id,
    overlaped_origin = overlaped_origin,
    overlaped_id = overlap_id
  )  
  list(
    canvas = canvas, canvas_origin = canvas_origin, canvas_id = canvas_id,
    overlap_id = overlap_id, patch_list = patch_list,
    cutset_global = cutset_global
  )
}
