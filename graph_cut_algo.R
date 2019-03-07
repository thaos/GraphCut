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

create_egdes <- function(idx, diff_patch){
  adj <- get_adjacent_gridpoints(idx, ni = nrow(diff_patch), nj = ncol(diff_patch))
  # adj <- adj[adj > idx]
  weight <- diff_patch[idx] + diff_patch[adj]
  data.frame(from = idx, to = adj, weight = weight)
}

add_source_to_left <- function(patch_id){
  idx_source <- 1 + (seq.int(ncol(patch_id)) - 1) * nrow(patch_id)
  idx_sink <-  seq.int(ncol(patch_id)) * nrow(patch_id) 
  arcs <- rbind( 
    data.frame(from = length(patch_id) + 1, to = idx_source, weight = 10E6),
    data.frame(from = length(patch_id) + 2, to = idx_sink, weight = 10E6)
  )
  return(arcs)
}

create_arcs <- function(diff_patch){
  arcs <- lapply(
    1:length(diff_patch),
    create_egdes,
    diff_patch = diff_patch
  ) %>%
    do.call(rbind, .) %>%
    subset(from < to)
}

# Add Source and Sink#
add_sink_and_source <- function(arcs, patch_id, add_source = add_source_to_left){
  arcs_with_ST <- rbind(
    arcs,
    add_source(patch_id)
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

local_to_global_arcs <- function(arcs_local, patch_id){
  arcs_global <- arcs_local
  arcs_global[, 1:2] <- patch_id[c(arcs_local[, 1:2])]
  return(arcs_global)
}

assemble_patch <- function(mincut, patch_A, patch_B, patch_id){
  patched <- patch_A
  tcut <- mincut$t.cut[mincut$t.cut <= length(patch_id)]
  patched[tcut] <- patch_B[tcut]
  return(patched)
}

assemble_opposite_patch <- function(mincut, patch_A, patch_B, patch_id){
  opposite_patched <- patch_B
  tcut <- mincut$t.cut[mincut$t.cut <= length(patch_id)]
  opposite_patched[tcut] <- patch_A[tcut]
  return(opposite_patched)
}

# adding patches to the canvas
update_canvas <- function(canvas, img_B, img_B_id, patched, patched_id){
  canvas[img_B_id] <- img_B
  canvas[patched_id] <- patched
  return(canvas)
}

update_opposite_canvas <- function(opposite_canvas, opposite_patched, patched_id){
  opposite_canvas[patched_id] <- opposite_patched
  return(opposite_canvas)
}

lines_seams <- function(cutset_global, canvas){
  n1_coord <- arrayInd(cutset_global[, 1], .dim = dim(canvas))
  n2_coord <- arrayInd(cutset_global[, 2], .dim = dim(canvas))
  ncolors = 10
  # wcuts <- cut(cutset_global[, 3], breaks = quantile(cutset_global[,3], probs = seq(0, 1, length.out = ncolors + 1))) 
  wcuts <- cut(cutset_global[, 3], breaks = ncolors) 
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

check_oldseam_in_patch <- function(cut, patch_id){
  (cut[1] %in% patch_id & cut[2] %in% patch_id)
}

compute_oldseam_weight1 <- function(cut, patch, patch_id, canvas, opposite_canvas){
  abs(canvas[cut[1]] - patch[patch_id == cut[1]]) + abs(opposite_canvas[cut[2]] - patch[patch_id == cut[2]])
}
compute_oldseam_weight2 <- function(cut, patch, patch_id, canvas, opposite_canvas){
  abs(canvas[cut[2]] - patch[patch_id == cut[2]]) + abs(opposite_canvas[cut[1]] - patch[patch_id == cut[1]])
}

add_oldcut <- function(cut, patch, patch_id, canvas, opposite_canvas, arcs){
  patch_id_local <- matrix(
    seq.int(length(patch_id)),
    nrow = nrow(patch_id),
    ncol = ncol(patch_id)
  )
  node1 <- patch_id_local[patch_id == cut[1]] 
  node2 <- patch_id_local[patch_id == cut[2]]
  iseam <- which(
    node1 == arcs[, 1] & node2 == arcs[, 2] |
      node1 == arcs[, 2] & node2 == arcs[, 1]
  )
  nodeseam_id <- max(arcs[, 1:2]) + 1
  arcs[iseam, 1:2] <- c(length(patch_id) + 2, nodeseam_id)
  arcs <- rbind(
    arcs,
    c(node1, nodeseam_id, compute_oldseam_weight1(cut = cut, patch = patch, patch_id = patch_id, canvas = canvas, opposite_canvas = opposite_canvas)),
    c(node2, nodeseam_id, compute_oldseam_weight2(cut = cut, patch = patch, patch_id = patch_id, canvas = canvas, opposite_canvas = opposite_canvas))
  )
  return(arcs)
}

update_graph_with_cuts <- function(cutset_global, patch, patch_id, canvas, opposite_canvas, arcs){
  for(i in seq.int(nrow(cutset_global))){
    if(check_oldseam_in_patch(cut = cutset_global[i, ], patch_id = patch_AB_id)){
      arcs <- add_oldcut(
        cut = cutset_global[i, ],
        patch = patch, patch_id = patch_id,
        canvas = canvas, opposite_canvas = opposite_canvas,
        arcs = arcs
      )
    }
  }
  return(arcs)
}

update_cutset_oneseamnode <- function(seam_node, arcs, cutset_local, cutset_global, patch_id){
  cutset_local <- rearrange_node(cutset_local)
  print(dim(cutset_global))
  cutset_global <- rearrange_node(cutset_global)
  print(dim(cutset_global))
  arcs_with_seamnode <- subset(arcs, from == seam_node | to == seam_node)
  arcs_with_seamnode <- reorder_nodes(arcs_with_seamnode)                                   
  edge_local_coord <- arcs_with_seamnode[1:2, 1]
  # Case with a cut between the sink and  a seam node
  icutset_local <- which(
    apply(cutset_local, 1, function(x) all(x == arcs_with_seamnode[3, 1:2]))
  )
  if(length(icutset_local) == 1){
    return(cutset_global)
  } 
  edge_global_coord <- patch_id[edge_local_coord]
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


remove_extra_nodes <- function(cutset_local, patch_id){
  irm <- apply(cutset_local[, 1:2], 1, function(x) any(x > length(patch_id)))
  return(cutset_local[!irm, ])
}

update_cutset_global <- function(arcs, cutset_local, cutset_global, patch_id) {
  seam_nodes <- sort(unique(unlist(arcs[, 1:2])))
  seam_nodes <- seam_nodes[seam_nodes > (length(patch_id) + 2)] 
  for (node in seam_nodes){
    cutset_global <- update_cutset_oneseamnode(
      seam_node = node,
      arcs = arcs,
      cutset_local = cutset_local,
      cutset_global = cutset_global,
      patch_id = patch_id
    )
  }
  cutset_global <- rbind(
    cutset_global,
    local_to_global_arcs(remove_extra_nodes(cutset_local, patch_id), patch_id)
  )
  return(cutset_global)
}
