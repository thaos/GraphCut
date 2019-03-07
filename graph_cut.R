library(magick)
library(magrittr)
library(igraph)
library(optrees)

img <- image_read("olives.gif")
plot(img)
img  %<>% 
  image_convert(type = 'grayscale') %>%
  image_data("gray") %>%
  "["(1,,) %>% 
  as.integer() %>%
  matrix(ncol = 128, nrow = 128)
image(img, col = grey.colors(256))

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

add_source_to_left <- function(diff_patch){
  idx_source <- 1 + (seq.int(ncol(diff_patch)) - 1) * nrow(diff_patch)
  idx_sink <-  seq.int(ncol(diff_patch)) * nrow(diff_patch) 
  arcs <- rbind( 
    data.frame(from = length(diff_patch) + 1, to = idx_source, weight = 10E6),
    data.frame(from = length(diff_patch) + 2, to = idx_sink, weight = 10E6)
  )
  return(arcs)
}

canvas <- canvas_origin <- opposite_canvas <- matrix(NA, ncol = 32, nrow = 32*2 - 8)
canvas_id <- matrix(1:length(canvas), ncol = ncol(canvas), nrow = nrow(canvas))
image(canvas_id)

img_A <- img[1:32, 1:32]
canvas[1:32, 1:32] <- img_A
canvas_origin[1:32, 1:32] <- 1
image(canvas, col = grey.colors(256), zlim = c(0, 256))

img_A_id <- canvas_id[1:32, 1:32]
patch_A <- canvas[25:32, 1:32]
patch_A_id <- canvas_id[25:32, 1:32]
patch_A_id_local <- matrix(1:length(patch_A), ncol = ncol(patch_A), nrow = nrow(patch_A))
plot(patch_A, patch_A_id_local)

# search for image B by scanning
ij_scan <- expand.grid(i = 1:(128 - 31), j = 1:(128 - 31))
ij_scan <- subset(ij_scan, !(i <= 32 & j <= 32)) 
plot(ij_scan)
compute_dist <- function(i, j){
  sum(abs(img_A[25:32, ] - img[i:(i+31), j:(j+31)][1:8,]))
}
dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j) 
imin <- which.min(dmat)
imin <- which(order(dmat) == 2)

img_B <- with(ij_scan[imin, ], img[i:(i+31), j:(j+31)])
img_B_id <- canvas_id[25:(25 + 31), 1:32]
patch_B <- img_B[1:8, 1:32]
diff_AB <- abs(patch_A - patch_B)

# create patch graph
# alway the same for patch of same size
# need only to be done once
arcs <- lapply(
  1:length(diff_AB),
  create_egdes,
  diff_patch = diff_AB
) %>%
  do.call(rbind, .) %>%
  subset(from < to)
graph <- graph_from_data_frame(arcs, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Add Source and Sink#
arcs_with_ST <- rbind(
  arcs,
  add_source_to_left(diff_AB)
)

graph <- graph_from_data_frame(arcs_with_ST, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Solve MaxFlow / ST-MinCut
mincut <- findMinCut(
  nodes = seq.int(length(diff_AB) + 2) ,
  arcs = as.matrix(arcs_with_ST),
  algorithm = "Ford-Fulkerson",
  source.node = length(diff_AB) + 1,
  sink.node = length(diff_AB) + 2,
  directed = FALSE
)

mincut

# keeping the seam cut cost
cutset_local <- mincut$cut.set

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
# expressing the seam cut with respect to the global node names
cutset_local <- rearrange_node(cutset_local)

local_to_global_arcs <- function(arcs_local, patch_id){
  arcs_global <- arcs_local
  arcs_global[, 1:2] <- patch_id[c(arcs_local[, 1:2])]
  return(arcs_global)
}
cutset_global <- local_to_global_arcs(cutset_local, patch_A_id)
cutset_global <- rearrange_node(cutset_global)

# getting the combined patch
patched <- patch_A
patched_origin <- matrix(1, nrow = nrow(patch_A), ncol = ncol(patch_A))
tcut <- mincut$t.cut[mincut$t.cut <= length(patch_A)]
patched[tcut] <- patch_B[tcut]
patched_origin[tcut] <- 2

assemble_patch <- function(mincut, patch_A, patch_B, patch_id){
  patched <- patch_A
  tcut <- mincut$t.cut[mincut$t.cut <= length(patch_id)]
  patched[tcut] <- patch_B[tcut]
  return(patched)
}
patched <- assemble_patch(mincut = mincut, patch_A = patch_A, patch_B = patch_B, patch_id = patch_A_id)

# opposite patch
patched_alt <- patch_B
patched_alt[tcut] <- patch_A[tcut]


assemble_opposite_patch <- function(mincut, patch_A, patch_B, patch_id){
  opposite_patched <- patch_B
  tcut <- mincut$t.cut[mincut$t.cut <= length(patch_id)]
  opposite_patched[tcut] <- patch_A[tcut]
  return(opposite_patched)
}
opposite_patched <- assemble_opposite_patch(mincut = mincut, patch_A = patch_A, patch_B = patch_B, patch_id = patch_A_id)

  
# adding patches to the canvas
update_canvas <- function(canvas, img_B, img_B_id, patched, patched_id){
  canvas[img_B_id] <- img_B
  canvas[patched_id] <- patched
  return(canvas)
}
canvas <- update_canvas(canvas = canvas, img_B = img_B, img_B_id = img_B_id, patched = patched, patched_id = patch_A_id) 

update_opposite_canvas <- function(opposite_canvas, opposite_patched, patched_id){
  opposite_canvas[patched_id] <- opposite_patched
  return(opposite_canvas)
}
opposite_canvas <- update_opposite_canvas(opposite_canvas = opposite_canvas, opposite_patched = opposite_patched, patched_id = patch_A_id)


canvas_origin[img_B_id] <- 2
canvas_origin[patch_A_id] <- patched_origin
par(mfrow = c(1, 1))
image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = canvas_origin
)

par(mfrow = c(1, 1))
image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = canvas,
  col = grey.colors(256), zlim = c(0, 256)
)

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
lines_seams(cutset_global = cutset_global, canvas = canvas)

# Repass with patch C
patch_AB <- canvas[25:32, 1:32]
patch_AB_id <- canvas_id[25:32, 1:32]
image(25:32, 1:32, patch_AB, col = grey.colors(256), zlim = c(0, 256))

# search for image C by scanning
ij_scan <- expand.grid(i = 1:(128 - 31), j = 1:(128 - 31))
ij_scan <- subset(ij_scan, !(i <= 32 & j <= 32)) 
plot(ij_scan)
compute_dist <- function(i, j){
  sum(abs(patch_AB - img[i:(i+31), j:(j+31)][1:8,]))
}
dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j) 
imin <- which.min(dmat)
# imin <- which(order(dmat) == 3)

img_C <- with(ij_scan[imin, ], img[i:(i+31), j:(j+31)])
img_C_id <- canvas_id[25:(25 + 31), 1:32]
patch_C <- img_C[1:8, 1:32]
diff_ABC <- abs(patch_AB - patch_C)
par(mfrow = c(3, 1))
image(patch_AB, col = rev(grey.colors(256)), zlim = c(0, 256))
image(patch_C, col = rev(grey.colors(256)), zlim = c(0, 256))
image(abs(diff_ABC), col = rev(grey.colors(256)), zlim = c(0, 256))
image(abs(patch_B - patch_C), col = rev(grey.colors(256)), zlim = c(0, 256))

# generate graph data
arcs_newpatch <- lapply(
  1:length(diff_ABC),
  create_egdes,
  diff_patch = diff_ABC
) %>%
  do.call(rbind, .) %>%
  subset(from < to)

# Add Source and Sink#
arcs_newpatch  <- rbind(
  arcs_newpatch ,
  add_source_to_left(diff_ABC)
)

reorder_nodes <- function(arcs){
  arcs <- arcs[order(arcs[, 1], arcs[, 2]), ]
}
arcs_newpatch <- reorder_nodes(arcs_newpatch)
cutset_global <- reorder_nodes(cutset_global)


par(mfrow = c(1, 1))
graph <- graph_from_data_frame(arcs_newpatch, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Check if old seam in patch

cutset_global 

cut_1 <- cutset_global[1, ]
check_oldseam_in_patch <- function(cut, patch_id){
  (cut[1] %in% patch_id & cut[2] %in% patch_id)
}
check_oldseam_in_patch(cut = cut_1, patch_id = patch_AB_id)


compute_oldseam_weight1 <- function(cut, patch, patch_id, canvas, opposite_canvas){
  abs(canvas[cut[1]] - patch[patch_id == cut[1]]) + abs(opposite_canvas[cut[2]] - patch[patch_id == cut[2]])
}

compute_oldseam_weight2 <- function(cut, patch, patch_id, canvas, opposite_canvas){
  abs(canvas[cut[2]] - patch[patch_id == cut[2]]) + abs(opposite_canvas[cut[1]] - patch[patch_id == cut[1]])
}

compute_oldseam_weight1(cut = cut_1, patch = patch_C, patch_id = patch_A_id, canvas = canvas, opposite_canvas = opposite_canvas)
compute_oldseam_weight2(cut = cut_1, patch = patch_C, patch_id = patch_A_id, canvas = canvas, opposite_canvas = opposite_canvas)

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
add_oldcut(
  cut = cut_1,
  patch = patch_C, patch_id = patch_AB_id,
  canvas = canvas, opposite_canvas = opposite_canvas,
  arcs = arcs_newpatch
) %>% 
  subset(from == 259 | to == 259)

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

arcs_newpatch_updated <- update_graph_with_cuts(
  cutset_global = cutset_global,
  patch = patch_C, patch_id = patch_AB_id,
  canvas = canvas, opposite_canvas = opposite_canvas,
  arcs = arcs_newpatch
)
arcs_newpatch_updated <- reorder_nodes(arcs_newpatch_updated)

par(mfrow = c(1, 1))
graph <- graph_from_data_frame(arcs_newpatch_updated, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))
 
# Solve MaxFlow / ST-MinCut
mincut <- findMinCut(
  nodes = sort(unique(unlist(arcs_newpatch_updated[, 1:2]))) ,
  arcs = as.matrix(arcs_newpatch_updated),
  algorithm = "Ford-Fulkerson",
  source.node = length(patch_AB) + 1,
  sink.node = length(patch_AB) + 2,
  directed = FALSE
)

mincut
cutset_local <- mincut$cut.set
cutset_local <- rearrange_node(cutset_local)
cutset_local <- reorder_nodes(cutset_local)
head(arcs_newpatch_updated)


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
update_cutset_oneseamnode(
  seam_node = 259,
  arcs = arcs_newpatch_updated,
  cutset_local = cutset_local,
  cutset_global = cutset_global,
  patch_id = patch_AB_id
)

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

patched <- assemble_patch(mincut = mincut, patch_A = patch_AB, patch_B = patch_C, patch_id = patch_AB_id)
par(mfrow = c(3, 1))
image(patch_AB, col = grey.colors(256), zlim = c(0, 256))
image(patch_C, col = grey.colors(256), zlim = c(0, 256))
image(patched, col = grey.colors(256), zlim = c(0, 256))


par(mfrow = c(2, 1))
image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = canvas,
  col = grey.colors(256), zlim = c(0, 256)
)
lines_seams(cutset_global = cutset_global, canvas = canvas)

# debug(update_cutset_oneseamnode)
cutset_global <- update_cutset_global(arcs_newpatch_updated, cutset_local = cutset_local, cutset_global =  cutset_global, patch_id = patch_AB_id)


opposite_patched <- assemble_opposite_patch(mincut = mincut, patch_A = patch_AB, patch_B = patch_C, patch_id = patch_AB_id)
canvas <- update_canvas(canvas = canvas, img_B = img_C, img_B_id = img_C_id, patched = patched, patched_id = patch_AB_id) 
opposite_canvas <- update_opposite_canvas(opposite_canvas = opposite_canvas, opposite_patched = opposite_patched, patched_id = patch_AB_id)
  

image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = canvas,
  col = grey.colors(256), zlim = c(0, 256)
)
lines_seams(cutset_global = cutset_global, canvas = canvas)

image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = opposite_canvas,
  col = grey.colors(256), zlim = c(0, 256)
)
lines_seams(cutset_global = cutset_global, canvas = canvas)
