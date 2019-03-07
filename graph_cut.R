library(magick)
library(magrittr)
library(igraph)
library(optrees)

source("graph_cut_algo.R")

img <- image_read("olives.gif")
plot(img)
img  %<>% 
  image_convert(type = 'grayscale') %>%
  image_data("gray") %>%
  "["(1,,) %>% 
  as.integer() %>%
  matrix(ncol = 128, nrow = 128)
image(img, col = grey.colors(256))

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

arcs <- create_arcs(diff_patch = diff_AB)
graph <- graph_from_data_frame(arcs, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Add Source and Sink#
arcs_with_ST <- add_sink_and_source(arcs, patch_A_id, add_source = add_source_to_left)

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

# expressing the seam cut with respect to the global node names
cutset_local <- rearrange_node(cutset_local)


cutset_global <- local_to_global_arcs(cutset_local, patch_A_id)
cutset_global <- rearrange_node(cutset_global)

# getting the combined patch
patched <- patch_A
patched_origin <- matrix(1, nrow = nrow(patch_A), ncol = ncol(patch_A))
tcut <- mincut$t.cut[mincut$t.cut <= length(patch_A)]
patched[tcut] <- patch_B[tcut]
patched_origin[tcut] <- 2


patched <- assemble_patch(mincut = mincut, patch_A = patch_A, patch_B = patch_B, patch_id = patch_A_id)

# opposite patch
opposite_patched <- assemble_opposite_patch(mincut = mincut, patch_A = patch_A, patch_B = patch_B, patch_id = patch_A_id)

  

canvas <- update_canvas(canvas = canvas, img_B = img_B, img_B_id = img_B_id, patched = patched, patched_id = patch_A_id) 
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
arcs_newpatch <- create_arcs(diff_patch = diff_ABC)

# Add Source and Sink#
arcs_newpatch  <- add_sink_and_source(arcs_newpatch, patch_AB_id)

arcs_newpatch <- reorder_nodes(arcs_newpatch)
cutset_global <- reorder_nodes(cutset_global)


par(mfrow = c(1, 1))
graph <- graph_from_data_frame(arcs_newpatch, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Check if old seam in patch

cutset_global 

cut_1 <- cutset_global[1, ]

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
