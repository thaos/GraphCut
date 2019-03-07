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

patch_list <- data.frame(name = "A", patch = I(list(patch_A)), patch_id = I(list(patch_A_id)))

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

add_patch_to_list <- function(patch_list, name, patch, patch_id){
  rbind(
    patch_list,
    data.frame(name = name, patch = I(list(patch)), patch_id = I(list(patch_id)))
  )
}
patch_list <- add_patch_to_list(patch_list = patch_list, name = "B", patch = patch_B, patch_id = patch_A_id)

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
patched <- assemble_patch(mincut = mincut, patch_A = patch_A, patch_B = patch_B, patch_id = patch_A_id)

assemble_patch_origin <- function(mincut, label_B,  patch_origin, patch_id){
  patched_origin <- patch_origin
  tcut <- mincut$t.cut[mincut$t.cut <= length(patch_id)]
  patched_origin[tcut] <- label_B
  return(patched_origin)
}

patch_origin <- patch_A_id
patch_origin[] <- "A"
patched_origin <- assemble_patch_origin(mincut = mincut, label_B = "B", patch_origin = patch_origin, patch_id = patch_A_id)  
image(matrix(as.numeric(factor(patched_origin)), ncol = ncol(canvas)))

canvas <- update_canvas(canvas = canvas, img_B = img_B, img_B_id = img_B_id, patched = patched, patched_id = patch_A_id) 

# adding patches to the canvas
update_canvas_origin <- function(canvas_origin, label, img_id, patched_origin = NULL, patched_id = NULL){
  canvas_origin[img_id] <- label
  if(!is.null(patched_origin) & !is.null(patched_id)){
    canvas_origin[patched_id] <- patched_origin
  }
  return(canvas_origin)
}

canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "A", img_id = img_A_id)
canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "B", img_id = img_B_id, patched_origin = patched_origin, patched_id = patch_A_id)

par(mfrow = c(1, 1))
image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = matrix(as.numeric(factor(canvas_origin)), ncol = ncol(canvas))
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

patch_list <- add_patch_to_list(patch_list = patch_list, name = "C", patch = patch_C, patch_id = patch_A_id)

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
  canvas = canvas, canvas_origin =  canvas_origin,
  patch_list = patch_list, arcs = arcs_newpatch
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
patched_origin <- assemble_patch_origin(mincut = mincut, label_B = "C", patch_origin = patched_origin, patch_id = patch_A_id)  
par(mfrow = c(2, 2))
image(patch_AB, col = grey.colors(256), zlim = c(0, 256))
image(patch_C, col = grey.colors(256), zlim = c(0, 256))
image(patched, col = grey.colors(256), zlim = c(0, 256))
image(matrix(as.numeric(factor(patched_origin)), ncol = ncol(canvas)))

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


canvas <- update_canvas(canvas = canvas, img_B = img_C, img_B_id = img_C_id, patched = patched, patched_id = patch_AB_id) 
canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "C", img_id = img_C_id, patched_origin = patched_origin, patched_id = patch_A_id)  

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
  z = matrix(as.numeric(factor(canvas_origin)), ncol = ncol(canvas)),
)
lines_seams(cutset_global = cutset_global, canvas = canvas)

