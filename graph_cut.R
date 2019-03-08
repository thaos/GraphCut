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



patch_A <- img[1:32, 1:32]
patch_A_id <- canvas_id[1:32, 1:32]
canvas <- update_canvas(canvas = canvas, patch = patch_A, patch_id = patch_A_id)
canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "1", patch_id = patch_A_id)
overlap_A <- canvas[25:32, 1:32]
overlap_A_id <- canvas_id[25:32, 1:32]
overlap_A_id_local <- matrix(1:length(overlap_A), ncol = ncol(overlap_A), nrow = nrow(overlap_A))
overlap_list <- data.frame(number = "1", overlap = I(list(overlap_A)), overlap_id = I(list(overlap_A_id)))

plot(overlap_A, overlap_A_id_local)
image(canvas, col = grey.colors(256), zlim = c(0, 256))

# search for image B by scanning
ij_scan <- expand.grid(i = 1:(128 - 31), j = 1:(128 - 31))
ij_scan <- subset(ij_scan, !(i <= 32 & j <= 32)) 
plot(ij_scan)
compute_dist <- function(i, j){
  sum(abs(patch_A[25:32, ] - img[i:(i+31), j:(j+31)][1:8,]))
}
dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j) 
imin <- which.min(dmat)
imin <- which(order(dmat) == 2)

patch_B <- with(ij_scan[imin, ], img[i:(i+31), j:(j+31)])
patch_B_id <- canvas_id[25:(25 + 31), 1:32]
overlap_B <- patch_B[1:8, 1:32]
diff_AB <- abs(overlap_A - overlap_B)


overlap_list <- add_overlap_to_list(overlap_list = overlap_list, number = 2, overlap = overlap_B, overlap_id = overlap_A_id)

# create overlap graph
# alway the same for overlap of same size
# need only to be done once

arcs <- create_arcs(diff_overlap = diff_AB)
graph <- graph_from_data_frame(arcs, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Add Source and Sink#
arcs_with_ST <- add_sink_and_source(arcs, overlap_A_id, add_source = add_source_to_left)

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


cutset_global <- local_to_global_arcs(cutset_local, overlap_A_id)
cutset_global <- rearrange_node(cutset_global)

# getting the combined overlap
overlaped <- assemble_overlap(mincut = mincut, overlap_A = overlap_A, overlap_B = overlap_B, overlap_id = overlap_A_id)


overlap_origin <- overlap_A_id
overlap_origin[] <- "A"
overlaped_origin <- assemble_overlap_origin(mincut = mincut, label = "2", overlap_origin = overlap_origin, overlap_id = overlap_A_id)  
image(matrix(as.numeric(factor(overlaped_origin)), ncol = ncol(canvas)))

canvas <- update_canvas(canvas = canvas, patch = patch_B, patch_id = patch_B_id, overlaped = overlaped, overlaped_id = overlap_A_id) 


canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "2", patch_id = patch_B_id, overlaped_origin = overlaped_origin, overlaped_id = overlap_A_id)

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

# Repass with overlap C
overlap_AB <- canvas[25:32, 1:32]
overlap_AB_id <- canvas_id[25:32, 1:32]
image(25:32, 1:32, overlap_AB, col = grey.colors(256), zlim = c(0, 256))

# search for image C by scanning
ij_scan <- expand.grid(i = 1:(128 - 31), j = 1:(128 - 31))
ij_scan <- subset(ij_scan, !(i <= 32 & j <= 32)) 
plot(ij_scan)
compute_dist <- function(i, j){
  sum(abs(overlap_AB - img[i:(i+31), j:(j+31)][1:8,]))
}
dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j) 
imin <- which.min(dmat)
# imin <- which(order(dmat) == 3)

patch_C <- with(ij_scan[imin, ], img[i:(i+31), j:(j+31)])
patch_C_id <- canvas_id[25:(25 + 31), 1:32]
overlap_C <- patch_C[1:8, 1:32]

overlap_list <- add_overlap_to_list(overlap_list = overlap_list, number = 3, overlap = overlap_C, overlap_id = overlap_A_id)

diff_ABC <- abs(overlap_AB - overlap_C)
par(mfrow = c(3, 1))
image(overlap_AB, col = rev(grey.colors(256)), zlim = c(0, 256))
image(overlap_C, col = rev(grey.colors(256)), zlim = c(0, 256))
image(abs(diff_ABC), col = rev(grey.colors(256)), zlim = c(0, 256))
image(abs(overlap_B - overlap_C), col = rev(grey.colors(256)), zlim = c(0, 256))



# generate graph data
arcs_newoverlap <- create_arcs(diff_overlap = diff_ABC)

# Add Source and Sink#
arcs_newoverlap  <- add_sink_and_source(arcs_newoverlap, overlap_AB_id)

arcs_newoverlap <- reorder_nodes(arcs_newoverlap)
cutset_global <- reorder_nodes(cutset_global)


par(mfrow = c(1, 1))
graph <- graph_from_data_frame(arcs_newoverlap, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Check if old seam in overlap

cutset_global 

cut_1 <- cutset_global[1, ]

arcs_newoverlap_updated <- update_graph_with_cuts(
  cutset_global = cutset_global,
  overlap = overlap_C, overlap_id = overlap_AB_id,
  canvas = canvas, canvas_origin =  canvas_origin,
  overlap_list = overlap_list, arcs = arcs_newoverlap
)
arcs_newoverlap_updated <- reorder_nodes(arcs_newoverlap_updated)

par(mfrow = c(1, 1))
graph <- graph_from_data_frame(arcs_newoverlap_updated, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))
 
# Solve MaxFlow / ST-MinCut
mincut <- findMinCut(
  nodes = sort(unique(unlist(arcs_newoverlap_updated[, 1:2]))) ,
  arcs = as.matrix(arcs_newoverlap_updated),
  algorithm = "Ford-Fulkerson",
  source.node = length(overlap_AB) + 1,
  sink.node = length(overlap_AB) + 2,
  directed = FALSE
)

mincut
cutset_local <- mincut$cut.set
cutset_local <- rearrange_node(cutset_local)
cutset_local <- reorder_nodes(cutset_local)
head(arcs_newoverlap_updated)



overlaped <- assemble_overlap(mincut = mincut, overlap_A = overlap_AB, overlap_B = overlap_C, overlap_id = overlap_AB_id)
overlaped_origin <- assemble_overlap_origin(mincut = mincut, label = "C", overlap_origin = overlaped_origin, overlap_id = overlap_A_id)  
par(mfrow = c(2, 2))
image(overlap_AB, col = grey.colors(256), zlim = c(0, 256))
image(overlap_C, col = grey.colors(256), zlim = c(0, 256))
image(overlaped, col = grey.colors(256), zlim = c(0, 256))
image(matrix(as.numeric(factor(overlaped_origin)), ncol = ncol(canvas)))

par(mfrow = c(2, 1))
image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = canvas,
  col = grey.colors(256), zlim = c(0, 256)
)
lines_seams(cutset_global = cutset_global, canvas = canvas)

# debug(update_cutset_oneseamnode)
cutset_global <- update_cutset_global(arcs_newoverlap_updated, cutset_local = cutset_local, cutset_global =  cutset_global, overlap_id = overlap_AB_id)


canvas <- update_canvas(canvas = canvas, patch = patch_C, patch_id = patch_C_id, overlaped = overlaped, overlaped_id = overlap_AB_id) 
canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "3", patch_id = patch_C_id, overlaped_origin = overlaped_origin, overlaped_id = overlap_A_id)  

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
  stopifnot(length(source_id) > 0)
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


add_newpatch <- function(
  xstart, ystart, xlength = 32, ylength = 32 ,
  canvas, canvas_origin, canvas_id,
  training_img,
  overlap_list = NULL,  cutset_global = NULL
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
  if(is.null(overlap_list)){
    overlap_list <- data.frame(
      number = 1,
      overlap = I(list(overlap_old)),
      overlap_id = I(list(overlap_id))
    )
  }
  overlap_list <- add_overlap_to_list(
    overlap_list = overlap_list, number = max(overlap_list$number) + 1,
    overlap = overlap_new, overlap_id = overlap_id
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
      overlap_list = overlap_list, arcs = arcs
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
    mincut = mincut, label = max(overlap_list$number) + 1,
    overlap_origin = overlaped_origin, overlap_id = overlap_id
  )  
  canvas <- update_canvas(
    canvas = canvas,
    patch = patch_new, patch_id = patch_id,
    overlaped = overlaped, overlaped_id = overlap_id
  ) 
  canvas_origin <- update_canvas_origin(
    canvas_origin = canvas_origin,
    label = max(overlap_list$number) + 1,
    patch_id = patch_id,
    overlaped_origin = overlaped_origin,
    overlaped_id = overlap_id
  )  
  list(
    canvas = canvas, canvas_origin = canvas_origin, canvas_id = canvas_id,
    overlap_id = overlap_id, overlap_list = overlap_list,
    cutset_global = cutset_global
  )
}
canvas <- canvas_origin <- opposite_canvas <- matrix(NA, ncol = 32, nrow = 32*2 - 8)
canvas_id <- matrix(1:length(canvas), ncol = ncol(canvas), nrow = nrow(canvas))
patch_A <- img[1:32, 1:32]
patch_A_id <- canvas_id[1:32, 1:32]
canvas <- update_canvas(canvas = canvas, patch = patch_A, patch_id = patch_A_id)
canvas_origin <- update_canvas_origin(canvas_origin = canvas_origin, label = "1", patch_id = patch_A_id)  
image(seq.int(nrow(canvas)), seq.int(ncol(canvas)), canvas, , zlim = c(0, 256), col = grey.colors(256))


new_canvas <- add_newpatch(
  xstart = 20, ystart = 15,
  xlength = 20, ylength = 10,
  canvas = canvas, canvas_origin =  canvas_origin, canvas_id = canvas_id,
  training_img = img, overlap_list = NULL, cutset_global = NULL
)
image(1:nrow(canvas), 1:ncol(canvas), new_canvas$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas$cutset_global, canvas = new_canvas$canvas)

new_canvas2 <- add_newpatch(
  xstart = 20, ystart = 23,
  xlength = 20, ylength = 10,
  canvas = new_canvas$canvas, canvas_origin =  new_canvas$canvas_origin, canvas_id = new_canvas$canvas_id,
  training_img = img, overlap_list = new_canvas$overlap_list, cutset_global = new_canvas$cutset_global
)
image(1:nrow(canvas), 1:ncol(canvas), new_canvas2$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas2$cutset_global, canvas = new_canvas2$canvas)
abline(v = 19.5)
abline(h = 22.5)

new_canvas3 <- add_newpatch(
  xstart = 20, ystart = 14,
  xlength = 20, ylength = 10,
  canvas = new_canvas2$canvas, canvas_origin =  new_canvas2$canvas_origin, canvas_id = new_canvas2$canvas_id,
  training_img = img, overlap_list = new_canvas2$overlap_list, cutset_global = new_canvas2$cutset_global
)
image(1:nrow(canvas), 1:ncol(canvas), new_canvas3$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas3$cutset_global, canvas = new_canvas3$canvas)
abline(v = 19.5)
abline(h = 22.5)
      
new_canvas4 <- add_newpatch(
  xstart = 20, ystart = 4,
  xlength = 20, ylength = 10,
  canvas = new_canvas3$canvas, canvas_origin =  new_canvas3$canvas_origin, canvas_id = new_canvas3$canvas_id,
  training_img = img, overlap_list = new_canvas3$overlap_list, cutset_global = new_canvas3$cutset_global
)
image(1:nrow(canvas), 1:ncol(canvas), new_canvas4$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas4$cutset_global, canvas = new_canvas4$canvas)
abline(v = 19.5)
abline(h = 22.5)

new_canvas5 <- add_newpatch(
  xstart = 20, ystart = 1,
  xlength = 20, ylength = 10,
  canvas = new_canvas4$canvas, canvas_origin =  new_canvas4$canvas_origin, canvas_id = new_canvas4$canvas_id,
  training_img = img, overlap_list = new_canvas4$overlap_list, cutset_global = new_canvas4$cutset_global
)
image(1:nrow(canvas), 1:ncol(canvas), new_canvas5$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas5$cutset_global, canvas = new_canvas5$canvas)
abline(v = 19.5)
abline(h = 22.5)

new_canvas6 <- add_newpatch(
  xstart = 30, ystart = 1,
  xlength = 27, ylength = 32,
  canvas = new_canvas5$canvas, canvas_origin =  new_canvas5$canvas_origin, canvas_id = new_canvas5$canvas_id,
  training_img = img, overlap_list = new_canvas5$overlap_list, cutset_global = new_canvas5$cutset_global
)
image(1:nrow(canvas), 1:ncol(canvas), new_canvas6$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas6$cutset_global, canvas = new_canvas6$canvas)
abline(v = 19.5)
abline(h = 22.5)
      