library(magick)
library(magrittr)
library(igraph)
img <- image_read("olives.gif")
plot(img)
img  %<>% 
  image_convert(type = 'grayscale') %>%
  image_data("gray") %>%
  "["(1,,) %>% 
  as.integer() %>%
  matrix(ncol = 128, nrow = 128)
image(img, col = grey.colors(256))

img_A <- img[1:32, 1:32]
image(img_A, col = grey.colors(256))
image(img_A[25:32, ], col = grey.colors(256))

ij_scan <- expand.grid(i = 1:(128 - 31), j = 1:(128 - 31))
ij_scan <- subset(ij_scan, !(i <= 32 & j <= 32)) 
plot(ij_scan)
compute_dist <- function(i, j){
  sum(abs(img_A[25:32, ] - img[i:(i+31), j:(j+31)][1:8,]))
}
dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j) 
imin <- which.min(dmat)

img_B <- with(ij_scan[imin, ], img[i:(i+31), j:(j+31)])
image(img_B, col = grey.colors(256))

image(abs(img_A[25:32, ] - img_B[1:8, ]), col = rev(grey.colors(256)))
sum(abs(img_A[25:32, ] - img_B[1:8, ]))


par(mfrow = c(2, 1))
image(img_A, col = grey.colors(256), zlim = c(0, 256))
image(img_B, col = grey.colors(256), zlim = c(0, 256))
image(img_A[25:32, ], col = grey.colors(256), zlim = c(0, 256))
image(img_B[1:8, ], col = grey.colors(256))

patch_A <- img_A[25:32, ]
patch_B <- img_B[1:8, ]

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
  adj <- ans$i + ni * (ans$j - 1)
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
diff_AB <- abs(patch_A - patch_B)
create_egdes(1, diff_patch = diff_AB)            
create_egdes(8, diff_patch = diff_AB)            
  
graph_df <- lapply(
  1:length(diff_AB),
  create_egdes,
  diff_patch = diff_AB
) %>%
  do.call(rbind, .) %>%
  subset(from < to)
graph <- graph_from_data_frame(graph_df, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

add_source_to_left <- function(diff_patch){
  idx_source <- 1 + (seq.int(ncol(diff_patch)) - 1) * nrow(diff_patch)
  idx_sink <-  seq.int(ncol(diff_patch)) * nrow(diff_patch) 
  graph_df <- rbind( 
    data.frame(from = "S", to = idx_source, weight = 10E6),
    data.frame(from = "T", to = idx_sink, weight = 10E6)
  )
  return(graph_df)
}

graph_with_ST_df <- rbind(
  graph_df,
  add_source_to_left(diff_AB)
)
graph_with_ST <- graph_from_data_frame(graph_with_ST, directed = FALSE)
plot(graph_with_ST)
edge_attr(graph)

min_cut(graph = graph_with_ST, source = "S", target = "T", capacity = E(graph_with_ST)$weight) 

graph_with_ST_df$from <- ifelse(graph_with_ST_df$from == "S", length(diff_AB) + 1, graph_with_ST_df$from)
graph_with_ST_df$from <- ifelse(graph_with_ST_df$from == "T", length(diff_AB) + 2, graph_with_ST_df$from)
graph_with_ST_df$from <- as.integer(graph_with_ST_df$from)

library(optrees)
mincut <- findMinCut(
  nodes = seq.int(length(diff_AB) + 2) ,
  arcs = as.matrix(graph_with_ST_df),
  algorithm = "Ford-Fulkerson",
  source.node = length(diff_AB) + 1,
  sink.node = length(diff_AB) + 2,
  directed = FALSE
)

par(mfrow = c(1, 1))
patched <- patch_A
patched[-mincut$t.cut] <- 0 
patched[mincut$t.cut] <- 1
patched <- matrix(patched, ncol = ncol(patch_A), nrow = nrow(patch_A))
image(patched)

par(mfrow = c(3, 1))
patched <- patch_A
tcut <- mincut$t.cut[mincut$t.cut <= length(patch_A)]
patched[tcut] <- patch_B[tcut]
image(patched, col = grey.colors(256), zlim = c(0, 256))
image(patch_A, col = grey.colors(256), zlim = c(0, 256))
image(patch_B, col = grey.colors(256), zlim = c(0, 256))

par(mfrow = c(3, 1))
combined <- rbind(img_A[-(25:32), ], patched, img_B[ -(1:8), ])
image(combined, col = grey.colors(256), zlim = c(0, 256))
image(rbind(img_A,  img_B[ -(1:8), ]), col = grey.colors(256), zlim = c(0, 256))
image(rbind(img_A[-(25:32), ],  img_B), col = grey.colors(256), zlim = c(0, 256))

cutset <- mincut$cut.set
cutset <- apply(
  cutset, 1, 
  function(x){
    ifelse(x[1] > x[2], x[c(2, 1, 3)], return(x))
  }
) %>%
  t()
colnames(cutset) <- colnames(mincut$cut.set)


canvas <- matrix(NA, ncol = 32, nrow = 32*2 - 8)
canvas_id <- matrix(1:length(canvas), ncol = ncol(canvas), nrow = nrow(canvas))
image(canvas_id)

canvas[1:32, 1:32] <- img_A
image(canvas, col = grey.colors(256), zlim = c(0, 256))

img_A <- canvas[1:32, 1:32]
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

img_B <- with(ij_scan[imin, ], img[i:(i+31), j:(j+31)])
img_B_id <- canvas_id[25:(25 + 31), 1:32]
patch_B <- img_B[1:8, 1:32]
diff_AB <- abs(patch_A - patch_B)

# create patch graph
# alway the same for patch of same size
# need only to be done once
graph_df <- lapply(
  1:length(diff_AB),
  create_egdes,
  diff_patch = diff_AB
) %>%
  do.call(rbind, .) %>%
  subset(from < to)
graph <- graph_from_data_frame(graph_df, directed = FALSE)
plot(graph, edge.width = exp(E(graph)$weight/max(E(graph)$weight)))

# Add Source and Sink#
graph_with_ST_df <- rbind(
  graph_df,
  add_source_to_left(diff_AB)
)
# Replacing S and T par node number (the two last nodes)
graph_with_ST_df$from <- ifelse(graph_with_ST_df$from == "S", length(diff_AB) + 1, graph_with_ST_df$from)
graph_with_ST_df$from <- ifelse(graph_with_ST_df$from == "T", length(diff_AB) + 2, graph_with_ST_df$from)
graph_with_ST_df$from <- as.integer(graph_with_ST_df$from)

# keeping the seam cut cost
cutset <- mincut$cut.set
cutset <- apply(
  cutset, 1, 
  function(x){
    ifelse(x[1] > x[2], x[c(2, 1, 3)], return(x))
  }
) %>%
  t()
colnames(cutset) <- colnames(mincut$cut.set)
# expressing the seam cut with respect to the global node names
cutset[, 1:2] <- patch_A_id[c(cutset[, 1:2])]

# getting the combined patch
patched <- patch_A
tcut <- mincut$t.cut[mincut$t.cut <= length(patch_A)]
patched[tcut] <- patch_B[tcut]

# adding patches to the canvas
canvas[img_B_id] <- img_B
canvas[patch_A_id] <- patched
par(mfrow = c(1, 1))
image(
  x = seq.int(nrow(canvas)),
  y = seq.int(ncol(canvas)),
  z = canvas,
  col = grey.colors(256), zlim = c(0, 256)
)

n1_coord <- arrayInd(cutset[, 1], .dim = dim(canvas))
n2_coord <- arrayInd(cutset[, 2], .dim = dim(canvas))
ncolors = 10
wcuts <- cut(cutset[, 3], breaks = quantile(cutset[,3], probs = seq(0, 1, length.out = ncolors + 1))) 
wcuts <- cut(cutset[, 3], breaks = 10) 
pts <- matrix(NA, ncol = 2, nrow = 2)
for(i in seq.int(nrow(cutset))){
  if(n1_coord[i, 1]  == n2_coord[i, 1]){
    pts[, 2] = mean( c(n1_coord[i, 2], n2_coord[i, 2]))
    pts[, 1] = n1_coord[i, 1] + c(-0.5, 0.5)
  } else {
    pts[, 1] = mean( c(n1_coord[i, 1], n2_coord[i, 1]))
    pts[, 2] = n1_coord[i, 2] + c(-0.5, 0.5)
  }
  lines(pts, col = rev(heat.colors(ncolors))[wcuts[i]], lwd = 2)
}

newpatch_center <- 644
newpatch_center_coordinate <- arrayInd(newpatch_center, .dim = dim(canvas))
xnewpatch <- seq(newpatch_center_coordinate[1] - 5, newpatch_center_coordinate[1] + 5)
ynewpatch <- seq(newpatch_center_coordinate[2] - 5, newpatch_center_coordinate[2] + 5)
newpatch <- canvas[ xnewpatch, ynewpatch]
newpatch_id <- canvas_id[ xnewpatch, ynewpatch]
image(xnewpatch, ynewpatch, newpatch, col = grey.colors(256), zlim = c(0, 256))

# find for closest patch in training image
ij_scan <- expand.grid(i = 1:(128 - 31), j = 1:(128 - 31))
plot(ij_scan)
compute_dist <- function(i, j){
  sum(abs(newpatch - img[i:(i+10), j:(j+10)]))
}
dmat <- mapply(compute_dist, i = ij_scan$i, j = ij_scan$j) 
imin <- which.min(dmat)

patch_C <- with(ij_scan[imin, ], img[i:(i+10), j:(j+10)])
par(mfrow = c(3, 1))
image(newpatch, col = rev(grey.colors(256)), zlim = c(0, 256))
image(patch_C, col = rev(grey.colors(256)), zlim = c(0, 256))
image(abs(newpatch - patch_C), col = rev(grey.colors(256)), zlim = c(0, 256))

# generate graph data
diff_ABC <- abs(newpatch - patch_C)
graph_df_newpath <- lapply(
  1:length(diff_ABC),
  create_egdes,
  diff_patch = diff_ABC
) %>%
  do.call(rbind, .) %>%
  subset(from < to)

