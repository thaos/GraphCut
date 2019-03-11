library(magick)
library(magrittr)
library(igraph)
library(optrees)

source("graph_cut_algo.R")

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

canvas <- canvas_origin <- matrix(NA, ncol = 32, nrow = 32*2 - 8)
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
  training_img = img, patch_list = NULL, cutset_global = NULL
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas$cutset_global, canvas = new_canvas$canvas)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas$patch_list)))
lines_seams(cutset_global = new_canvas$cutset_global, canvas = new_canvas$canvas)

new_canvas2 <- add_newpatch(
  xstart = 20, ystart = 23,
  xlength = 20, ylength = 10,
  canvas = new_canvas$canvas, canvas_origin =  new_canvas$canvas_origin, canvas_id = new_canvas$canvas_id,
  training_img = img, patch_list = new_canvas$patch_list, cutset_global = new_canvas$cutset_global
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas2$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas2$cutset_global, canvas = new_canvas2$canvas)
abline(v = 19.5)
abline(h = 22.5)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas2$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas2$patch_list)))
lines_seams(cutset_global = new_canvas2$cutset_global, canvas = new_canvas2$canvas)

new_canvas3 <- add_newpatch(
  xstart = 20, ystart = 14,
  xlength = 20, ylength = 10,
  canvas = new_canvas2$canvas, canvas_origin =  new_canvas2$canvas_origin, canvas_id = new_canvas2$canvas_id,
  training_img = img, patch_list = new_canvas2$patch_list, cutset_global = new_canvas2$cutset_global
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas3$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas3$cutset_global, canvas = new_canvas3$canvas)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas3$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas3$patch_list)))
lines_seams(cutset_global = new_canvas3$cutset_global, canvas = new_canvas3$canvas)
      
new_canvas4 <- add_newpatch(
  xstart = 20, ystart = 4,
  xlength = 20, ylength = 10,
  canvas = new_canvas3$canvas, canvas_origin =  new_canvas3$canvas_origin, canvas_id = new_canvas3$canvas_id,
  training_img = img, patch_list = new_canvas3$patch_list, cutset_global = new_canvas3$cutset_global
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas4$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas4$cutset_global, canvas = new_canvas4$canvas)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas4$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas4$patch_list)))
lines_seams(cutset_global = new_canvas4$cutset_global, canvas = new_canvas4$canvas)


new_canvas5 <- add_newpatch(
  xstart = 20, ystart = 1,
  xlength = 20, ylength = 10,
  canvas = new_canvas4$canvas, canvas_origin =  new_canvas4$canvas_origin, canvas_id = new_canvas4$canvas_id,
  training_img = img, patch_list = new_canvas4$patch_list, cutset_global = new_canvas4$cutset_global
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas5$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas5$cutset_global, canvas = new_canvas5$canvas)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas5$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas5$patch_list)))
lines_seams(cutset_global = new_canvas5$cutset_global, canvas = new_canvas5$canvas)

new_canvas6 <- add_newpatch(
  xstart = 30, ystart = 1,
  xlength = 27, ylength = 32,
  canvas = new_canvas5$canvas, canvas_origin =  new_canvas5$canvas_origin, canvas_id = new_canvas5$canvas_id,
  training_img = img, patch_list = new_canvas5$patch_list, cutset_global = new_canvas5$cutset_global
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas6$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas6$cutset_global, canvas = new_canvas6$canvas)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas6$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas6$patch_list)))
lines_seams(cutset_global = new_canvas6$cutset_global, canvas = new_canvas6$canvas)

new_canvas7 <- add_newpatch(
  xstart = 10, ystart = 5,
  xlength = 16, ylength = 24,
  canvas = new_canvas6$canvas, canvas_origin =  new_canvas6$canvas_origin, canvas_id = new_canvas6$canvas_id,
  training_img = img, patch_list = new_canvas6$patch_list, cutset_global = new_canvas6$cutset_global
)
par(mfrow = c(2, 1))
image(1:nrow(canvas), 1:ncol(canvas), new_canvas7$canvas, zlim = c(0, 256), col = grey.colors(256))
lines_seams(cutset_global = new_canvas7$cutset_global, canvas = new_canvas6$canvas)
image(1:nrow(canvas), 1:ncol(canvas), matrix(as.numeric(new_canvas7$canvas_origin), ncol = ncol(canvas)), col = rainbow(nrow(new_canvas7$patch_list)))
lines_seams(cutset_global = new_canvas7$cutset_global, canvas = new_canvas7$canvas)
