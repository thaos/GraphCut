library(ncdf4)
library(fields)
colorTable <- designer.colors (32, c( " blue "," grey 90", " red "))
nc <- nc_open("t2_erai_clim_1979_2008.nc")
tas_ref <- ncvar_get(nc, "t2")
tas_ref <- tas_ref[, ncol(tas_ref):1]
tas_ref <- rbind(tas_ref[454:480,167:221],tas_ref[1:81,167:221])
nc_close(nc)
image.plot(tas_ref, col = colorTable, zlim = c(200, 300))
tas_ref <- as.vector(tas_ref)

nc <- nc_open("tas_cnrm_clim_1979_2008.nc")
tas_lab1 <- ncvar_get(nc, "tas")
tas_lab1 <- tas_lab1[, ncol(tas_lab1):1]
tas_lab1 <- rbind(tas_lab1[454:480,167:221],tas_lab1[1:81,167:221])
nc_close(nc)

nc <- nc_open("tas_ipsl_clim_1979_2008.nc")
tas_lab2 <- ncvar_get(nc, "tas")
tas_lab2 <- tas_lab2[, ncol(tas_lab2):1]
tas_lab2 <- rbind(tas_lab2[454:480,167:221],tas_lab2[1:81,167:221])
nc_close(nc)

nc <- nc_open("tas_mpi_clim_1979-01-01,2008-12-31.nc")
tas_lab3 <- ncvar_get(nc, "tas")
tas_lab3 <- tas_lab3[, ncol(tas_lab3):1]
tas_lab3 <- rbind(tas_lab3[454:480,167:221],tas_lab3[1:81,167:221])
nc_close(nc)

# tas_lab <- cbind(c(tas_lab1), c(tas_lab2), c(tas_lab3), c(tas_ref))
tas_lab <- cbind(c(tas_lab1), c(tas_lab2), c(tas_lab3))
image.plot(tas_lab1, col = colorTable, zlim = c(200, 300))
image.plot(tas_lab2, col = colorTable, zlim = c(200, 300))
image.plot(tas_lab3, col = colorTable, zlim = c(200, 300))
LABELS <- ncol(tas_lab)

create_smoothcost <- function(labtovalue){
  function(x, y, labx, laby, scale){
    # return( 1/scale * (labx != laby)  * (labtovalue[x, labx] - labtovalue[y, laby])^2)
    return(0.5/scale * (abs(labtovalue[x, labx] - labtovalue[x, laby]) + abs(labtovalue[y, labx] - labtovalue[y, laby])))
  }
}

create_datacost <- function(labtovalue, refvalue){
  function(x, labx, scale){
    return( 1/scale * (labtovalue[x, labx] - refvalue[x])^2)
  }
}

get_neighbour_up <- function(x, height, width){
  ans <- x - 1
  if((x %% height) == 1){
    ans <- NULL
  }
  return(ans)
}

get_neighbour_down <- function(x, height, width){
  ans <- x + 1
  if((x %% height) == 0){
    ans <- NULL
  }
  return(ans)
}

get_neighbour_left <- function(x, height, width){
  ans <- x - width
  if(((x - 1) %/% height) == 0){
    ans <- NULL
  }
  return(ans)
}

get_neighbour_right <- function(x, height, width){
  ans <- x + width
  if(((x - 1) %/% height) == (width - 1)){
    ans <- NULL
  }
  return(ans)
}
get_neighbour_left(4, 3, 3)

# UP DOWN LEFT RIGHT DATA



init_MRF <- function(labtovalue, refvalue, height, width){
  LABELS <- ncol(labtovalue)
  datacost <- create_datacost(labtovalue, refvalue)
  smoothcost <- create_smoothcost(labtovalue)
  # init <- round(runif(nrow(labtovalue), min = 0.5, max = LABELS + 0.5))
  # init <- rep(LABELS, nrow(labtovalue))
  init <- rep(1, nrow(labtovalue))
  MRF <- list(
    LABELS = LABELS,
    height = height, width = width,
    theta = init,
    best = init,
    energy_best = Inf,
    datacost = datacost,
    smoothcost = smoothcost
  )
  return(MRF)
}

iterate_MRF_Gibbs <- function(MRF, node_p, scale){
  width <- MRF$width
  height = MRF$height
  smoothcost <- MRF$smoothcost
  datacost <- MRF$datacost
  probs <- numeric(MRF$LABELS)
  nup <- get_neighbour_up(node_p, height, width)
  ndown <- get_neighbour_down(node_p, height, width)
  nleft <- get_neighbour_left(node_p, height, width)
  nright <- get_neighbour_right(node_p, height, width)
  for(lab in 1:LABELS){
    probs[lab] <- probs[lab] + datacost(node_p, lab, scale)
    if(!is.null(nup))
      probs[lab] <- probs[lab] + smoothcost(node_p, nup, lab, MRF$theta[nup], scale)
    if(!is.null(ndown))
      probs[lab] <- probs[lab] + smoothcost(node_p, ndown, lab, MRF$theta[ndown], scale)
    if(!is.null(nleft))
      probs[lab] <- probs[lab] + smoothcost(node_p, nleft, lab, MRF$theta[nleft], scale)
    if(!is.null(nright))
      probs[lab] <- probs[lab] + smoothcost(node_p, nright, lab, MRF$theta[nright], scale)
  }
  probs_backup <- probs
  probs <- -probs
  # print(min(probs))
  # probs <- probs - max(probs)
  probs <- exp(probs)
  if(all(probs == 0)) probs <- rep(1, MRF$LABELS)
  probs <- probs / sum(probs)
  # print(probs)
  multinorm_draw <- rmultinom(1, 1, probs)
  # if(MRF$theta[node_p] != 4 & which(c(multinorm_draw) == 1) != 4 & scale < 0.80){
  #   print(node_p)
  #   print(probs_backup)
  #   print(probs)
  #   print(multinorm_draw)
  #   print("-------------------------")
    # browser()
  # }
  MRF$theta[node_p] <- which(c(multinorm_draw) == 1)
  # MRF$theta[node_p] <- which.max(probs)
  return(MRF)
}

simulated_annealing <- function(MRF, scales, period){
  energy <- numeric(length(scales) * period + 1)
  i <- 1
  energy[i] <- compute_energy(MRF, scale = 1)
  width <- MRF$width
  height = MRF$height
  for (scale in scales){
    cat("\n", i, "/", length(scales) * period + 1, "\n")
    print(energy[i])
    for (time in 1:period){
      i <- i + 1
      cat(".")
      for(node_p in  1:(width * height)){
      # for(node_p in  sample.int(width * height)){
        MRF <- iterate_MRF_Gibbs(MRF, node_p, scale)
        # gc(verbose = FALSE)
      }
      energy[i] <- compute_energy(MRF, scale = 1)
      if(energy[i] < MRF$energy_best){
        MRF$energy_best <- energy[i]
        MRF$best <- MRF$theta
      }
    }
    mat = matrix(MRF$theta, ncol = height)
    image.plot(mat)
  }
  MRF$energy <- energy
  return(MRF)
}


compute_energy <- function(MRF, scale){
  width = MRF$width;
  height = MRF$height;
  smoothcost <- MRF$smoothcost
  datacost <- MRF$datacost
  theta <- MRF$theta
  # Energy
  energy = 0;
  for(node_p in seq_along(MRF$best)){
    cur_label = MRF$theta[node_p]
    # Data cost
    energy = energy + datacost(node_p, cur_label, scale)
    nup <- get_neighbour_up(node_p, height, width)
    ndown <- get_neighbour_down(node_p, height, width)
    nleft <- get_neighbour_left(node_p, height, width)
    nright <- get_neighbour_right(node_p, height, width)
    if(!is.null(nup))
      energy = energy + smoothcost(node_p, nup, cur_label, MRF$theta[nup], scale)
    if(!is.null(ndown))
      energy = energy + smoothcost(node_p, ndown, cur_label, MRF$theta[ndown], scale)
    if(!is.null(nleft))
      energy = energy + smoothcost(node_p, nleft, cur_label, MRF$theta[nleft], scale)
    if(!is.null(nright))
      energy = energy + smoothcost(node_p, nright, cur_label, MRF$theta[nright], scale)
  }
  return(energy)
}
height <- nrow(tas_lab1)
width <-  ncol(tas_lab1)
MRF <- init_MRF(tas_lab, tas_ref, height = height, width = width)
print(compute_energy(MRF, scale = 1))
MRF <- simulated_annealing(MRF, scales = seq(1, 1E-2, length.out = 100),  period = 500)

plot(MRF$energy)
mat = matrix(MRF$theta, ncol = height)
image.plot(mat)
mat = matrix(MRF$best, ncol = height)
image.plot(mat)
mat = matrix(sapply(seq_along(MRF$best), function(i) tas_lab[i, MRF$best[i]]), ncol = width)
image.plot(mat, col = colorTable, zlim = c(200, 300))
image.plot(mat - tas_ref, col = colorTable, zlim = c(-7, 7))
