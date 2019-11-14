library(ncdf4)
library(fields)
LABELS = 3
height = 100
width = 100
dataref <- matrix(1, nrow = height, ncol = width)
dataref[1:height > 25 & 1:height < 50, 1:width > 25 & 1:width < 50] <- 2
dataref[1:height > 50 & 1:height < 75, 1:width > 50 & 1:width < 75] <- 3
datanoise <- dataref
noisypixels <- sample.int(height * width, 2000)
datanoise[noisypixels] <- round(runif(length(noisypixels), min = 0.5, max = 3.5))
image.plot(datanoise)
datanoise <- as.vector(datanoise)

create_smoothcost <- function(datanoise){
  function(labx, laby, scale){
    return(1 * (labx != laby) / scale)
  }
}

create_datacost <- function(datanoise){
  function(x, labx, scale){
    return(1 * (labx != datanoise[x]) / scale)
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



init_MRF <- function(datanoise, LABELS){
  datacost <- create_datacost(datanoise)
  smoothcost <- create_smoothcost(datanoise)
  MRF <- list(
    LABELS = LABELS,
    height = height, width = width,
    theta = datanoise,
    best = datanoise,
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
      probs[lab] <- probs[lab] + smoothcost(lab, MRF$theta[nup], scale)
    if(!is.null(ndown))
      probs[lab] <- probs[lab] + smoothcost(lab, MRF$theta[ndown], scale)
    if(!is.null(nleft))
      probs[lab] <- probs[lab] + smoothcost(lab, MRF$theta[nleft], scale)
    if(!is.null(nright))
      probs[lab] <- probs[lab] + smoothcost(lab, MRF$theta[nright], scale)
  }
  # print(probs)
  probs <- -probs
  probs <- probs - max(probs)
  probs <- exp(probs)
  # print(probs)
  if(all(probs == 0)) probs <- rep(1, MRF$LABELS)
  probs <- probs / sum(probs)
  multinorm_draw <- rmultinom(1, 1, probs)
  # print(multinorm_draw)
  # print("-------------------------")
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
      energy = energy + smoothcost(cur_label, MRF$theta[nup], scale)
    if(!is.null(ndown))
      energy = energy + smoothcost(cur_label, MRF$theta[ndown], scale)
    if(!is.null(nleft))
      energy = energy + smoothcost(cur_label, MRF$theta[nleft], scale)
    if(!is.null(nright))
      energy = energy + smoothcost(cur_label, MRF$theta[nright], scale)
  }
  return(energy)
}

MRF <- init_MRF(datanoise, LABELS = LABELS)
print(compute_energy(MRF, scale = 1))
MRF <- simulated_annealing(MRF, scale = 1 * seq(1, 0.01, -0.01), period = 10)
mat = matrix(MRF$best, ncol = height)
image.plot(mat)
mat = matrix(MRF$theta, ncol = height)
image.plot(mat)
plot(MRF$energy)
