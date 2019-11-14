library(ncdf4)
library(fields)
LABELS = 3
height = 100
width = 100
dataref <- matrix(1, nrow = height, ncol = width)
dataref[1:height > 25 & 1:height < 50, 1:width > 25 & 1:width < 50] <- 2
dataref[1:height > 50 & 1:height < 75, 1:width > 50 & 1:width < 75] <- 3
datanoise <- dataref
noisypixels <- sample.int(height * width, 1000)
datanoise[noisypixels] <- round(runif(length(noisypixels), min = 0.5, max = 3.5))
image.plot(datanoise)
datanoise <- as.vector(datanoise)

create_smoothcost <- function(datanoise){
  function(labx, laby){
    return(1 * (labx != laby))
  }
}

create_datacost <- function(datanoise){
  function(x, labx){
    return(1 * (labx != datanoise[x]))
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
    "grid" = array(0, dim = c(length(datanoise), 5, LABELS)),
    height = height, width = width,
    best = datanoise,
    datacost = datacost,
    smoothcost = smoothcost
  )
  for(node_p in 1:(MRF$height *  MRF$width)){
    for(xp in 1:LABELS){
      MRF$grid[node_p, 5, xp] <- datacost(node_p, xp)
    }
  }
  return(MRF)
}

send_msg <- function(MRF, node_q, direction){
  LABELS <- dim(MRF$grid)[3]
  new_msg <- numeric(LABELS)
  width <- MRF$width
  height = MRF$height
  smoothcost <- MRF$smoothcost
  node_p <- switch(direction,
         get_neighbour_up(node_q, height, width),
         get_neighbour_down(node_q, height, width),
         get_neighbour_left(node_q, height, width),
         get_neighbour_right(node_q, height, width)
  )
  if(is.null(node_p)) return(MRF)
  # if(node_q == 5) browser()
  for(xq in 1:LABELS) {
    min_val <- Inf
    for(xp in 1:LABELS) {
      mpq <-  MRF$grid[node_p, 5, xp]
      mpq <- mpq + smoothcost(xp, xq)
      #Exclude the incoming message direction that we are sending to
      if(direction != 1) mpq <- mpq + MRF$grid[node_p, 1, xp]
      if(direction != 2) mpq <- mpq + MRF$grid[node_p, 2, xp]
      if(direction != 3) mpq <- mpq + MRF$grid[node_p, 3, xp]
      if(direction != 4) mpq <- mpq + MRF$grid[node_p, 4, xp]

      min_val = min(min_val, mpq)
    }
    new_msg[xq] = min_val
  }
  MRF$grid[node_q, direction, ] = new_msg
  return(MRF)
}

BP <- function(MRF, direction){
  stopifnot(direction %in% 1:4)
  for(node_p in seq_along(MRF$best)){
    MRF <- send_msg(MRF, node_p, direction);
  }
  # print(MRF$grid[20,,])
  return(MRF)
}

compute_energy <- function(MRF){
  width = MRF$width;
  height = MRF$height;
  smoothcost <- MRF$smoothcost
  # Energy
  energy = 0;
  for(node_p in seq_along(MRF$best)){
    cur_label = MRF$best[node_p]
    # Data cost
    energy = energy + MRF$grid[node_p, 5, cur_label];
    nup <- get_neighbour_up(node_p, height, width)
    ndown <- get_neighbour_down(node_p, height, width)
    nleft <- get_neighbour_left(node_p, height, width)
    nright <- get_neighbour_right(node_p, height, width)
    if(!is.null(nup))
      energy = energy + smoothcost(cur_label, MRF$best[nup])
    if(!is.null(ndown))
      energy = energy + smoothcost(cur_label, MRF$best[ndown])
    if(!is.null(nleft))
      energy = energy + smoothcost(cur_label, MRF$best[nleft])
    if(!is.null(nright))
      energy = energy + smoothcost(cur_label, MRF$best[nright])
  }
  return(energy)
}

MAP <- function(MRF){
  #Finds the MAP assignment as well as calculating the energy
  #MAP assignment
  LABELS <- dim(MRF$grid)[3]
  for(node_p in seq_along(MRF$best)) {
    best = -Inf
    for(xp in 1:LABELS) {
      bp <-  -sum(MRF$grid[node_p, , xp])
       if(bp > best) {
        best = bp;
        MRF$best[node_p] = xp;
      }
    }
  }
  return(MRF);
}

MRF <- init_MRF(datanoise, LABELS = LABELS)
print(compute_energy(MRF))
BP_ITERATIONS <- 10
for(i in 1:BP_ITERATIONS) {
  # if (i == 2) debug(send_msg)
  for(direction in 1:4){
    MRF <- BP(MRF, direction)
    # browser()
  }
  MRF <- MAP(MRF)
  cat("iteration ", (i), "/", BP_ITERATIONS, ", energy = ", compute_energy(MRF), "\n")
}
mat = matrix(MRF$best, ncol = height)
image.plot(mat)

