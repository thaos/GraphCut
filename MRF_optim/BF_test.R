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

tas_lab <- cbind(c(tas_lab1), c(tas_lab2), c(tas_lab3))#, c(tas_ref) + rnorm(length(tas_ref), sd = 0.1))
image.plot(tas_lab1, col = colorTable, zlim = c(200, 300))
image.plot(tas_lab2, col = colorTable, zlim = c(200, 300))
image.plot(tas_lab3, col = colorTable, zlim = c(200, 300))
LABELS <- ncol(tas_lab)

create_smoothcost <- function(labtovalue){
  function(x, y, labx, laby){
    # return( 1/1 * (labx != laby)  * (labtovalue[x, labx] - labtovalue[y, laby])^2)
    return( 1/1 * (labtovalue[x, labx] - labtovalue[x, laby])^2 + (labtovalue[y, labx] - labtovalue[y, laby])^2)
  }
}

create_datacost <- function(labtovalue, refvalue){
  function(x, labx){
    return((labtovalue[x, labx] - refvalue[x])^2)
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
get_neighbour_left(4, 2, 2)

# UP DOWN LEFT RIGHT DATA



init_MRF <- function(labtovalue, refvalue){
  LABELS <- ncol(labtovalue)
  datacost <- create_datacost(labtovalue, refvalue)
  smoothcost <- create_smoothcost(labtovalue)
  MRF <- list(
    "grid" = array(0, dim = c(length(refvalue), 5, LABELS)),
    height = nrow(tas_lab1), width = ncol(tas_lab1),
    best = rep(2, nrow(tas_lab1) * ncol(tas_lab1)),
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
  height = MRF$height;
  smoothcost <- MRF$smoothcost
  node_p <- switch(direction,
         get_neighbour_up(node_q, height, width),
         get_neighbour_down(node_q, height, width),
         get_neighbour_left(node_q, height, width),
         get_neighbour_right(node_q, height, width)
  )
  if(is.null(node_p)) return(MRF)
  for(xq in 1:LABELS) {
    min_val <- Inf
    for(xp in 1:LABELS) {
      mpq <-  MRF$grid[node_p, 5, xp]
      mpq <- mpq + smoothcost(node_p, node_q, xp, xq)
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
    # if (node_p == 100) print(MRF$grid[node_p,,])
  }
  return(MRF);
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
      energy = energy + smoothcost(node_p, nup, cur_label, MRF$best[nup])
    if(!is.null(ndown))
      energy = energy + smoothcost(node_p, ndown, cur_label, MRF$best[ndown])
    if(!is.null(nleft))
      energy = energy + smoothcost(node_p, nleft, cur_label, MRF$best[nleft])
    if(!is.null(nright))
      energy = energy + smoothcost(node_p, nright, cur_label, MRF$best[nright])
  }
  return(energy)
}

MRF <- init_MRF(tas_lab, tas_ref)
print(compute_energy(MRF))
BP_ITERATIONS <- 100
energy_best <- Inf
for(i in 1:BP_ITERATIONS) {
  # if (i == 2) debug(send_msg)
  for(direction in 4:1){
    MRF <- BP(MRF, direction)
  }
  MRF <- MAP(MRF)
  energy <- compute_energy(MRF)
  cat("iteration ", (i), "/", BP_ITERATIONS, ", energy = ", energy, "\n")
  if(energy < energy_best) {
    energy_best <- energy
    MRF_best <- MRF
  }
}
mat = matrix(sapply(seq_along(MRF_best$best), function(i) tas_lab[i, MRF_best$best[i]]), ncol = ncol(tas_lab1))
image.plot(mat, col = colorTable, zlim = c(200, 300))
image.plot(mat - tas_ref, col = colorTable, zlim = c(-10, 10))
mat = matrix(MRF_best$best, ncol = ncol(tas_lab1))
image.plot(mat, nlevel = LABELS)
