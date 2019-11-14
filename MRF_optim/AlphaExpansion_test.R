library(ncdf4)
library(fields)
library(Rcpp)

setwd("~/2_Code/Thao/GraphCut/MRF_optim/")
sourceCpp("../lib/gco-v3.0/example_rcpp.cpp", rebuild = TRUE)

colorTable <- designer.colors (32, c( " blue "," grey 90", " red "))
nc <- nc_open("t2_erai_clim_1979_2008.nc")
tas_ref <- ncvar_get(nc, "t2")
tas_ref <- tas_ref[, ncol(tas_ref):1]
tas_ref <- rbind(tas_ref[454:480,167:221],tas_ref[1:81,167:221])
nc_close(nc)
# image.plot(tas_ref, col = colorTable, zlim = c(200, 300))


nc <- nc_open("tas_cnrm_clim_1979_2008.nc")
tas_lab1 <- ncvar_get(nc, "tas")
tas_lab1 <- tas_lab1[, ncol(tas_lab1):1]
tas_lab1 <- rbind(tas_lab1[454:480,167:221],tas_lab1[1:81,167:221])
nc_close(nc)
image.plot(tas_lab1, col = colorTable, zlim = c(260, 300))

nc <- nc_open("tas_ipsl_clim_1979_2008.nc")
tas_lab2 <- ncvar_get(nc, "tas")
tas_lab2 <- tas_lab2[, ncol(tas_lab2):1]
tas_lab2 <- rbind(tas_lab2[454:480,167:221],tas_lab2[1:81,167:221])
nc_close(nc)
image.plot(tas_lab2, col = colorTable, zlim = c(260, 300))

nc <- nc_open("tas_mpi_clim_1979-01-01,2008-12-31.nc")
tas_lab3 <- ncvar_get(nc, "tas")
tas_lab3 <- tas_lab3[, ncol(tas_lab3):1]
tas_lab3 <- rbind(tas_lab3[454:480,167:221],tas_lab3[1:81,167:221])
nc_close(nc)
image.plot(tas_lab3, col = colorTable, zlim = c(260, 300))

tas_labs <- array(c(tas_lab1, tas_lab2, tas_lab3), dim = c(dim(tas_lab1), 3))
bias <- tas_labs - c(tas_ref)

width = ncol(tas_ref)
height = nrow(tas_ref)
nlabs = dim(tas_labs)[3]
# debug(alphaExpansion_wrapper)
ae_cut = alphaExpansion_wrapper(bias, niter = 10)
mat2format = matrix(tas_labs, ncol = nlabs)
mat = matrix(sapply(seq.int(width * height), function(i){mat2format[i, ae_cut$labels[i] + 1]}), ncol = width, nrow = height) 
image.plot(tas_ref, col = colorTable, zlim = c(260, 300), main = "ref")
image.plot(mat, col = colorTable, zlim = c(260, 300), main = "graphcut")
image.plot(mat - tas_ref, col = colorTable, zlim = c(-10, 10), main = "residuals")
mat = matrix(ae_cut$labels + 1, ncol = ncol(tas_lab1))
image.plot(mat, nlevel = nlabs)
