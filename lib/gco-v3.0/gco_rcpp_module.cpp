#include "GCoptimization_rcpp.h"
#include <Rcpp.h>
// RCPP_MODULE(PACE){
//   class_<Base>("Base")
//   .constructor<double>()
//   .field("x", &Base::x) 
//   ;
//   class_<Derived>("Derived")
//     .derives<Base>("Base")
//     .constructor<int>()
//     .field("y", &Derived::y)
//   ;
RCPP_MODULE(gco) {
  using namespace Rcpp;
  // we expose the class std::vector<double> as "vec" on the R side
  class_<GCoptimization>("GCoptimization")
    .method( "setDataCost", &GCoptimization::setDataCostRcpp )
    .method( "setSmoothCost", &GCoptimization::setSmoothCostRcpp )
    .method( "expansion", &GCoptimization::expansion )
    .method( "swap", &GCoptimization::swap )
    .method( "compute_energy", &GCoptimization::compute_energy )
    .method( "giveDataEnergy", &GCoptimization::giveDataEnergy )
    .method( "giveSmoothEnergy", &GCoptimization::giveSmoothEnergy )
    .method( "whatLabel", &GCoptimization::whatLabelRcpp )
    .method( "setLabel", &GCoptimization::setLabel )
    ;
  class_<GCoptimizationGridGraph>("GCoptimizationGridGraph")
    // exposing constructors
    .derives<GCoptimization>("GCoptimization")
    .constructor<GCoptimization::SiteID, GCoptimization::SiteID, GCoptimization::LabelID>()
    ;
}

/*** R
library(RcppXPtrUtils)

nlabs = 3
width = height = 8
ref = lab1 = lab2 = lab3 = matrix(0, nrow = height, ncol = width)
ref[1:3, 1:3] =  1
lab1[1:4, 1:4] =  1
lab2[1:2, 1:2] =  1
lab3[1:3, 1:3] =  1
lab3[8, 8] = 1

labs = array(c(lab1, lab2, lab3), dim = c(height, width, nlabs))
bias = labs - c(ref)
print(bias)
bias = c(aperm(bias, c(2, 1, 3)))

gco <- new(GCoptimizationGridGraph, width, height, 3)
ptrDataCost <- cppXPtr(
  code = 'GCoptimization::EnergyTermType dataFn(int p, int l, Rcpp::List extraData)
{
  int numPix = extraData["numPix"];
  NumericVector data = extraData["data"];
  return( 1.0 * pow(data[p + numPix * l], 2) );
}', 
  includes = c("#include <math.h>", "#include <Rcpp.h>", '#include "GCoptimization_rcpp.h"')
)
ptrSmoothCost <- cppXPtr(
  code = 'GCoptimization::EnergyTermType smoothFn(int p1, int p2, int l1, int l2, Rcpp::List extraData)
{
  int numPix = extraData["numPix"];
  NumericVector data = extraData["data"];
  float cost = abs(data[p1 + numPix * l1]-data[p1 + numPix * l2]) + abs(data[p2 + numPix * l1] - data[p2 + numPix * l2]) ;
  return(0.5 * cost);
}', 
  includes = c("#include <math.h>", "#include <Rcpp.h>", '#include "GCoptimization_rcpp.h"')
)

gco$setDataCost(ptrDataCost, list(numPix = width * height, data = bias))
gco$setSmoothCost(ptrSmoothCost, list(numPix = width * height, data = bias))
cat("energy before optimization: ", gco$compute_energy(), "\n with labels: \n")
for(j in 1:height){
  for(i in 1:width){
    cat(gco$whatLabel((i - 1) + width * (j - 1)), " ")
  }
  cat("\n")
}
# gco$expansion(-1)
gco$swap(-1)
cat("energy after optimization: ", gco$compute_energy(), "\n with labels: \n")
for(j in 1:height){
  for(i in 1:width){
    cat(gco$whatLabel((i - 1) + width * (j - 1)), " ")
  }
  cat("\n")
}
*/
