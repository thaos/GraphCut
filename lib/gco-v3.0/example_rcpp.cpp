//////////////////////////////////////////////////////////////////////////////
// Example illustrating the use of GCoptimization.cpp
//
/////////////////////////////////////////////////////////////////////////////
//
//  Optimization problem:
//  is a set of sites (pixels) of width 10 and hight 5. Thus number of pixels is 50
//  grid neighborhood: each pixel has its left, right, up, and bottom pixels as neighbors
//  7 labels
//  Data costs: D(pixel,label) = 0 if pixel < 25 and label = 0
//            : D(pixel,label) = 10 if pixel < 25 and label is not  0
//            : D(pixel,label) = 0 if pixel >= 25 and label = 5
//            : D(pixel,label) = 10 if pixel >= 25 and label is not  5
// Smoothness costs: V(p1,p2,l1,l2) = min( (l1-l2)*(l1-l2) , 4 )
// Below in the main program, we illustrate different ways of setting data and smoothness costs
// that our interface allow and solve this optimizaiton problem

// For most of the examples, we use no spatially varying pixel dependent terms.
// For some examples, to demonstrate spatially varying terms we use
// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with
// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "GCoptimization.h"
#include <Rcpp.h>
using namespace Rcpp;

struct ForEnergyFn{
	int numPix;
	NumericVector data;
};


float smoothFn(int p1, int p2, int l1, int l2, void *data)
{
	ForEnergyFn *myData = (ForEnergyFn *) data;
	int numPix = myData->numPix;
  // float cost = pow(myData->data[p1 + numPix * l1]-myData->data[p1 + numPix * l2], 2) + pow(myData->data[p2 + numPix * l1] - myData->data[p2 + numPix * l2], 2) ;
  float cost = abs(myData->data[p1 + numPix * l1]-myData->data[p1 + numPix * l2]) + abs(myData->data[p2 + numPix * l1] - myData->data[p2 + numPix * l2]) ;
  return(0.5 * cost);
}

float dataFn(int p, int l, void *data)
{
	ForEnergyFn *myData = (ForEnergyFn *) data;
	int numPix = myData->numPix;
	return( 1.0 * pow(myData->data[p + numPix * l], 2) );
}



////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood structure is assumed
//
// [[Rcpp::export]]
List alphaExpansion(
    NumericVector bias,
    int width, int height, int num_labels,
    int niter
){
  Rprintf("\n --- In AlphaExpansion ---\n");
  int num_pixels = width * height;
  float energy;
	IntegerVector result = IntegerVector(num_pixels);   // stores result of optimization
	try{
    Rprintf("\n --- Initializing Graph Starting ---\n");
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width,height,num_labels);
		// set up the needed data to pass to function for the data costs
		ForEnergyFn toEnergyFn;
		toEnergyFn.data = bias;
		toEnergyFn.numPix = num_pixels;

		// data term comes from function pointer
		gc->setDataCost(&dataFn,&toEnergyFn);

		// smoothness comes from function pointer
		gc->setSmoothCost(&smoothFn, &toEnergyFn);
    Rprintf("\n --- Initializing Graph Done ---\n");
    Rprintf("\n ------------------------------------ \n");
		Rprintf("Before optimization energy is %f \n", (float) gc->compute_energy());
		Rprintf("with pixel labels: \n\n");
    for(int i = 0; i < width; i++){
      for(int j = 0; j < height; j++){
        Rprintf("%i \t", gc->whatLabel(i + width * j));
      }
      Rprintf("\n");
    }
    Rprintf("\n ------------------------------------ \n");
		gc->expansion(niter);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		// gc->swap(niter);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		energy = (float) gc->compute_energy();
    Rprintf("\n ------------------------------------ \n");
    Rprintf("After optimization energy is %f \n", energy);
		Rprintf("with pixel labels: \n\n");
    for(int i = 0; i < width; i++){
      for(int j = 0; j < height; j++){
        result[i + width * j] = gc->whatLabel(i + width * j);
        Rprintf("%i \t", result[i + width * j]);
      }
      Rprintf("\n");
    }
    Rprintf("\n ------------------------------------ \n");
  }
	catch (GCException e){
		e.Report();
	}
  return List::create(
      Rcpp::Named("labels", result),
      Rcpp::Named("energy", energy)
  );

}
////////////////////////////////////////////////////////////////////////////////

/*** R
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


alphaExpansion_wrapper = function(bias_array, niter = 2){
  dimension = dim(bias_array)
  width = dimension[2]
  height = dimension[1]
  num_labels = dimension[3]
  bias = c(aperm(bias_array, c(2, 1, 3)))
  ae_res = alphaExpansion(
    bias = bias,
    width = width,
    height = height,
    num_labels = num_labels,
    niter =  niter
  )
  ae_res$labels = matrix(ae_res$labels, ncol = width, nrow = height, byrow = TRUE)
  return(ae_res)
}
ae_cut = alphaExpansion_wrapper(bias)
print(ref)
mat2format = matrix(labs, ncol = nlabs)
print(matrix(sapply(seq.int(width * height), function(i){mat2format[i, ae_cut$labels[i] + 1]}), ncol = width, nrow = height))

*/
