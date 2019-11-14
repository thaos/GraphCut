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


struct ForDataFn{
	int numPix;
	float *data;
};

struct ForSmoothFn{
	int numPix;
	float *data;
};

float smoothFn(int p1, int p2, int l1, int l2, void *data)
{
	ForDataFn *myData = (ForDataFn *) data;
	int numPix = myData->numPix;
  float cost = pow(myData->data[p1 + numPix * l1]-myData->data[p1 + numPix * l2], 2) + pow(myData->data[p2 + numPix * l1] - myData->data[p2 + numPix * l2], 2) ;
  return(0.5 * cost);
}

float dataFn(int p, int l, void *data)
{
	ForDataFn *myData = (ForDataFn *) data;
	int numPix = myData->numPix;
	return( 1.0 * pow(myData->data[p + numPix * l], 2) );
}



////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood structure is assumed
//
void GridGraph_DfnSfn(int width,int height,int num_pixels,int num_labels)
{

	int *result = new int[num_pixels];   // stores result of optimization
  float *ref = new float[num_pixels];
  float *labs = new float[num_pixels * num_labels];
  float *bias = new float[num_pixels * num_labels];
  for(int i = 0; i < width; i++){  
    for(int j = 0; j < height; j++){
      if(i < 3 && j < 3){
        ref[i + width * j] = 1;
      } 
      else {
        ref[i + width * j] = 0;
      }
      for(int k = 0; k < num_labels; k++){
        if(k == 0){
          if(i < 2 && j < 2){
            labs[i + width * j + num_pixels * k] = 1;
          } 
          else {
            labs[i + width * j + num_pixels * k] = 0;
          }
        }
        if (k == 1) {
          if(i < 4 && j < 4){
            labs[i + width * j + num_pixels * k ] = 1;
          } 
          else {
            labs[i + width * j + num_pixels * k] = 0;
          }
        }
        if (k == 2){
          if(i < 3 && j < 3){
            labs[i + width * j + num_pixels * k] = 1;
          } 
          else {
            labs[i + width * j + num_pixels * k] = 0;
          }
          if( i == 7 && j == 7){ 
            labs[i + width * j + num_pixels * k] = 1;
          }
        }
        bias[i + width * j + num_pixels * k] = labs[i + width * j + num_pixels * k] - ref[i + width * j];
      }
    }
  }
  printf("\n Ref: \n");
  for(int i = 0; i < width; i++){  
    for(int j = 0; j < height; j++){
      printf("%.0f", ref[i + width * j]);
    }
    printf("\n");
  }
  printf("\n lab1: \n");
  for(int i = 0; i < width; i++){  
    for(int j = 0; j < height; j++){
      printf("%.0f", labs[i + width * j]);
    }
    printf("\n");
  }
  printf("\n lab2: \n");
  for(int i = 0; i < width; i++){  
    for(int j = 0; j < height; j++){
      printf("%.0f", labs[i + width * j + num_pixels]);
    }
    printf("\n");
  }
  printf("\n lab3: \n");
  for(int i = 0; i < width; i++){  
    for(int j = 0; j < height; j++){
      printf("%.0f", labs[i + width * j + 2*num_pixels]);
    }
    printf("\n");
  }

	try{
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width,height,num_labels);

		// set up the needed data to pass to function for the data costs
		ForDataFn toDataFn;
		toDataFn.data = bias;
		toDataFn.numPix = num_pixels;

		ForSmoothFn toSmoothFn;
		toSmoothFn.data = bias;
		toSmoothFn.numPix = num_pixels;

    printf("\n bias1: \n");
    for(int i = 0; i < width; i++){  
      for(int j = 0; j < height; j++){
        printf("%.0f", dataFn(i + width * j, 0, &toDataFn));
        // printf("%f \t", bias[i + width * j]);
      }
      printf("\n");
    }
    printf("\n bias2: \n");
    for(int i = 0; i < width; i++){  
      for(int j = 0; j < height; j++){
        printf("%.0f", dataFn(i + width * j, 1, &toDataFn));
      }
      printf("\n");
    }
    printf("\n bias3: \n");
    for(int i = 0; i < width; i++){  
      for(int j = 0; j < height; j++){
        printf("%.0f", dataFn(i + width * j, 2, &toDataFn));
      }
      printf("\n");
    }
		gc->setDataCost(&dataFn,&toDataFn);

		// smoothness comes from function pointer
		gc->setSmoothCost(&smoothFn, &toSmoothFn);

		printf("\nBefore optimization energy is %Lf",gc->compute_energy());
		printf("\n Pixel Labels: \n");
    for(int i = 0; i < width; i++){  
      for(int j = 0; j < height; j++){
        printf("%d", gc->whatLabel(i + width * j));
      }
      printf("\n");
    }
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		printf("\nAfter optimization energy is %Lf",gc->compute_energy());
		printf("\n Pixel Labels: \n");

    for(int i = 0; i < width; i++){  
      for(int j = 0; j < height; j++){
        result[i + width * j] = gc->whatLabel(i + width * j);
        printf("%d", result[i + width * j]);
      }
      printf("\n");
    }
		printf("\n Output: \n");

    for(int i = 0; i < width; i++){  
      for(int j = 0; j < height; j++){
        printf("%.0f", labs[i + width * j + num_pixels * result[i + width * j]]);
      }
      printf("\n");
    }
    delete gc;
  }
	catch (GCException e){
		e.Report();
	}

	delete [] result;
	delete [] bias;
	delete [] labs;
	delete [] ref;

}
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	int width = 80;
	int height = 80;
	int num_pixels = width*height;
	int num_labels = 3;


	// smoothness and data costs are set up using functions
	GridGraph_DfnSfn(width,height,num_pixels,num_labels);
	
	printf("\n  Finished %d (%d) clock per sec %d",clock()/CLOCKS_PER_SEC,clock(),CLOCKS_PER_SEC);

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////

