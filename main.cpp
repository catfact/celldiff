/*
 *  main.cpp
 *  cellDiffusion
 *
 *  Created by Ezra Buchla on 10/6/11.
 *
 */

//#include <iostream>
#include <cstdio>
#include "CellModel.hpp"

// cube width
static const u64 n=8;

static void printState(CellModel* model)
{
  for(u32 i=0; i<n; i++)
  {
    for(u32 j=0; j<n; j++)
    {
      for(u32 k=0; k<n; k++) {
        printf("%d ", model->cells[(i*n*n) + (j*n) + k]->state);
      }
      printf("\n");
    }
    printf("\n");
  }  
}  

int main (const int argc, const char** argv)
{
  CellModel model(
		  n,      // cube width
		  0.2,    // p drug
		  0.4,    // p excipient
		  0.4,    // p polymer
		  0.001,  // cell size
		  0.001,  // time step
		  7e-6,   // drug diffusion constant
		  7e-6,   // excipient diffusion constant
		  47u     // RNG seed
		  );

  model.setup();
  printState(&model);
  
  int step = 0;
  while(1)
  {
    getchar();
    step++;
    model.iterate();
    printState(&model);
    printf("finished iteration %d\n", step);
  }
}
