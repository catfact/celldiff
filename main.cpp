/*
 *  main.cpp
 *  cellDiffusion
 *
 *  Created by Ezra Buchla on 10/6/11.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <ncurses.h>
#include "CellModel.hpp"

static u32 n = 64;

static void printFrame(CellModel* model, u64 frame) {
  u32 i = frame;
  u64 idx;
   for(u32 j=0; j<n; j++)
    {
      // print states
      for(u32 k=0; k<n; k++) {
      idx = i*n*n + j*n + k;
        //printf("%d ", model->cells[(i*n*n) + (j*n) + k]->state);
        attron(COLOR_PAIR(model->cells[idx]->state + 1));
        mvprintw(2+j, k * 2, "%d ", model->cells[idx]->state);
        attroff(COLOR_PAIR(model->cells[idx]->state + 1));
      }
    }
    refresh();
}

int main (const int argc, const char** argv) {

  u32 count;
  u32 frameCount, frameNum;
  f32 pd, pe, pp;
  f32 releaseRatio;  
  
  //////// ncurses stuff
  initscr();
  raw();
  noecho();


  
//////// ncurses stuff
  
  initscr();
  
  raw();
  noecho();

mvprintw(4, 0, "celldiff n count pDrug pExcip pPoly frameInterval frameSlice");
  ////// parse arguments
  // max iterations count
  if (argc < 2) {
    count =  100;
  } else {
    count = (u64)atoi(argv[1]);
  }
 
 // cube length
  if (argc < 3) {
   n = 32;
  } else {
    n = (u64)atof(argv[2]);
  }

// drug ratio
  if (argc < 4) {
    pd= 0.25;
  } else {
    pd = (f32)atof(argv[3]);
  }

// excipient ratio
  if (argc < 5) {
    pe = 0.25;
  } else {
    pe = (f32)atof(argv[4]);
  }

// polymer ratio
  if (argc < 6) {
    pp = 0.5f;
  } else {
    pp = (f32)atof(argv[5]);
  }

// animation frame interval
  if (argc < 7) {
    frameCount = 1;
  } else {
    frameCount = (u32)atoi(argv[6]);
  }

// animation frame slice
  if (argc < 8) {
    frameNum = n >> 1;
  } else {
    frameNum = (u32)atoi(argv[7]);
  }


if(has_colors() == FALSE) {
   endwin();
   printf("Your terminal does not support color\n");
   return 1;
  }
  
  start_color();
  
	init_pair(1, COLOR_BLACK, COLOR_WHITE);
	init_pair(2, COLOR_BLACK, COLOR_RED);
	init_pair(3, COLOR_BLACK, COLOR_GREEN);
	init_pair(4, COLOR_BLACK, COLOR_YELLOW);
	init_pair(5, COLOR_BLACK, COLOR_BLUE);
	init_pair(6, COLOR_BLACK, COLOR_MAGENTA);
	init_pair(7, COLOR_BLACK, COLOR_CYAN);
	init_pair(8, COLOR_WHITE, COLOR_BLACK);

  CellModel model(
		  n,      // cube width
		  pd,    // p drug
		  pe,    // p excipient
		  pp,    // p polymer
		  0.001,  // cell size
		  0.001,  // time step
		  7e-6,   // drug diffusion constant
		  7e-6,   // excipient diffusion constant
		  47u     // RNG seed
		  );

//  printf("setting up model... cube width %i, pd: %f, pe: %f, pp: %f\n\n", (int)n, pd, pe, pp);
	  mvprintw(0, 0, "cube width %i, pd: %f, pe: %f, pp: %f\n\n", (int)n, pd, pe, pp);
	  mvprintw(1, 0, "initializing...");
	refresh();	
  model.setup();



    mvprintw(4, 0, "(usage: celldiff count n pDrug pExcip pPoly frameInt frameSlice) \n\n");

  mvprintw(1, 0, "performing %d iterations on %d cells. press any key to continue...\n\n", count, n*n*n);
  refresh();
	getchar();
	mvprintw(1, 0, "                                                                            ");
	refresh();


  int step = 0;
	u32 frameStep = 0;
  while(step < count)
  {
    step++;
    frameStep++;
    if(frameStep == frameCount) {
    printFrame(&model, frameNum);
    frameStep = 0;
    }
    
    
    releaseRatio = model.iterate();
    mvprintw(1, 0, "finished iteration %d of %d, released ratio: %f, intial mass: %f", step, count, releaseRatio, model.drugMassTotal);
    refresh();
  }
  printFrame(&model, frameNum); 
 	
    mvprintw(2, 0, "finished simulation. press any key to quit...                                                ");
    refresh();
 	getchar(); 
  endwin();
  return 1;
}
