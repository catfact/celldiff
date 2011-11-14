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

//======= defines
#define HALT_NO_CHANGE 1
#define HALT_MAX_ITERATIONS 2

//=======  vars
static const f64 noChangeMassThresh = 0.00001f;
static const u32 noChangeCountThresh = 500;
static u32 n = 64;
// current iteration count
static  u32 count;
static  u32 frameCount, frameNum;
// concentrations
static  f64 pd, pe, pp, rt;

//============== function declarations
int main(const int argc, const char** argv);

static void parse_args(const int argc, const char** argv);
static void start_graphics(void);
static void end_graphics(void);
static void printFrame(CellModel* model, u32 frame);


//============== function definitions

int main (const int argc, const char** argv) {
	
  // released drug mass 
  f64 released[2] = {-1000.f, 0.f};
  // change in released mass
  f64 dr;
  u64 noChangeCount = 0;

  parse_args(argc, argv);
  start_graphics();
	
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
	
  mvprintw(0, 0, "cube width %i, pd: %f, pe: %f, pp: %f\n\n", (int)n, pd, pe, pp);
  mvprintw(1, 0, "initializing...");
  refresh();	
  model.setup();
	
  mvprintw(1, 0, "(usage: celldiff count cubeLength pPoly releaseThresh frameInt frameSlice) \n\n");
  mvprintw(2, 0, "cell memory is %d bytes\n\n", model.numCells * sizeof(Cell));
  mvprintw(3, 0, "performing %d iterations on %d cells. press any key to continue...\n\n", count, n*n*n);
	
  refresh();
  getchar();
  mvprintw(1, 0, "                                                                            ");
  mvprintw(2, 0, "                                                                            ");
  mvprintw(3, 0, "                                                                            ");
  mvprintw(4, 0, "                                                                            ");
  refresh();
	
	
  int step = 0;
  u32 frameStep = 0;
  u8 halt = 0;
  //  while((step < count) && (releaseRatio < rt) )
  while(halt == 0)    
{
      step++;
      frameStep++;

      if ( step == count ) {
	halt = HALT_MAX_ITERATIONS;
      }

      if(frameStep == frameCount) {
	printFrame(&model, frameNum);
	frameStep = 0;
      }
			
		
      released[1] = model.iterate();
      dr = released[1] - released[0];
      released[0] = released[1];
  
      if(dr < noChangeMassThresh) {
	if (released[1] > 0.01) {
	  noChangeCount++;
	}
      }
      else {
	noChangeCount = 0;
      }
      if(noChangeCount == noChangeCountThresh) {
	halt = HALT_NO_CHANGE;
      }

      mvprintw(1, 0, "finished iteration %d of %d, released %f of %f, ratio %f", step, count, released[1], model.drugMassTotal, released[1] / model.drugMassTotal);
		
      refresh();
 } // end main loop
	
	
	
  printFrame(&model, frameNum); 
  
  switch(halt) {
  case HALT_MAX_ITERATIONS:
mvprintw(2, 0, "finished maximum iterations; simulation halted. press any key to quit...              ");
    break;
  case HALT_NO_CHANGE:

    mvprintw(2, 0, "released mass appears stable; simulation halted . press any key to quit...           ");
 break;
  }
  
  refresh();
  getchar();

  end_graphics(); 
  return 1;
}

//-------- graphics init
static void start_graphics(void) {
  initscr();
  raw();
  noecho();
	
  if(has_colors() == FALSE) {
    endwin();
    printf("Your terminal does not support color!\n");
  }
      
  start_color();
	
  init_pair(1, COLOR_BLACK, COLOR_WHITE);
  init_pair(2, COLOR_BLACK, COLOR_CYAN);
  init_pair(3, COLOR_BLACK, COLOR_GREEN);
  init_pair(4, COLOR_BLACK, COLOR_YELLOW);
  init_pair(5, COLOR_BLACK, COLOR_BLUE);
  init_pair(6, COLOR_BLACK, COLOR_MAGENTA);
  init_pair(7, COLOR_BLACK, COLOR_RED);
  init_pair(8, COLOR_WHITE, COLOR_BLACK);
}

//------ graphics de-init
static void end_graphics(void) {
  endwin();
}

//------ parse arguments
static void parse_args(const int argc, const char** argv) {
  ////// parse arguments
  // max iterations count
  if (argc < 2) {
    count =  100;
  } else {
    count = (u32)atoi(argv[1]);
  }
	
  // cube length
  if (argc < 3) {
    n = 32;
  } else {
    n = (u32)atof(argv[2]);
  }
 
  // polymer ratio
  if (argc < 4) {
    pp = 0.5f;
  } else {
    pp = (f64)atof(argv[3]);
  }
  
  // drug concentration is fixed
  pd  = 0.1;
  pe = 1.f - pd - pp;

  // release threshold
  if (argc < 5) {
    rt = 0.5f;
  } else {
    rt = (f64)atof(argv[4]);
  }
	
  // animation frame interval
  if (argc < 6) {
    frameCount = 1;
  } else {
    frameCount = (u32)atoi(argv[5]);
  }
	
  // animation frame slice
  if (argc < 7) {
    frameNum = n >> 1;
  } else {
    frameNum = (u32)atoi(argv[6]);
  }
}

static void printFrame(CellModel* model, u32 frame) {
  u32 i = frame;
  u32 idx;

  for(u32 j=0; j<n; j++)
    {
      // print states
      for(u32 k=0; k<n; k++) {
	idx = i*n*n + j*n + k;
	attron(COLOR_PAIR(model->cells[idx]->state + 1));
	if (model->cells[idx]->state == eStateWet) {
	  mvprintw(6+j, k * 2, "%d.1 ", (int)(model->cells[idx]->concentration[eStateDrug] * 9.5));
	} else {
	  mvprintw(6+j, k * 2, "%d ", model->cells[idx]->state);
	}
	attroff(COLOR_PAIR(model->cells[idx]->state + 1));
      }
    }
  refresh();
}
