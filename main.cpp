/*
 *  main.cpp
 *  cellDiffusion
 *
 *  Created by Ezra Buchla on 10/6/11.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <ncurses.h>
#include <getopt.h>

#include "CellModel.hpp"

using namespace std;

//======= defines
#define HALT_NO_CHANGE 1
#define HALT_MAX_ITERATIONS 2

//=======  variables
static f64 noChangeMassThresh = 0.00001;
static u32 noChangeCountThresh = 100;
// length of cube
static u32 n = 32;
// current iteration count
static u32 iterationCount = 0;
// current animation frame step
static u32 frameStep = 0;
// animation frame period
static u32 framePeriod = 1;
// current index of animation slice
static u32 frameNum = 0;
// concentrations
static f64 pd=0.1, pe=0.5, pp=0.4;
// cylinder height
static f64 h=0.23;
// RNG seed
static u32 seed = 47;
// release curve output file path
static string releasedPath;
// state output file path
static string statePath;
// state output period (0 == no output)
static u32 statePeriod = 0;
// state output step counter 
static u32 stateStep = 0;
// ascii output toggle
static u32 asciiout = 1;

//============== function declarations
int main(const int argc, char* const* argv);

static int parse_args(const int argc, char* const* argv);
static void start_graphics(void);
static void end_graphics(void);
static void printFrame(CellModel* model, u32 frame);

//============== function definitions

//------ main
int main (const int argc, char* const* argv) {
	
  // released drug mass 
  f64 released[2] = {-1000.f, 0.f};
  // change in released mass
  f64 dr;
  // how many iterations with no change
  u64 noChangeCount = 0;

  // time stuff
  time_t rawtime;
  struct tm * ptm;
  ostringstream timetag;
  time ( &rawtime );
  ptm = gmtime ( &rawtime );
  timetag << ptm->tm_year+1900 << "_" << ptm->tm_mon+1 << "_" << ptm->tm_mday << "_" << ptm->tm_hour << "_" << ptm->tm_min;

  // set default variables
  releasedPath = "diff_release_" + timetag.str() + ".txt";
  statePath = "diff_state_" + timetag.str() + ".txt";
  iterationCount = 100;

	
	frameNum = 6;
	
  // return something if --help passed?
  int parsed = parse_args(argc, argv);

  if(parsed) {
    // print help message and return
    //  return 0;
  }

  
  start_graphics();
  
  // finish setting up variables
  pd = 0.1;
  pe = 1.0 - pd - pp;

  FILE* releasedOut = fopen(releasedPath.c_str(), "w");
  if (releasedOut == NULL) {
    printf("error opening release curve output file, exiting!\n");
    return 1;
  }

  FILE* stateOut;
  if (statePeriod > 0) {
    stateOut = fopen(statePath.c_str(), "w");
    if (stateOut == NULL) {
      printf("error opening state output file, exiting!\n");
      return 1;
    }
  }

  start_graphics();

  CellModel model(
		  n,      // cube width,
		  h,     // cylinder height
		  pd,    // p drug
		  pe,    // p excipient
		  pp,    // p polymer
		  0.001,  // cell size
		  0.001,  // time step
		  7e-6,   // drug diffusion constant
		  7e-6,   // excipient diffusion constant
		  seed     // RNG seed
		  );
	
  mvprintw(0, 0, "cube width %i, pd: %f, pe: %f, pp: %f\n\n", (int)n, pd, pe, pp);
  mvprintw(1, 0, "initializing...");
  refresh();	
  model.setup();
	
  mvprintw(1, 0, "(usage: celldiff count cubeLength pPoly releaseThresh frameInt frameSlice) \n\n");
  mvprintw(2, 0, "cell memory is %d bytes\n\n", model.numCells * sizeof(Cell));
  mvprintw(3, 0, "performing %d iterations on %d cells. press any key to continue...\n\n", iterationCount, n*n*n);
	
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

  fprintf(releasedOut, "0.0");

  while(halt == 0)    {
    step++;
    frameStep++;
    stateStep++;

    if ( step == iterationCount ) {
      halt = HALT_MAX_ITERATIONS;
    }

    if( frameStep == framePeriod ) {
      printFrame(&model, frameNum);
      frameStep = 0;
    }

    if( stateStep == statePeriod ) {
      stateStep = 0;
    }
	  
	  
	  if( (stateStep == 1) && (statePeriod != 0) ) {
		  // print model state data
		  u64 cell;
		  for(cell = 0; cell<model.numCells; cell++) {
			  fprintf(stateOut, "\n%i", model.cells[cell]->state);
			  fprintf(stateOut, "\t%f", model.cells[cell]->concentration[0]);
			  fprintf(stateOut, "\t%f", model.cells[cell]->concentration[1]);
		  }

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

    const double r = released[1] / model.drugMassTotal;
    mvprintw(1, 0, "iteration %d of %d, released %f of %f, ratio %f", step, iterationCount, released[1], model.drugMassTotal, r);

    fprintf(releasedOut, "\n%f", r);
		
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

  fclose(releasedOut);
  if(statePeriod > 0) {
    fclose(stateOut);
  }
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
int parse_args(const int argc, char* const* argv) {

  static struct option long_options[] = {
      {"cubelength", required_argument, 0, 'n'}, 
      {"maxiterations", required_argument, 0, 'c'},
      {"polymer", required_argument, 0, 'p'},
      {"cylinderheight", required_argument, 0, 'h'},
      {"releasedfile", required_argument, 0, 'r'},
      {"statefile", required_argument, 0, 's'},
      {"stateperiod", required_argument, 0, 't'},
      {"nochangecount", required_argument, 0, 'd'},
      {"seed", required_argument, 0, 'e'},
      {"asciiperiod", required_argument, 0, 'a'},
      {0, 0, 0, 0}
    };

  int opt = 0;
  int opt_idx = 0;
  while (1) {
    opt = getopt_long(argc, argv, "n:c:p:h:r:s:t:d:e:a:",
			 long_options, &opt_idx);
    if (opt == -1) { break; }

    switch(opt) {
    case 'n':
      n = atoi(optarg);
      break;
    case 'c':
      iterationCount = atoi(optarg);
      break;
    case 'p' :
      pp = atof(optarg);
      break;
    case 'h' :
      h = atof(optarg);
      break;
    case 'r' :
      releasedPath = optarg;
      break;
    case 's':
     statePath = optarg;
     break;
    case 't':
      statePeriod = atoi(optarg);
      break;
    case 'd':
      noChangeCountThresh = atoi(optarg);
      break;
    case 'e':
      seed = atoi(optarg);
      break;
    case 'a':
      framePeriod = atoi(optarg);
      break;
    default:
      break;
    }
  }
}

// draw an animation frame
void printFrame(CellModel* model, u32 slice) {
  u32 i = slice;
  u32 idx;

  for(u32 j=0; j<n; j++) {
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
