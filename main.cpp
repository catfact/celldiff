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
#include <cstdarg>

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
static f64 pd=0.1, pp=0.4;
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
// dissolution probability scale (denominator == dissprob * numneighbors)
static f64 dissprob = 1;
// dissolution prob / polymer correlation factor
static f64 disspolyscale = 0.0;
// width of polymer shell
static u32 polyShellWidth = 1;
// polymer-shell "imbalance factor"
static f64 polyShellBalance = 1.0;
// boundary decay factor
static f64 boundDiff = 0.9;
// no-graphics mode
static u8 nographics = 0;
// drug diffusion rate;
static f64 drugdiff = 0.000001;
// excipient diffusion rate;
static f64 exdiff = 0.000001;
// dissolution time scaling;
static f64 dissScale = 1.0;
// cell size (in meters)
static f64 cellsize = 0.001;
//  time step size (in seconds)
static f64 timestep = 0.001;
// ncurses window pointer
static WINDOW* win;

//============== function declarations
int main(const int argc, char* const* argv);

static int parse_args(const int argc, char* const* argv);
static void start_graphics(void);
static void end_graphics(void);
static void printFrame(CellModel* model, u32 frame);
static void print(const int x, const int y, const char* fmt, ...);

//============== function definitions
void print(const int x, const int y, const char* fmt, ...) {
  
  va_list args;
  va_start(args,fmt);
  if(nographics) {
    vprintf(fmt, args);
    printf("\n");
  } else {
    move(x, y);
    vw_printw(win, fmt, args);
    wrefresh(win);
  }
  va_end(args);
}

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
	
  // return something if --help passed?
  int parsed = parse_args(argc, argv);
  
  if(parsed) {
    // print help message and return
    //  return 0;
  }
  
  frameNum = n >> 1; // show center slice
  
  // finish setting up variables
  pd = 0.1;
  
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
  
  if(nographics) {} else { start_graphics(); }
  
  CellModel model(
                  n,      // cube width,
                  h,     // cylinder height
                  pd,    // p drug
                  pp,    // p polymer
                  cellsize,  // cell size
                  timestep,  // time step
                  drugdiff,   // drug diffusion constant
                  exdiff,          // excipient diffusion constant
                  seed,          // RNG seed,
                  dissprob,    // dissolution probability scale,
                  disspolyscale, // dissolution probability -> polymer correlation factor,
                  polyShellWidth,       // shell width
                  polyShellBalance,  // shell balance
                  boundDiff,       // boundary diffusino factor (exponential)
                  dissScale          // dissolution time scaling 
                  );
	
  print(0, 0, "cube width %i, pd: %f, pp: %f", (int)n, pd, pp);
  if(nographics) {} else { print(1, 0, "initializing..."); }
  model.setup();
  
  print(2, 0, "cell memory is %d bytes", model.numCells * sizeof(Cell));
  if(nographics) {
    print(3, 0, "performing %d iterations on %d cells.", iterationCount, n*n*n);
  } else {
    print(3, 0, "performing %d iterations on %d cells. press any key to continue...", iterationCount, n*n*n);
    getchar();
    print(1, 0, "                                                                            ");
    print(2, 0, "                                                                            ");
    print(3, 0, "                                                                            ");
    print(4, 0, "                                                                            ");
	}
  
  int step = 0;
  u32 frameStep = 0;
  u8 halt = 0;
  
  fprintf(releasedOut, "0.0\t0.0");
  
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
    print(1, 0, "iteration %d of %d, released %f of %f, ratio %f", step, iterationCount, released[1], model.drugMassTotal, r);
    fprintf(releasedOut, "\n%f\t%f", model.dt * (float)step, r);
    
  } // end main loop
  
  printFrame(&model, frameNum); 
  
  switch(halt) {
    case HALT_MAX_ITERATIONS:
      if(nographics) {
        print(2, 0, "finished maximum iterations; simulation halted.");
      } else {
        print(2, 0, "finished maximum iterations; simulation halted. press any key to quit...              ");
        getchar();
      }
      break;
    case HALT_NO_CHANGE:
      if(nographics) {
        print(2, 0, "released mass appears stable; simulation halted.");
      } else {
        print(2, 0, "released mass appears stable; simulation halted. press any key to quit...           ");
        getchar();
      }
      break;
  }
  
  
  fclose(releasedOut);
  if(statePeriod > 0) {
    fclose(stateOut);
  }
  if (nographics) { } else { end_graphics(); }
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
  init_pair(9, COLOR_WHITE, COLOR_RED);
  win = newwin(10, 64, 0, 0);
}

//------ graphics de-init
static void end_graphics(void) {
  delwin(win);
  endwin();
}

//------ parse arguments
int parse_args(const int argc, char* const* argv) {
  
  static struct option long_options[] = {
    {"cubelength",        required_argument, 0, 'n'}, 
    {"maxiterations",		  required_argument, 0, 'c'},
    {"polymerratio",      required_argument, 0, 'p'},
    {"drugratio",         required_argument, 0, 'g'},
    {"cylinderheight",	  required_argument, 0, 'h'},
    {"releasedfile",		  required_argument, 0, 'r'},
    {"statefile",         required_argument, 0, 's'},
    {"stateperiod",       required_argument, 0, 't'},
    {"nochangecount",		  required_argument, 0, 'd'},
    {"seed",              required_argument, 0, 'e'},
    {"asciiperiod",       required_argument, 0, 'a'},
    {"dissdenom",         required_argument, 0, 'o'},
    {"disspolyscale",		  required_argument, 0, 'i'},
    {"polyshellwidth",	  required_argument, 0, 'w'},
    {"polyshellbalance",  required_argument, 0, 'b'},
    {"boundarydiffusion", required_argument, 0, 'f'},
    {"nographics",        required_argument, 0, 'x'},
    {"drugdiffusionrate", required_argument, 0, 'u'},
    {"exdiffusionrate",   required_argument, 0, 'k'},
    {"dissolutionrate",   required_argument, 0, 'l'},
    {"cellsize",          required_argument, 0, 'y'},
    {"timestep",          required_argument, 0, 'z'},
    {0, 0, 0, 0}
  };
  
  int opt = 0;
  int opt_idx = 0;
  while (1) {
    opt = getopt_long(argc, argv, "n:c:p:g:h:r:s:t:d:e:a:o:i:w:b:f:x:u:k:l:y:z:",
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
      case 'g' :
        pd = atof(optarg);
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
      case 'o':
        dissprob = atof(optarg);
        break;
      case 'i':
        disspolyscale = atof(optarg);
        break;
      case 'w':
        polyShellWidth = atoi(optarg);
        break;
      case 'b':
        polyShellBalance = atof(optarg);
        break;
      case 'f':
        boundDiff = atof(optarg);
        break;
      case 'x':
        nographics = atoi(optarg);
        break;
      case 'u':
        drugdiff = atof(optarg);
        break;
      case 'k':
        exdiff = atof(optarg);
        break;
      case 'l':
        dissScale = atof(optarg);
        break;
      case 'y':
        cellsize = atof(optarg);
        break;
      case 'z':
        timestep = atof(optarg);
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
      if ((model->cells[idx]->state == eStateWet) || (model->cells[idx]->state == eStateBound) ) {
        mvprintw(6+j, k * 2, "%d0", (int)(model->cells[idx]->concentration[eStateDrug] * 99.0));
      } else {
        // if(model->cells[idx]->state == eStateDummy) {
        //   mvprintw(6+j, k * 2, "XX", model->cells[idx]->state);
        // } else {
        mvprintw(6+j, k * 2, "%d ", model->cells[idx]->state);
        // }
      }
      attroff(COLOR_PAIR(model->cells[idx]->state + 1));
    }
  }
  refresh();
}
