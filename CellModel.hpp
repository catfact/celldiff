/*
 *  CellModel.h
 *  cellDiffusion
 *
 *  Created by Ezra Buchla on 10/6/11.
 *
 */

#ifndef _CELLDIFF_CELLMODEL_H
#define _CELLDIFF_CELLMODEL_H_

#define USE_BOOST 0

#if USE_BOOST
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#else
#include <cstdlib>
#endif

#include "types.h"

//======= defines

//------- neighbors
// use diagonal neighbors?
#define DIAG_NEIGHBORS 0
// neighbors per cell
#if DIAG_NEIGHBORS
#define NUM_NEIGHBORS 26
#else
#define NUM_NEIGHBORS 6
#define NUM_NEIGHBORS_R 0.16666666666666666
#endif 

//======= types
// enumeration of cell states
enum eCellState {
  eStateDrug      = 0,    // DON'T CHANGE THESE FIRST FOUR VALUES (stupid hack reasons)
  eStateEx        = 1,
  eStateDissDrug  = 2,   // dissolving drug
  eStateDissEx    = 3,   // dissolving excipient
  eStateWet       = 4,
  eStatePoly      = 5,
  eStateVoid      = 6,
  eStateBound     = 7,
  eStateDummy
};

//======= classes
class Cell {
public:
  // c-tor
  Cell(u32 idx); 
  // d-tor
  ~Cell();
public:
  // state
  eCellState state; 
  // index
  u32 idx;
  // array of neighbor indices
  u32 neighborIdx[NUM_NEIGHBORS];
  // concentration of drug, excipient
  f64 concentration[2];
  // counter for gradual dissolution
  u16 dissCount;
  // maximum count for gradual dissolution
  u16 dissSteps;
  // dissolution increment
  f64 dissInc;
  // diffusion multiplier
  f64 diffMul;
	// dissolution probability for this cell (function of NPN)
	f64 dissProb;
  
};

class CellModel {
public:
  CellModel(
            u32 n,
            f64 h,
            f64 pDrug,
            f64 pPoly,
            f64 cellW,
            f64 ddrug=7e-6,
            f64 dex=7e-6,
            u32 seed=47u,
            f64 dissprobdrug=1.0,
            f64 dissprobex=1.0,
            u32 shellWidth=1,
            f64 polyshellbalance = 1.0,
            f64 bounddiffrate = 0.02,
            f64 dissratescale=1.0,
	    u8 compress=1
            );
  ~CellModel(void);
  // set initial values, tablet shape, etc
  void setup(void);
  // advance time in the model by one step
  f64 iterate(void);
private:
  ///// more setup funxtions...
  // initial distribution of particles
  void distribute(void);
  // compression step
  void compress(void);
  // find cells that need processing
  void findCellsToProcess(void);
  // paopulate neighbor index array for a given cell
  void findNeighbors(Cell* cell);
  // decide whether to dissolve given cell; return new state
  eCellState dissolve(const Cell* const cell);
  // continue dissolving this cell; return new states
  eCellState continueDissolve(const Cell* const cell);
  // calculate diffusion on this cell
  void diffuse(const Cell* const cell);
  // calculate the current mass of drug remaining 
  // FIXME: this is rather inefficient
  void calcDrugMass(void);
  // index / coordinates conversion
  u32 subToIdx(const u32 x, const u32 y, const u32 z);
  void idxToSub(u32 idx, u32* pX, u32* pY, u32* pZ);
  // set state of a 2x2x2 block of cells
  void setBlockState(const u32 idx, eCellState state);
  // set state of a single cell
  void setCellState(const u32 idx, eCellState state);
  // random number generation
  f64 getRand(void);
	
public: // FIXME: many of these could be privatized
	// cell type distribution
  //  u32 nDrug;
  // u32 nEx;
  // u32 nPoly;
  f64 pPoly;
  f64 pDrug;
  // diffusion constants
  f64 dDrug;
  f64 dEx;
  // time to diffuse length of cell
  f64 dt;
  // length of cell
  f64 cellLength;
  // number of cells on each side of space
  u32 cubeLength;
  u32 cubeLength2;
  // shell width
  u32 shellN;
  // cylinder height as fraction of cube height
  f64 cylinderHeight;
  // number of total cells
  u32 numCells;
  // initial total mass of drug
  f64 drugMassTotal;
  // current mass of drug, including diffused concentrations
  f64 drugMass;
  // trapped drug mass (never changes)
  f64 trappedDrugMass;
	// dissolution probability (drug)
	f64 dissProbDrug;
	// dissolution probability (excipient)
	f64 dissProbEx;
  // boundary diffusion rate
  f64 boundDiff;
  // shell width in cells
  u32 wShell;
  // polymer-shell "imbalance factor"
  f64 pShellBalance;
  // dissolution rate scaling factor for all cells
  f64 dissratescale;
  // compressino flag
  u8 compressFlag;
  // diffusion weights given number of polymer neighbors
  static const f64 diffNMul[7];
  // dissolution steps given number of polymer neighbors
  static const f64 dissNSteps[7];
  // flattened array of all cells
  Cell**        cells;         
  // copy for updating after iteration
  Cell**        cellsUpdate;
  // cells-to-process (drug, excip, water, diffusing, or immediate boundary) 
  u32* cellsToProcess;
  u32 numCellsToProcess;
	
  //====== random number stuff
#if USE_BOOST
  // randomization algorithm
  typedef boost::mt19937 rng_t; // mersenne twister
  // distribution (maps algorithm to data type and range)
  typedef boost::uniform_real<f64> dist_t;
  rng_t   rngEngine;
  dist_t  rngDist;
  // generator (functor) for streams of values
  boost::variate_generator<rng_t, dist_t>* rngGen;
#endif
};

#endif // header guard
