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
#endif 

//------- dissolution constants
// MATLAB and the paper characterize this as some function of neighbors-with-polymer...
// but in the MATLAB code it is constant
// (is it supposed to be == number of neighbors ?)
#define DISS_DENOM 6
// FIXME: this gradual dissolution is described in the paper,
// but in the MATLAB program it seems to always occur in a single iteration.
// so i'm making an arbitrary assumption here:
// number of dissolution steps
// #define DISS_STEPS 5
// amount by which drug/excip concentration increases per dissolution step
// #define DISS_INC 0.2

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
  eStateBound     = 7
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
  f32 concentration[2];
  // counter for gradual dissolution
  u16 dissCount;
  // maximum count for gradual dissolution
  u16 dissMax;
  // dissolution increment
  f32 dissInc;

};

class CellModel {
public:
  CellModel(u32 n,
	    f32 pDrug,
	    f32 pEx,
	    f32 pPoly,
	    f32 cellW,
	    f32 dT,
	    f32 ddrug=7e-6,
	    f32 dex=7e-6,
	    u32 seed=47u);
  ~CellModel(void);
  // set initial values, tablet shape, etc
  void setup(void);
  // advance time in the model by one step
  f32 iterate(void);
private:
  // populate neighbor index array for a given cell
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
  // random number generation
  f32 getRand(void);
	
public: // FIXME: many of these could be privatized
	// cell type distribution
  f32 pDrug;
  f32 pEx;
  f32 pPoly;
  // diffusion constants
  f32 dDrug;
  f32 dEx;
  // time to diffuse length of cell
  f32 dt;
  // length of cell
  f32 cellLength;
  // number of cells on each side of space
  u32 cubeLength;
  u32 cubeLength2;
  // number of total cells
  u32 numCells;
  // initial total mass of drug
  f32 drugMassTotal;
  // current mass of drug, including diffused concentrations
  f32 drugMass;
  // a common intermediate multiplier
  f32 dt_l2;
  // diffusion weights given number of polymer neighbors
  static const f32 diffNMul[7];
  // dissolution steps given number of polymer neighbors
  static const f32 dissNSteps[7];
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
  typedef boost::uniform_real<f32> dist_t;
  rng_t   rngEngine;
  dist_t  rngDist;
  // generator (functor) for streams of values
  boost::variate_generator<rng_t, dist_t>* rngGen;
#endif
};

#endif // header guard
