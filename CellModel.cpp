/*
 *  CellModel.cpp
 *  celldiff
 *
 *  Created by Ezra Buchla on 10/6/11.
 */

#if USE_BOOST
#include <boost/thread.hpp>
#endif

#include <cstdio>
#include <cassert>
#include "CellModel.hpp"

//================================================================
//================================================================
//======= Cell
Cell::Cell(u32 i) {
  state = eStateDummy;
  idx = i;
  concentration[0] = 0.0;
  concentration[1] = 0.0;
}

Cell::~Cell() {
}

//================================================================
//================================================================
// ===== CellModel


///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
///!!!!!!!!!!! tweakable: 
//// constant diffusion multiplier given count of poly neighbors
const f64 CellModel::diffNMul[7] = {
  // 1.0, 0.9, 0.1, 0.01, 0.001, 0.0, 0.0
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
};

//// constant dissolution steps given count of poly neighbors
const f64 CellModel::dissNSteps[7] = {	
	10, 10, 10, 10, 10, 10, 10
};
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// index <-> coordinate conversion
u32 CellModel::subToIdx(const u32 x,
                        const u32 y, 
                        const u32 z) {
  return cubeLength*(z*cubeLength + y) + x;
}

void CellModel::idxToSub(u32 idx, 
                         u32* pX, 
                         u32* pY, 
                         u32* pZ) {
  *pX = idx % cubeLength;
  *pY = (idx / cubeLength) % cubeLength;
  *pZ = (idx / cubeLength2) % cubeLength;
}


// populate neighbor index array
void CellModel::findNeighbors(Cell* cell) {
  u32 x, y, z;
  idxToSub(cell->idx, &x, &y, &z);
  if ( (x==0) || (y==0) || (z==0)
      || (x > (cubeLength-2))
      || (y > (cubeLength-2))
      || (z > (cubeLength-2)) ) {
    // these are cells on the extreme boundary of the cube...
    // hopefully their neighbor indices should not be used
    for(u8 i=0; i<NUM_NEIGHBORS; i++) {
      cell->neighborIdx[i] = 0;
    }
  } else { 
#if DIAG_NEIGHBORS
		//...
#else
    cell->neighborIdx[0] = subToIdx(x+1,  y,    z);
    cell->neighborIdx[1] = subToIdx(x-1,  y,    z);
    cell->neighborIdx[2] = subToIdx(x,    y+1,  z);
    cell->neighborIdx[3] = subToIdx(x,    y-1,  z);
    cell->neighborIdx[4] = subToIdx(x,    y,    z+1);
    cell->neighborIdx[5] = subToIdx(x,    y,    z-1);
  }
#endif
}

//------ c-tor
CellModel::CellModel(
                     u32 n,
                     f64 h,
                     f64 pdrug,
                     f64 ppoly,
                     f64 cl,
                     f64 dT,
                     f64 ddrug,
                     f64 dex,
                     u32 seed,
                     f64 dprob,
                     f64 dpolyscale,
                     u32 shellwidth,
                     f64 polyshellbalance,
                     f64 bounddiffrate,
                     f64 dissScale
                     ) :
#if USE_BOOST
rngEngine(), rngDist(0.f, 1.f),
#endif
cubeLength(n),
cylinderHeight(h),
cellLength(cl),
dt(dT),
dDrug(ddrug),
dEx(dex),
drugMassTotal(0.0),
trappedDrugMass(0.0),
drugMass(0.0),
numCellsToProcess(0),
wShell(shellwidth),
pPoly(ppoly),
pDrug(pdrug),
pShellBalance(polyshellbalance),
boundDiff(bounddiffrate),
dissratescale(dissScale)
{
  
  cubeLength2 = cubeLength * cubeLength;
  numCells = cubeLength * cubeLength * cubeLength;
  dt_l2 = dt / (cellLength * cellLength);
  
  cells =		new Cell* [numCells];
  cellsUpdate =		new Cell* [numCells];
  cellsToProcess =	new u32 [numCells];
  
  for(u32 i=0; i<numCells; i++) {
    cells[i] =	        new Cell(i);
    cellsUpdate[i] =	new Cell(i);
  }
  // seed the random number engine
#if USE_BOOST
  rngEngine.seed(seed);
  rngGen = new boost::variate_generator<rng_t, dist_t> (rngEngine, rngDist);
#else
  srand((u32)seed);
#endif
	dissprob = dprob / NUM_NEIGHBORS;
	disspolyscale = dpolyscale;
}

//------ d-tor
CellModel::~CellModel() {
  for(u32 i=0; i<numCells; i++) {
    delete(cells[i]);
    delete(cellsUpdate[i]);
  }
  delete cells;
  delete cellsUpdate;
  delete cellsToProcess;
#if USE_BOOST
  delete rngGen;
#endif
}

//------- dissolve
eCellState CellModel::dissolve(const Cell* const cell) {
  u8 nw = 0;      // number of wet neighbors
  u8 np = 0;      // number of polymer neighbors 
  f64 sumC = 0.f; // sum of neighbor concentrations
  
	// count the wet/boundary neighbors
  for(u8 i = 0; i < NUM_NEIGHBORS; i++) {
    if ((cells[cell->neighborIdx[i]]->state == eStateWet) || (cells[cell->neighborIdx[i]]->state == eStateBound)) {
      nw++;
    }
  }
  // return early if there are no wet neighbors
  if (nw == 0) {
    return cellsUpdate[cell->idx]->state;
  }
	
  // compare dry-neighbor states with this cell's state
  for(u8 i = 0; i < NUM_NEIGHBORS; i++) {
    if (cells[cell->neighborIdx[i]]->state == eStateWet) {
      sumC += cells[cell->neighborIdx[i]]->concentration[cell->state];
    }
  }
  
  // dissolve randomly
  //  if (getRand() < ((nw - sumC) / (float)DISS_DENOM)) {
  if (getRand() < ((nw - sumC) * cell->dissProb)) {
    if (cell->state == eStateDrug) {
      cellsUpdate[cell->idx]->state = eStateDissDrug;
      cellsUpdate[cell->idx]->dissCount = 0;
    }
    if (cell->state == eStateEx) {
      cellsUpdate[cell->idx]->state = eStateDissEx;
      cellsUpdate[cell->idx]->dissCount = 0;
    }
    if (cell->state == eStateVoid) {
      // FIXME: should void cells "dissolve" gradually?
      cellsUpdate[cell->idx]->state = eStateWet;
    }		
  }
  return cellsUpdate[cell->idx]->state;
}


// continue dissolution for partially-wetted cells
eCellState CellModel::continueDissolve(const Cell* const cell) {
  cellsUpdate[cell->idx]->dissCount++;
  // FIXME: (?) careful, this concentration index is a nasty enum hack
  cellsUpdate[cell->idx]->concentration[cell->state - 2] = cells[cell->idx]->concentration[cell->state - 2] + cell->dissInc;
  if(cellsUpdate[cell->idx]->dissCount >= cell->dissSteps) {
    cellsUpdate[cell->idx]->state = eStateWet;
  }
  return cellsUpdate[cell->idx]->state;
}


// calculate diffusion for fully-dissolved cells
void CellModel::diffuse(const Cell* const cell) {
  f64 cMeanDrug = 0.f;
  f64 cMeanEx = 0.f;
  u8 nw = 0;
	
  for(u8 i=0; i<NUM_NEIGHBORS; i++) {
    if ((cells[cell->neighborIdx[i]]->state == eStateWet) || (cells[cell->neighborIdx[i]]->state == eStateBound)) {
      nw++;
      cMeanDrug += cells[cell->neighborIdx[i]]->concentration[eStateDrug];
      cMeanEx += cells[cell->neighborIdx[i]]->concentration[eStateEx];
    }
  }
  // no wet neighbors => no effect
  if (nw == 0) { return; }
	
  cMeanDrug /= nw;
  cMeanEx /= nw;
	
  
  //   cellsUpdate[cell->idx]->concentration[eStateDrug] = cell->concentration[eStateDrug]
  //   + (dDrug * nw / cellLength2 * (cMeanDrug - cell->concentration[eStateDrug]) * dt);
  
  //    cellsUpdate[cell->idx]->concentration[eStateEx] = cell->concentration[eStateEx]
  //   + (dEx * nw / cellLength2 * (cMeanEx - cell->concentration[eStateEx]) * dt);
  
  // refactored:
  const f64 tmp = nw * dt_l2;
  const f64 drugDiff = (dDrug * tmp * (cMeanDrug - cell->concentration[eStateDrug] * cell->diffMul));  
	
  cellsUpdate[cell->idx]->concentration[eStateDrug] = cell->concentration[eStateDrug] + drugDiff;
	
  cellsUpdate[cell->idx]->concentration[eStateEx] = cell->concentration[eStateEx]
  + (dEx * tmp * (cMeanEx - cell->concentration[eStateEx]  * cell->diffMul));
}




//---------- iterate!!
f64 CellModel::iterate(void) {
  Cell* cell;
  ///////// DEBUG
  u8 dum=0;
  //////////
  
  /////// TODO: incorporate threading engine. debugging single-threaded version first...
	
  for (u32 i=0; i<numCellsToProcess; i++) {
    cell = cells[cellsToProcess[i]];
    // FIXME: processing order of these cases could be optimized
    switch(cell->state) {
        case eStatePoly:
         // shouldn't get here!
         // polymer cells: no change
        dum++;
         break;
      case eStateWet:
        diffuse(cell);
        break;
      case eStateVoid:
        // void cells: deissolve (FIXME?)
        dissolve(cell);
        break;
      case eStateEx:
      case eStateDrug:
        // drug or excipient: 
        dissolve(cell);
        break;
      case eStateDissDrug:
      case eStateDissEx:
        continueDissolve(cell);
        break;
      case eStateBound:
        diffuse(cell);
	//// TODO: optimize
        cell->concentration[0] *= boundDiff; // exponential decay
        cell->concentration[1] *= boundDiff;
        // denormal and saturate low
        if(cell->concentration[0] < 0.000000000001) cell->concentration[0] = 0.0;
        if(cell->concentration[1] < 0.000000000001) cell->concentration[1] = 0.0;
        break;
      default:
        break;
    }
  }
  
  ///// TODO: join udpate threads here
	
  // update the cell data
  // FIXME: memcpy() in this function is eating 20% of CPU time.
  // should be able to just swap pointers... 
  // (why doesn't this work??)
  /*
   Cell** cellsTmp = cells;
   cells = cellsUpdate;
   cellsUpdate = cellsTmp; 
  */
  
  for(u32 i=0; i<numCellsToProcess; i++) {
    *(cells[cellsToProcess[i]]) = *(cellsUpdate[cellsToProcess[i]]);
  }
  
	
  ///// TODO: join copy threads here
	
  calcDrugMass();
  
  return drugMassTotal - drugMass;
  //  return drugMass;
}

void CellModel::calcDrugMass(void) {
  // calculate current drug mass
  // FIXME: this is the slow way to do it.
  // better to update during the diffusion step, and save a loop
  // tryingfor accuracy first.
  Cell* cell;	
  drugMass = trappedDrugMass;
	
  for (u32 i=0; i<numCellsToProcess; i++) {
    cell = cells[cellsToProcess[i]];
    switch(cell->state) {
      case eStateDrug:
        drugMass += 1.0;
        break;
      case eStateDissDrug:
        // drugMass += 1.f - (cell->dissInc * cell->dissCount) + cell->concentration[eStateDrug];
        drugMass += 1.0;
        break;
      case eStateWet:
        drugMass += cell->concentration[eStateDrug];
        break;
      case eStateDissEx:
        break;
      default:
        break;
    }
  }
}

/// random number generation
f64 CellModel::getRand(void) {
#if USE_BOOST
  return (*rngGen)(); 
#else
  return (float)rand() / (float)RAND_MAX;
#endif
}
