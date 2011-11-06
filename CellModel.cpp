/*
 *  CellModel.cpp
 *  cellDiffusion
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
Cell::Cell(u64 i) {
  idx = i; 
}

Cell::~Cell() {
}

//================================================================
//================================================================
// ===== CellModel

// index <-> coordinate conversion
u64 CellModel::subToIdx(const u64 x,
                        const u64 y, 
                        const u64 z) {
  //  return z*n*n + y*n + x;
  return cubeLength*(z*cubeLength + y) + x;
}

void CellModel::idxToSub(u64 idx, 
                         u64* pX, 
                         u64* pY, 
                         u64* pZ) {
  *pX = idx % cubeLength;
  *pY = (idx / cubeLength) % cubeLength;
  *pZ = (idx / cubeLength2) % cubeLength;;
}


// populate neighbor index array
void CellModel::findNeighbors(Cell* cell) {
  u64 x, y, z;
  idxToSub(cell->idx, &x, &y, &z);
  // boundary cells are treated differently,
  // and their neighbor indices are irrelevant
  if ( (x==0) || (y==0) || (z==0)
      || (x > (cubeLength-2))
      || (y > (cubeLength-2))
      || (z > (cubeLength-2)) ) {
    for(u8 i=0; i<NUM_NEIGHBORS; i++) {
      cell->neighborIdx[i] = 0;
    }
  } else { 
#if DIAG_NEIGHBORS
    
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
                     u64 n,
                     f32 pdrug,
                     f32 pex,
                     f32 ppoly,
                     f32 cl,
                     f32 dT,
                     f32 ddrug,
                     f32 dex,
                     u32 seed
) :
#if USE_BOOST
rngEngine(), rngDist(0.f, 1.f),
#endif
cubeLength(n),
pDrug(pdrug),
pEx(pex),
pPoly(ppoly),
cellLength(cl),
dt(dT),
dDrug(ddrug),
dEx(dex),
drugMassTotal(0.f),
drugMass(0.f)
{
  cubeLength2 = cubeLength * cubeLength;
  numCells = cubeLength * cubeLength * cubeLength;
  dt_l2 = dt / (cellLength * cellLength);
  cells = new Cell* [numCells];
  cellsUpdate = new Cell* [numCells];
  for(u64 i=0; i<numCells; i++) {
    cells[i] = new Cell(i);
    findNeighbors(cells[i]);
    cellsUpdate[i] = new Cell(i);
  }
  // seed the random number engine
#if USE_BOOST
  rngEngine.seed(seed);
  rngGen = new boost::variate_generator<rng_t, dist_t> (rngEngine, rngDist);
#endif
}

//------ d-tor
CellModel::~CellModel() {
  for(u64 i=0; i<numCells; i++) {
    delete(cells[i]);
    delete(cellsUpdate[i]);
  }
  delete cells;
  delete cellsUpdate;
#if USE_BOOST
  delete rngGen;
#endif
}

//------- set up tablet shape
void CellModel::setup(void) {
  u64 x, y, z;
  u64 cX;
  u64 cY;
  u64 r2;
  f32 rnd;
  u64 cubeR2;
  
	cX = cubeLength >> 1;
	cY = cX;
  cubeR2 = cX * cX;
  
  for (u64 i=0; i<numCells; i++) {
    idxToSub(i, &x, &y, &z);
    //    printf("%d, %d, %d, %d\n", i, x, y, z);
    if ((x < 1) 
        || (y < 1) 
        || (z < 1)
        || (x > (cubeLength-2))
        || (y > (cubeLength-2))
        || (z > (cubeLength-2)) ) {
      // boundary cells
      cells[i]->state = eStateBound;
      cells[i]->concentration[0] = 0.f;
      cells[i]->concentration[1] = 0.f;
    }
    else {
      r2 = (x-cX)*(x-cX) + (y-cY)*(y-cY);
      if(r2 < cubeR2) {
        
        rnd = getRand();
        if(rnd < pDrug) {
          cells[i]->state = eStateDrug;
          drugMassTotal += 1.f;
        } else if (rnd < (pDrug + pEx)) {
          cells[i]->state = eStateEx;
        } else if (rnd < (pDrug + pEx + pPoly)) {
          cells[i]->state = eStatePoly;
      	} 
        
      } else {
        cells[i]->state = eStateBound;
        cells[i]->concentration[0] = 0.f;
        cells[i]->concentration[1] = 0.f;
      }
    }
  }
  
  // create copy for new data
  for(u64 i=0; i<numCells; i++) {
    *(cellsUpdate[i]) = *(cells[i]);
  }
  
  // TODO: tablet compression
  
  drugMass = drugMassTotal;
}

//------- dissolve
eCellState CellModel::dissolve(const Cell* const cell) {
  u8 nw = 0;      // number of wet neighbors
  f32 sumC = 0.f; // sum of neighbor concentrations
  
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
    if (cells[cell->neighborIdx[i]]->state == cell->state) {
      sumC += cells[cell->neighborIdx[i]]->concentration[cell->state];
    }
  }      
  
  // dissolve randomly
  if(getRand() < ((nw - sumC) / DISS_DENOM)) {
    if (cell->state == eStateDrug) {
      cellsUpdate[cell->idx]->state = eStateDissDrug;
    }
    if (cell->state == eStateEx) {
      cellsUpdate[cell->idx]->state = eStateDissEx;
    }
    cellsUpdate[cell->idx]->dissCount = 0;
  }
  return cellsUpdate[cell->idx]->state;
}


// continue dissolution for partially-wetted cells
eCellState CellModel::continueDissolve(const Cell* const cell) {
  cellsUpdate[cell->idx]->dissCount++;
  // FIXME: (?) careful, this concentration index is a nasty enum hack
  cellsUpdate[cell->idx]->concentration[cell->state - 2] += DISS_INC;
  if(cellsUpdate[cell->idx]->dissCount == DISS_STEPS) {
    cellsUpdate[cell->idx]->state = eStateWet;
  }
  return cellsUpdate[cell->idx]->state;
}

// calculate diffusion for fully-dissolved cells
void CellModel::diffuse(const Cell* const cell) {
  f32 cMeanDrug = 0.f;
  f32 cMeanEx = 0.f;
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
  
  /*
   cellsUpdate[cell->idx]->concentration[eStateDrug] = cell->concentration[eStateDrug]
   + (dDrug * nw / cellLength2 * (cMeanDrug - cell->concentration[eStateDrug]) * dt);
   
   cellsUpdate[cell->idx]->concentration[eStateEx] = cell->concentration[eStateEx]
   + (dEx * nw / cellLength2 * (cMeanEx - cell->concentration[eStateEx]) * dt);
   */ // refactored:
  
  const f32 tmp = nw * dt_l2;
  const f32 drugDiff = (dDrug * tmp * (cMeanDrug - cell->concentration[eStateDrug]));
 drugMass += drugDiff;
  
  cellsUpdate[cell->idx]->concentration[eStateDrug] = cell->concentration[eStateDrug] + drugDiff;
  
  cellsUpdate[cell->idx]->concentration[eStateEx] = cell->concentration[eStateEx]
  + (dEx * tmp * (cMeanEx - cell->concentration[eStateEx]));
}

//---------- iterate!!
f32 CellModel::iterate(void) {
  eCellState newState;
  
  /////// TODO: incorporate threading engine. debugging single-threaded version first...
  
  for (u64 i=0; i<numCells; i++) {
    switch(cells[i]->state) {
      case eStatePoly:
        continue;
        break;
      case eStateVoid:
        /// FIXME: what actually is supposed to happen to empty cells???
        continue;
        break;
      case eStateEx:
      case eStateDrug:
        dissolve(cells[i]);
        break;
      case eStateDissDrug:
      case eStateDissEx:
        continueDissolve(cells[i]);
        break;
      case eStateWet:
        diffuse(cells[i]);
        break;
      default:
        break;
    }
  }
    
  ///// TODO: join udpate threads here
  
  // update the cell data
  for(u64 i=0; i<numCells; i++) {
    *(cells[i]) = *(cellsUpdate[i]);
  }
  
  ///// TODO: join copy threads here
  return (drugMassTotal - drugMass) / drugMassTotal;
}

/// random number generation
f32 CellModel::getRand(void) {
#if USE_BOOST
  return (*rngGen)(); 
#else
  return (float)rand() / (float)RAND_MAX;
#endif
}
