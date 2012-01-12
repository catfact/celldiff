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
Cell::Cell(u32 i) {
  idx = i;
  concentration[0] = 0.f;
  concentration[1] = 0.f;
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
  //  1, 50, 100, 100, 200, 400, 800
  1, 1, 1, 1, 1, 1, 1
	// 50, 50, 50, 50, 50, 50, 50
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
  *pZ = (idx / cubeLength2) % cubeLength;;
}


// populate neighbor index array
void CellModel::findNeighbors(Cell* cell) {
  u32 x, y, z;
  idxToSub(cell->idx, &x, &y, &z);
  // boundary cells; their neighbor indices are irrelevant
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
                     u32 n,
		     f64 h,
                     f64 pdrug,
                     f64 pex,
                     f64 ppoly,
                     f64 cl,
                     f64 dT,
                     f64 ddrug,
                     f64 dex,
                     u32 seed
		     ) :
#if USE_BOOST
  rngEngine(), rngDist(0.f, 1.f),
#endif
  cubeLength(n),
  cylinderHeight(h),
  pDrug(pdrug),
  pEx(pex),
  pPoly(ppoly),
  cellLength(cl),
  dt(dT),
  dDrug(ddrug),
  dEx(dex),
  drugMassTotal(0.0),
  trappedDrugMass(0.0),
  drugMass(0.0),
  numCellsToProcess(0)
{
  cubeLength2 = cubeLength * cubeLength;
  numCells = cubeLength * cubeLength * cubeLength;
  dt_l2 = dt / (cellLength * cellLength);
  cells =				new Cell* [numCells];
  cellsUpdate =		new Cell* [numCells];
  cellsToProcess =	new u32 [numCells];
  for(u32 i=0; i<numCells; i++) {
    cells[i] =			new Cell(i);
    cellsUpdate[i] =	new Cell(i);
  }
  // seed the random number engine
#if USE_BOOST
  rngEngine.seed(seed);
  rngGen = new boost::variate_generator<rng_t, dist_t> (rngEngine, rngDist);
#endif
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

//------- set up tablet shape, find neighbors, and populate
void CellModel::setup(void) {
  u32 cX;
  u32 cY;
  u32 r2;
  f64 rnd;
  u32 cubeR2;
  u32 idx;
  u32 nIdx;
  u32 nIdx2;
  u8 l, m, n;
  u32 numH;

  cX = cubeLength >> 1;
  cY = cX;
  cubeR2 = (cX-1) * (cX-1);
  eCellState swapstate;

  numH = cylinderHeight * cubeLength;
	
  // offset coordinates to get diagonals in 2x2x2
  const u8 diags[4][3]	= { {0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0} };
  // complementary diagonals
  const u8 diagsNot[4][3]	= { {1, 1, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  u8 diag = 0;
	
  //-------- with compression
	
  // loop over alternate [x, y, z] indices
  // process each 2x2x2 based on state result
  for(u32 i=0; i<cubeLength; i += 2) {
    for(u32 j=0; j<cubeLength; j += 2) {
      for(u32 k=0; k<cubeLength; k += 2) {
	idx = subToIdx(i, j, k);
	if ((i < 1) 
	    || (j < 1) 
	    || (k < 1)
	    || (i > (cubeLength-3))
	    || (j > (cubeLength-3))
	    || (k > (cubeLength-3)) ) {
	  // boundary cells
	  // fill the whole 2x2x2
	  for(l=0; l<2; l++) {
	    for(m=0; m<2; m++) {
	      for(n=0; n<2; n++) {
		nIdx = subToIdx(i+l, j+m, k+n);
		cells[nIdx]->state = eStateBound;
	      }
	    }
	  }
	} // end edge branch
	else {
	  r2 = (i-cX+1)*(i-cX+1) + (j-cY+1)*(j-cY+1);
	  if((r2 < cubeR2) && (k <= (numH))) {
						
	    rnd = getRand();
	    if(rnd < pDrug) {
	      // drug cells: fill diagonals
	      for(diag = 0; diag<4; diag++) {
		nIdx = subToIdx(i+diags[diag][0], j+diags[diag][1], k+diags[diag][2]);
		cells[nIdx]->state = eStateDrug;
		drugMassTotal += 1.f;
	      }
	      // opposite diagonals are empty
	      for(diag = 0; diag<4; diag++) {
		nIdx = subToIdx(i+diagsNot[diag][0], j+diagsNot[diag][1], k+diagsNot[diag][2]);
		cells[nIdx]->state = eStateVoid;
	      }
							
	    } else if (rnd < (pDrug + pEx)) {
	      // excipient cells: fill diagonals
	      for(diag = 0; diag<4; diag++) {
		u32	nIdx = subToIdx(i+diags[diag][0], j+diags[diag][1], k+diags[diag][2]);
		cells[nIdx]->state = eStateEx;
	      }
	      // opposite diagonals are empty
	      for(diag = 0; diag<4; diag++) {
		nIdx = subToIdx(i+diagsNot[diag][0], j+diagsNot[diag][1], k+diagsNot[diag][2]);
		cells[nIdx]->state = eStateVoid;
	      }
	    } else if (rnd < (pDrug + pEx + pPoly)) {
	      cells[idx]->state = eStatePoly;
	      // poly cells: fill the whole 2x2x2
	      for(l=0; l<2; l++) {
		for(m=0; m<2; m++) {
		  for(n=0; n<2; n++) {
		    nIdx = subToIdx(i+l, j+m, k+n);
		    cells[nIdx]->state = eStatePoly;
		  }
		}
	      }
	    } 
						
	  } // end in-cylinder branch
	  else { 
	    // boundary cells
	    // fill the whole 2x2x2
	    for(l=0; l<2; l++) {
	      for(m=0; m<2; m++) {
		for(n=0; n<2; n++) {
		  nIdx = subToIdx(i+l, j+m, k+n);
		  cells[nIdx]->state = eStateBound;
		}
	      }
	    }
	  }
	} // end not-edge branch
      } // end k loop
    } // end j loop
  } // end i loop
	
  /// 2nd loop: each polymer 2x2x2 randomly swaps diagonal complements with neighbors
  for(u32 i=0; i<cubeLength; i += 2) {
    for(u32 j=0; j<cubeLength; j += 2) {
      for(u32 k=0; k<cubeLength; k += 2) {
	idx = subToIdx(i, j, k);
	if( (cells[idx]->state == eStatePoly) && (cells[idx+1]->state == eStatePoly) ) {
	  u32 nIdxBase[3] = {i, j, k};
	  // only want to swap with neighbor meta-cells that are not poly...
	  // array for storing neighbor idx's which meet this criterion.
	  u64 notPolyN[6] = { 0, 0, 0, 0, 0, 0 };
	  u8 numNotPolyN = 0;
	  u8 swapN;

	  for (u8 n=0; n < NUM_NEIGHBORS; n++) {
	    switch(n) { 
	    case 0:
	      nIdxBase[0] = i+2;
	      nIdxBase[1] = j;
	      nIdxBase[2] = k;
	      break;
	    case 1:
	      nIdxBase[0] = i-2;
	      nIdxBase[1] = j;
	      nIdxBase[2] = k;
	      break;
	    case 2:
	      nIdxBase[0] = i;
	      nIdxBase[1] = j+2;
	      nIdxBase[2] = k;
	      break;
	    case 3:
	      nIdxBase[0] = i;
	      nIdxBase[1] = j-2;
	      nIdxBase[2] = k;
	      break;
	    case 4:
	      nIdxBase[0] = i;
	      nIdxBase[1] = j;
	      nIdxBase[2] = k+2;
	      break;
	    case 5:
	      nIdxBase[0] = i;
	      nIdxBase[1] = j;
	      nIdxBase[2] = k-2;
	      break;
	    default:
	      nIdxBase[0] = i;
	      nIdxBase[1] = j;
	      nIdxBase[2] = k;
	      break;
	    }
	    nIdx = subToIdx(nIdxBase[0], nIdxBase[1], nIdxBase[2]);
	    if( (cells[nIdx]->state != eStatePoly)
		&& (cells[nIdx+1]->state != eStatePoly)
		&& (cells[nIdx+1]->state != eStateBound)) {
	      // this neighbor is not polymer, and has not been swapped, so add it to the swappable list
	      notPolyN[numNotPolyN] = nIdx;
	      numNotPolyN++;
	    }
	  }
	  
	  // continue loop if there are no non-poly neighors,
	  // otherwise randomly choose between them and swap diagonals
	  if (numNotPolyN == 0) { continue; } 
	  else {
	    swapN = (u8)(getRand() * ((f64)numNotPolyN  - 0.5f));
	    idxToSub(notPolyN[swapN], &(nIdxBase[0]), &(nIdxBase[1]), &(nIdxBase[2]));
	    for(diag = 0; diag<4; diag++) {
	      nIdx = subToIdx(i+diags[diag][0], j+diags[diag][1], k+diags[diag][2]);
	      nIdx2 = subToIdx(nIdxBase[0]+diagsNot[diag][0], nIdxBase[1]+diagsNot[diag][1], nIdxBase[2]+diagsNot[diag][2]);
	      swapstate = cells[nIdx2]->state;
	      if (swapstate != eStateBound) {
		cells[nIdx2]->state = cells[nIdx]->state;
		cells[nIdx]->state = swapstate;
	      }
	    }
	  }
	}
      }
    }
  }
	
	
  /// 3rd loop: find cells to process
  u8 proc = 0;
  eCellState tmpState;
  for(u32 i=0; i<numCells; i++) {
    findNeighbors(cells[i]);
    /// add to processCells
    switch (cells[i]->state) {
    case eStateDrug:
    case eStateEx:
    case eStateVoid:
      proc = 1;
      // printf("{ %d, %d } ; ", (int)numCellsToProcess, (int)cells[i]->idx);
      break;
    case eStateBound:
      // want to process boundary cells only if they adjoin a non-boundary
      for(u8 ni = 0; ni<NUM_NEIGHBORS; ni++) {
	tmpState = cells[cells[i]->neighborIdx[ni]]->state;
	proc |= (tmpState != eStatePoly && tmpState != eStateBound);
      }
    break;
    default:
      break;
    }
    if(proc) { 	
      // find neighbors-with-polymer count
      u8 np = 0;
      for(u8 nb=0; nb<NUM_NEIGHBORS; nb++) {
        if ( cells[cells[i]->neighborIdx[nb]]->state == eStatePoly ) {
          np++;
        }
      }
      
      // don't need to process if cell is trapped by polymer
      if (np > 5) {
        if(cells[i]->state == eStateDrug) {
          trappedDrugMass += 1.0;
        }
        continue;
      }
      
      cellsToProcess[numCellsToProcess] = cells[i]->idx;
      numCellsToProcess++;
      
      cells[i]->diffMul = diffNMul[np];
      cells[i]->dissSteps = dissNSteps[np];
      cells[i]->dissInc = 1.0 / (f64)(cells[i]->dissSteps);
      ///// DEBUG
      //    int dum=0;
      // dum++;
    }
  }
    
  // initialize the update data
  for(u32 i=0; i<numCells; i++) {
    *(cellsUpdate[i]) = *(cells[i]);
  }
  drugMass = drugMassTotal;
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
  if (getRand() < ((nw - sumC) / (float)DISS_DENOM)) {
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
  /////// TODO: incorporate threading engine. debugging single-tvhreaded version first...
	
  for (u32 i=0; i<numCellsToProcess; i++) {
    cell = cells[cellsToProcess[i]];
    // FIXME: processing order of these cases could be optimized
    switch(cell->state) {
      /*
    case eStatePoly:
      // shouldn't get here!
      // polymer cells: no change
      break;
	   */
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
			/*
    case eStateBound:
      cell->concentration[0] = 0.f;
      cell->concentration[1] = 0.f;
			 */
    default:
      break;
    }
  }
    
  ///// TODO: join udpate threads here
	
  // update the cell data
  // FIXME: memcpy() in this function is eating 20% of CPU time.
  // should be able to just swap pointers somehow...

  // (this is not cutting it)
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
  // better to update during the diffusion step?
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
      drugMass += 1.f - (cell->dissInc * cell->dissCount) + cell->concentration[eStateDrug];
      // drugMass += 1.0;
      break;
    case eStateWet:
    case eStateDissEx:
      drugMass += cell->concentration[eStateDrug];
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
