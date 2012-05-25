/*
 *  CellModelSetup.cpp
 *  celldiff
 *
 *  i'm putting CellModel initialization methods in their own file, just for readability
 *
 *  Created by Ezra Buchla on 03/05/2012
 */
#include <vector>
#include <algorithm>
#include "CellModel.hpp"

using namespace std;


//// top-level setup function. initializes cell type data
void CellModel::setup(void) {
  //////////// DEBUG
    
  u32 stateCount[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  u32 cellsInTablet = 0;
  f64 drugCountRatio, polyCountRatio;
  f64 exCountRatio, voidCountRatio;
  u32 cellsAdjoining = 0;
  u8 dum=0;
  
  ////////////
  
  
  ///---- DISTRIBUTE
  this->distribute();
  
  ////////// DEBUG
  /*
  // another loop over all cells to verify final cell distribution count
  for(u8 s=0; s<8; s++) { stateCount[s] = 0; }
  for (u32 n=0; n<numCells; n++) {
  stateCount[cells[n]->state]++;
  }
  // easier-to-read totals
  cellsInTablet   = stateCount[eStateDrug] + stateCount[eStateEx] + stateCount[eStatePoly] + stateCount[eStateVoid];
  drugCountRatio  = (f64)stateCount[eStateDrug] / (f64)cellsInTablet;
  polyCountRatio  = (f64)stateCount[eStatePoly] / (f64)cellsInTablet;
  exCountRatio    = (f64)stateCount[eStateEx] / (f64)cellsInTablet;
  voidCountRatio  = (f64)stateCount[eStateVoid] / (f64)cellsInTablet;
  // debugger hook
  dum++;
  */
  /////////////////
  
  //------- COMPRESS
  this->compress();
	
  // offset coordinates to get diagonals in 2x2x2
  const u8 diags[4][3]	= { {0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0} };
  // complementary diagonals
  const u8 diagsNot[4][3]	= { {1, 1, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  u8 diag = 0;
	
  /// find cells to process
  numCellsToProcess = 0;
  u8 proc = 0;
  eCellState tmpState;
  
  drugMassTotal = 0.0;
  
  for(u32 i=0; i<numCells; i++) {
    findNeighbors(cells[i]);
    proc=0;
    switch (cells[i]->state) {
    case eStatePoly:
      proc = 0;
      break;
    case eStateDrug:
      proc = 1;
      drugMassTotal += 1.0;
      break;
    case eStateEx:
      proc = 1;
      break;
    case eStateVoid:
      proc = 1;
      // printf("{ %d, %d } ; ", (int)numCellsToProcess, (int)cells[i]->idx);
      break;
    case eStateBound:
      // want to process boundary cells only if they adjoin a non-boundary, non-poly
      for(u8 ni = 0; ni<NUM_NEIGHBORS; ni++) {
	tmpState = cells[cells[i]->neighborIdx[ni]]->state;
	proc |= ((tmpState == eStateDrug) || (tmpState == eStateEx) || (tmpState == eStateVoid));
      }
      ////// DEBUG:
      // if(proc) { cellsAdjoining++; }
      //  dum++;
      /////////
      break;
    default:
      proc = 0;
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
      cells[i]->dissSteps = (u32)((f64)dissNSteps[np] * dissratescale);
      cells[i]->dissInc = 1.0 / (f64)(cells[i]->dissSteps);
			
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // scale this cell's dissolution probabilty by NPN.. 
      // something like this?
      const f64 npscale = ((NUM_NEIGHBORS - np) / NUM_NEIGHBORS);
      cells[i]->dissProb = (dissprob * npscale * disspolyscale) + (dissprob * (1.0 - disspolyscale));
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
  }
	
  // initialize the update data
  for(u32 i=0; i<numCells; i++) {
    *(cellsUpdate[i]) = *(cells[i]);
  }
  drugMass = drugMassTotal;
  
  ////////// DEBUG
  /*  
  // another loop over all cells to verify final cell distribution count
  for(u8 s=0; s<8; s++) { stateCount[s] = 0; }
  for (u32 n=0; n<numCells; n++) {
    stateCount[cells[n]->state]++;
  }
  // easier-to-read totals
  cellsInTablet   = stateCount[eStateDrug] + stateCount[eStateEx] + stateCount[eStatePoly] + stateCount[eStateVoid];
  drugCountRatio  = (f64)stateCount[eStateDrug] / (f64)cellsInTablet;
  polyCountRatio  = (f64)stateCount[eStatePoly] / (f64)cellsInTablet;
  exCountRatio    = (f64)stateCount[eStateEx] / (f64)cellsInTablet;
  voidCountRatio  = (f64)stateCount[eStateVoid] / (f64)cellsInTablet;
  // debugger hook
  dum++;
  */
  /////////////////
  
}



//// cell type distribution (on 8-cell blocks)
void CellModel::distribute(void) {
  //  unsigned long int* polyIdx = new unsigned long int[nPoly]; // array of polymer cell idx
  //  unsigned long int* drugIdx = new unsigned long int[nDrug]; // array of drug cell idx
  
  u32 i, j, k;      // temp cartesian coordinates
  u32 cX, cY, cZ;   // center of cylinder
  u32 cubeR2;       // squared radius of cylinder
  u32 shellR2;      // squared radius at the inner edge of the shell
  u32 r2;           // temp squared radius
  u32 dcX, dcY;     // temp distance from center, cartesian
  u32 idx;          // temp idx
  u32 boundH;             // cylinder height in cells
  u32 loBoundH, hiBoundH; // cylinder vertical bounds
  u32 loShellH, hiShellH; // shell vertical bounds
  
  u32 nPolyBlocks, nDrugBlocks;
  
  vector<u32> shellIdx;   // cell idx's adjoining boundary
  vector<u32> tabletIdx;  // cell idx's inside cylinder
  vector<u32> polyIdx;    
  vector<u32> drugIdx;    
  vector<u32> exIdx;   
  
  // cylinder dimensions
  cX = cubeLength >> 1;
  cY = cX;
  cZ = cX;
  cubeR2 = (cX-1) * (cX-1);
  shellR2 = (cX-1-(wShell*2)) * (cX-1-(wShell*2));
  
  // number of cells in cylinder height
  boundH = cylinderHeight * cubeLength;
  loBoundH = cZ - (boundH >> 1);
  hiBoundH = cZ + (boundH >> 1);
  loShellH = loBoundH + (wShell*2);
  hiShellH = hiBoundH - (wShell*2);
  
  // divide the space: interior, exterior, shell
  for(u32 i=0; i<cubeLength; i += 2) {
    for(u32 j=0; j<cubeLength; j += 2) {
      for(u32 k=0; k<cubeLength; k += 2) {
        idx = subToIdx(i, j, k);
        dcX = i-cX+1;
        dcY = j-cY+1;
        r2 = (dcX * dcX) + (dcY * dcY);
        if ((i < 1) 
            || (j < 1) 
            || (k < 1)
            || (i > (cubeLength-3))
            || (j > (cubeLength-3))
            || (k < loBoundH)
            || (k > hiBoundH)
            || (r2 >= cubeR2) ) {
 
	  //        if(r2 >= cubeR2) {
          // exterior cells
          this->setBlockState(idx, eStateBound);
        } // end outside-cylinder
        else { 
          if((r2 > shellR2) || (k < loShellH) || (k > hiShellH)) {
            shellIdx.push_back(idx);
          } else {
            tabletIdx.push_back(idx);
          }
        } // end outside-cylinder test
      } // k
    } // j
  } // i
  
  // shuffle the shell and tablet idx's 
  random_shuffle(shellIdx.begin(), shellIdx.end());
  random_shuffle(tabletIdx.begin(), tabletIdx.end());
  
  ///// distribute polymer cells
  u32 nBlocks = shellIdx.size() + tabletIdx.size();
  nPolyBlocks = (u32)((f64)nBlocks * pPoly);
  nDrugBlocks = (u32)((f64)nBlocks * pDrug);
  
  if(nPolyBlocks > (nBlocks - nDrugBlocks)) {
    nPolyBlocks = nBlocks - nDrugBlocks;
  }
  
  // use shellBalance as deviation from expected distribution in shell
  u32 nPolyInShell = (u32)( (f64)shellIdx.size() * pPoly * pShellBalance);
  // use shellBalance as proportion of total polymer
  //  u32 nPolyInShell = (u32)((float)nPolyBlocks * pShellBalance);
  // use shellBalance as proportion of total shell size
  // u32 nPolyInShell = shellIdx.size() * pShellBalance;
  
  nPolyInShell = min((unsigned int)nPolyInShell, (unsigned int)(shellIdx.size()));
  nPolyInShell = min((unsigned int)nPolyInShell, (unsigned int)nPolyBlocks);
  if (nPolyBlocks > tabletIdx.size()) {
    nPolyInShell = max(nPolyInShell, nPolyBlocks - tabletIdx.size());
  }
  u32 nPolyInTablet = nPolyBlocks - nPolyInShell;
  nPolyInTablet = min((unsigned int)nPolyInTablet, (unsigned int)(tabletIdx.size()));
  
  
  ///// DEBUG
  /*
  if ((nPolyInShell + nPolyInTablet) != nPolyBlocks) {
    int dum = 0;
    dum++;
  }
  */
  /////////////
  
  u32 n;
  for(n=0; n < nPolyInShell; n++) {
    polyIdx.push_back(shellIdx.back());
    shellIdx.pop_back();
  }

  for(n=0; n < nPolyInTablet; n++) {
    polyIdx.push_back(tabletIdx.back());
    tabletIdx.pop_back();
  }
  
  // reassign all remaining shell idx's to tablet
  const u32 shells = shellIdx.size();
  for(n=0; n<shells; n++) {
    tabletIdx.push_back(shellIdx.back());
    shellIdx.pop_back();
  }
	
  // re-shuffle
  random_shuffle(tabletIdx.begin(), tabletIdx.end());
	
  //// distribute drug cells (shell and tablet);
  for(n=0; n<nDrugBlocks; n++) {
    drugIdx.push_back(tabletIdx.back());
    tabletIdx.pop_back();
  }
  //// everything else becomes excipient
  const u32 tablets = tabletIdx.size();
  for( n=0; n<tablets; n++) {
    exIdx.push_back(tabletIdx.back());
    tabletIdx.pop_back();
  }
  
  // iterate over index vectors and fill 2x2x2 blocks
  vector<u32>::iterator it;
  for(it=polyIdx.begin(); it != polyIdx.end(); it++) {
    this->setBlockState(*it, eStatePoly);
  }
  for(it=drugIdx.begin(); it != drugIdx.end(); it++) {
    this->setBlockState(*it, eStateDrug);
  }
  for(it=exIdx.begin(); it != exIdx.end(); it++) {
    this->setBlockState(*it, eStateEx);
  }
}

/// compression step
void CellModel::compress(void) {
  u32 idx;
  u32 nIdx;
  u32 nIdx2;
  u8 diag;
  eCellState swapstate;
  
  // offset coordinates to get diagonals in 2x2x2
  const u8 diags[4][3]	= { {0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0} };
  // complementary diagonals
  const u8 diagsNot[4][3]	= { {1, 1, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  
  // each polymer block randomly swaps diagonal complements with neighbors
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
}

// set state of a 2x2x2 block of cells
void CellModel::setBlockState(const u32 idx, eCellState state) {
  u32 i, j, k;
  u32 nIdx;
  u8 l, m, n, diag;
  idxToSub(idx, &i, &j, &k);
  
  static const u8 diags[4][3]	= { {0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0} };
  // complementary diagonals
  static const u8 diagsNot[4][3]	= { {1, 1, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
	
  // force border cells to assume boundary state
  if ( i==0 || j==0 || k==0
       || i>(cubeLength-2) || j>(cubeLength-2) || k>(cubeLength-2) ) {
    state = eStateBound;
  }
  
  switch(state) {
  case eStatePoly:
  case eStateBound:
    for(l=0; l<2; l++) {
      for(m=0; m<2; m++) {
	for(n=0; n<2; n++) {
	  nIdx = subToIdx(i+l, j+m, k+n);
	  cells[nIdx]->state = state;
	}
      }
    }
    break;
  case eStateDrug:
    // drug cells: fill diagonals
    for(diag = 0; diag<4; diag++) {
      nIdx = subToIdx(i+diags[diag][0], j+diags[diag][1], k+diags[diag][2]);
      cells[nIdx]->state = eStateDrug;
    }
    // opposite diagonals are empty
    for(diag = 0; diag<4; diag++) {
      nIdx = subToIdx(i+diagsNot[diag][0], j+diagsNot[diag][1], k+diagsNot[diag][2]);
      cells[nIdx]->state = eStateVoid;
    }
    break;
  case eStateEx:
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
    break;
  }
}
  
