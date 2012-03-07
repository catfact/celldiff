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

void CellModel::setup(void) {
  this->distribute();
  this->compress();
  
}

void CellModel::distribute(void) {
//  unsigned long int* polyIdx = new unsigned long int[nPoly]; // array of polymer cell idx
//  unsigned long int* drugIdx = new unsigned long int[nDrug]; // array of drug cell idx

  u32 i, j, k;  // temp cartesian coordinates
  u32 cX, cY;   // center of cylinder
  u32 cubeR2;   // squared radius of cylinder
  u32 edgeR2;   // squared radius at the inner edge of the shell
  u32 r2;       // temp squared radius
  u32 dcX, dcY; // temp distance from center, cartesian
  u32 idx;      // temp idx
  u32 numH;     // cylinder height (in cells)
  u32 edgeH;     // edge height
  
  vector<u32> shellIdx;   // cell idx's adjoining boundary
  vector<u32> tabletIdx;  // cell idx's inside cylinder
  vector<u32> polyIdx;    
  vector<u32> drugIdx;    
  vector<u32> exIdx;   
  
  // cylinder dimensions
  cX = cubeLength >> 1;
  cY = cX;
  cubeR2 = (cX-1) * (cX-1);
  edgeR2 = (cX-1-wShell) * (cX-1-wShell);
  
  // number of cells in cylinder height
  numH = cylinderHeight * cubeLength;
  edgeH = numH - wShell;
  
  // divide the space: interior, exterior, edge
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
            || (k > (cubeLength-3))
          || (r2 >= cubeR2) ) {
          // exterior cells
          this->setBlockState(idx, eStateBound);
        } // end outside-cylinder
        else { 
          if((r2 > edgeR2) || (k < wShell) || (k > edgeH)) {
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
  // in shell
  nPoly = 
  u32 nPolyInShell = (u32)( (f64)(shellIdx.size()) * (f64)nPoly / (f64)(shellIdx.size() + tabletIdx.size()) * pShellB );
  u32 n;
  for(n=0; n < nPolyInShell; n++) {
    polyIdx.push_back(shellIdx.back());
    shellIdx.pop_back();
  }
  // in tablet
  for(n=nPolyInShell; n < nPoly; n++) {
    polyIdx.push_back(tabletIdx.back());
    tabletIdx.pop_back();
  }
  // reassign all remaining shell idx's to tablet
  for(n=0; n<shellIdx.size(); n++) {
    tabletIdx.push_back(shellIdx.back());
    shellIdx.pop_back();
  }
  //// distribute drug cells (shell and tablet);
  for(n=0; n<nDrug; n++) {
    drugIdx.push_back(tabletIdx.back());
    tabletIdx.pop_back();
  }
  //// the rest of the in-tablet idx's become excipient
  for( n=0; n<tabletIdx.size(); n++) {
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
  static const u8 diags[4][3]	= { {0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0} };
  // complementary diagonals
  static const u8 diagsNot[4][3]	= { {1, 1, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  
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
void CellModel::setBlockState(const u32 idx, const eCellState state) {
  u32 i, j, k;
  u32 nIdx;
  u8 l, m, n;
  idxToSub(idx, &i, &j, &k);
  for(l=0; l<2; l++) {
    for(m=0; m<2; m++) {
      for(n=0; n<2; n++) {
        nIdx = subToIdx(i+l, j+m, k+n);
        cells[nIdx]->state = state;
      }
    }
  }
}

