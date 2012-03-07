/*
 *  CellModelSetup.cpp
 *  celldiff
 *
 *  i'm putting CellModel initialization methods in their own file, just for readability
 *
 *  Created by Ezra Buchla on 03/05/2012
 */

#include <vector>
#include "CellModel.hpp"

void CellModel::setup(void) {
  this->distribute();
  this->compress();
}

void CellModel::distribute(void) {
  u32* polyIdx = new u32[nPoly]; // array of polymer cell idx
  u32* drugIdx = new u32[nDrug]; // array of drug cell idx
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
  
  std::vector<u32> shellIdx; // vector of cell idx's adjoining boundary
  std::vector<u32> tabletIdx; // vector of all cell idx's inside cylinder
  
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
          setBlockState(i, j, k, eStateBound);
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
  u32 nPolyInShell = (u32)( (f64)(shellIdx.size()) * (f64)nPoly / (f64)(shellIdx.size() + tabletIdx.size()) * pShellB );
  u32 pi;
  for(pi=0; pi < nPolyInShell; pi++) {
    polyIdx[pi] = shellIdx.back();
    shellIdx.pop_back();
  }
  // in tableta
  for(pi=nPolyInShell; pi < nPoly; pi++) {
    polyIdx[pi] = tabletIdx.back();
    tabletIdx.pop_back();
  }
  // reassign all remaining shell idx's to tablet
  for(u32 n=0; n<shellIdx.size(); n++) {
    tabletIdx.push_back(shellIdx.back());
    shellIdx.pop_back();
  }
  //// distribute drug cells (shell and tablet);
  //// all remaining shell+tablet idx's get excipient
  delete polyIdx;
  delete drugIdx;
  delete excipIdx;
  
}

void CellModel::compress(void) {
}

// utility to set state of a 2x2x2 block of cells
void CellModel::setBlockState(const u32 i, const u32 j, const u32 k, const eCellState state) {
  u32 nIdx;
  u8 l, m, n;
  for(l=0; l<2; l++) {
    for(m=0; m<2; m++) {
      for(n=0; n<2; n++) {
        nIdx = subToIdx(i+l, j+m, k+n);
        cells[nIdx]->state = eStateBound;
      }
    }
  }
}