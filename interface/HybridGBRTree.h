
#ifndef EGAMMAOBJECTS_HybridGBRTree
#define EGAMMAOBJECTS_HybridGBRTree

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBRForest                                                            //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// Designed to be built from TMVA-trained trees, but could also be      //
// generalized to otherwise-trained trees, classification,              //
//  or other boosting methods in the future                             //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

// The decision tree is implemented here as a set of two arrays, one for
// intermediate nodes, containing the variable index and cut value, as well
// as the indices of the 'left' and 'right' daughter nodes.  Positive indices
// indicate further intermediate nodes, whereas negative indices indicate
// terminal nodes, which are stored simply as a vector of regression responses


#include <vector>
#include <map>
#include <stdio.h>
#include <cmath>
#include "Rtypes.h"

  class HybridGBRTree {

    public:

       HybridGBRTree() {}
       virtual ~HybridGBRTree();
              
       //double GetResponse(const float* vector) const;
       double GetResponse(int termidx) const { return fResponses[termidx]; }
       int TerminalIndex(const float *vector) const;
              
       std::vector<float> &Responses() { return fResponses; }       
       const std::vector<float> &Responses() const { return fResponses; }
       
       std::vector<unsigned char> &CutIndices() { return fCutIndices; }
       const std::vector<unsigned char> &CutIndices() const { return fCutIndices; }
       
       std::vector<float> &CutVals() { return fCutVals; }
       const std::vector<float> &CutVals() const { return fCutVals; }
       
       std::vector<int> &LeftIndices() { return fLeftIndices; }
       const std::vector<int> &LeftIndices() const { return fLeftIndices; } 
       
       std::vector<int> &RightIndices() { return fRightIndices; }
       const std::vector<int> &RightIndices() const { return fRightIndices; }
       
       std::vector<std::vector<std::pair<float,float> > > &Limits() { return fLimits; }
       const std::vector<std::vector<std::pair<float,float> > > &Limits() const { return fLimits; }
       
    protected:              
	std::vector<unsigned char> fCutIndices;
	std::vector<float> fCutVals;
	std::vector<int> fLeftIndices;
	std::vector<int> fRightIndices;
	std::vector<float> fResponses;  
        std::vector<std::vector<std::pair<float,float> > > fLimits;

    private:

      ClassDef(HybridGBRTree,2)        
        
  };


//_______________________________________________________________________
inline int HybridGBRTree::TerminalIndex(const float* vector) const {
  
  int index = 0;
  
  unsigned char cutindex = fCutIndices[0];
  float cutval = fCutVals[0];
  
  while (true) {
    //printf("cutindex = %i, cutval = %5f, val = %5f\n",cutindex,cutval,vector[cutindex]);
    if (vector[cutindex] > cutval) {
      //printf("right\n");
      index = fRightIndices[index];
    }
    else {
      //printf("left\n");
      index = fLeftIndices[index];
    }
    
    if (index>0) {
      cutindex = fCutIndices[index];
      cutval = fCutVals[index];
    }
    else {
      //printf("terminal\n");
      return (-index);
    }
    
  }
  

}

#endif
