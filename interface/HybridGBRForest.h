
#ifndef EGAMMAOBJECTS_HybridGBRForest
#define EGAMMAOBJECTS_HybridGBRForest

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// HybridGBRForest                                                            //
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

#include <vector>
#include "HybridGBRTree.h"
#include <math.h>
#include <stdio.h>
#include "Rtypes.h"

  class HybridGBRForest {

    public:

       HybridGBRForest() {}      
       HybridGBRForest(int ntargets);
       virtual ~HybridGBRForest();
       
       int NTargets() const { return fTrees.size(); }
       
       //void GetResponse(const float* vector) const;
       //double GetResponse(int idx) const { return fResponses[idx]; }
       
       double GetResponse(const float* vector, int idx) const;
       
       void SetInitialResponse(int idx, double response) { fInitialResponse[idx] = response; }
       
       std::vector<std::vector<HybridGBRTree> > &Trees() { return fTrees; }
       const std::vector<std::vector<HybridGBRTree> > &Trees() const { return fTrees; }
       
    protected:
      std::vector<double> fInitialResponse;
      //mutable std::vector<double> fResponses;
      std::vector<std::vector<HybridGBRTree> > fTrees;  
      
    private:

      ClassDef(HybridGBRForest,2)             
      
  };

//_______________________________________________________________________
// inline void HybridGBRForest::GetResponse(const float* vector) const {
//   for (unsigned int itgt=0; itgt<fResponses.size(); ++itgt) {
//     fResponses[itgt] = fInitialResponse[itgt];
//     for (std::vector<HybridGBRTree>::const_iterator it=fTrees[itgt].begin(); it!=fTrees[itgt].end(); ++it) {
//       int termidx = it->TerminalIndex(vector);
//       fResponses[itgt] += it->GetResponse(termidx);
//     }    
//   }
// }

//_______________________________________________________________________
inline double HybridGBRForest::GetResponse(const float* vector, int idx) const {
  double response = fInitialResponse[idx];
  for (std::vector<HybridGBRTree>::const_iterator it=fTrees[idx].begin(); it!=fTrees[idx].end(); ++it) {
    int termidx = it->TerminalIndex(vector);
    response += it->GetResponse(termidx);
  }    
  return response;
}


#endif
