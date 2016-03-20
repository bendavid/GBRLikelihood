
#ifndef GBRLIKELIHOOD_MCGBRForest
#define GBRLIKELIHOOD_MCGBRForest

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MCGBRForest                                                            //
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
#include "MCGBRTreeD.h"
#include <math.h>
#include <stdio.h>
#include "Rtypes.h"

  class MCGBRForest {
    public:
       typedef MCGBRTreeD TreeT;

       MCGBRForest();   
       virtual ~MCGBRForest();
       
       int NTargets() const { return fTrees.size(); }
       
       double GetResponse(const float* vector) const;
       double GetResponse(const float* vector, double &responsemin, double &response3) const;
       
       bool ResponseGEQ(const float* vector, double thres, double &response) const;
       
       double InitialResponse() const { return fInitialResponse; }
       void SetInitialResponse(double response) { fInitialResponse = response; }
       
       double InitialResponseMin() const { return fInitialResponseMin; }
       void SetInitialResponseMin(double response) { fInitialResponseMin = response; }       
       
       std::vector<MCGBRTreeD> &Trees() { return fTrees; }
       const std::vector<MCGBRTreeD> &Trees() const { return fTrees; }
       
    protected:
      double fInitialResponse;
      double fInitialResponseMin;
      std::vector<MCGBRTreeD> fTrees;  
      
    private:

      ClassDef(MCGBRForest,1)       
      
  };

//_______________________________________________________________________
inline double MCGBRForest::GetResponse(const float* vector) const {
  double response = fInitialResponse;
  for (std::vector<MCGBRTreeD>::const_iterator it=fTrees.begin(); it!=fTrees.end(); ++it) {
    int termidx = it->TerminalIndex(vector);
    response += it->GetResponse(termidx);
  }    
  return response;
}

//_______________________________________________________________________
inline double MCGBRForest::GetResponse(const float* vector, double &responsemin, double &response3) const {
  double response = fInitialResponse;
//   responsemin = std::max(0.,fInitialResponse);
  responsemin = fInitialResponseMin;
  response3 = 0.;
  for (std::vector<MCGBRTreeD>::const_iterator it=fTrees.begin(); it!=fTrees.end(); ++it) {
    int termidx = it->TerminalIndex(vector);
    double treeresponse = it->GetResponse(termidx);
//     response += it->GetResponse(termidx);
    response += treeresponse;
//     absresponse += std::max(0.,treeresponse);
//     responsemin += std::abs(treeresponse);
    responsemin += it->GetResponseMin(termidx);
    response3 += it->GetResponse3(termidx);
  }    
  return response;
}


//_______________________________________________________________________
inline bool MCGBRForest::ResponseGEQ(const float* vector, double thres, double &response) const {
  response = fInitialResponse;
//   printf("response = %5e, thres = %5e\n",response,thres);
  if (response<thres) return false;
  for (std::vector<MCGBRTreeD>::const_iterator it=fTrees.begin(); it!=fTrees.end(); ++it) {
    int termidx = it->TerminalIndex(vector);
    response += it->GetResponse(termidx);
    if (response<thres) return false;
  }
  return true;
}

#endif
