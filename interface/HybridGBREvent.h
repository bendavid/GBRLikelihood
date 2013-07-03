
#ifndef EGAMMAOBJECTS_HybridGBREvent
#define EGAMMAOBJECTS_HybridGBREvent

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// HybridGBREvent                                                             //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// This is a helper class for GBRTrainer to store  needed information   //
// in memory and facilitate sorting according to target or input        //
// variables.                                                           //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <stdio.h>
#include <algorithm>
#include "TVectorD.h"  
#include "TMatrixDSym.h"
  
  class HybridGBREvent {

    public:

       HybridGBREvent(int nvars, int ntargets = 1, int nparms=0); 
       HybridGBREvent(int nvars, int ntargets, const HybridGBREvent &other, int nparms=0) :
        fValVector(other.fValVector),
        fParmMatrix(other.fParmMatrix),
        fPdfVal(other.fPdfVal),
        fTarget(other.fTarget),
        fTransTarget(other.fTransTarget),
        fWeight(other.fWeight),
        fWeightedTransTarget(other.fWeightedTransTarget),
        fWeightedTransTarget2(other.fWeightedTransTarget2),
        fInverseGenPdf(other.fInverseGenPdf),
        fClass(other.fClass),
        fCurrentNode(other.fCurrentNode),
        fIsPrimary(other.fIsPrimary) {
         
          fVars = new float[nvars];
          fVarsAlt = new float[nvars];
          fQuantiles = new int[nvars];
          fTargets = new float[ntargets];
          fTransTargets = new float[ntargets];
          fSmoothedTargets = new float[ntargets];
          fDerivatives = new double[nparms];
	  fDerivatives2 = new double[nparms];
	  fCurrentNodes =  new unsigned int[ntargets];
	  
          for (int ivar=0; ivar<nvars; ++ivar) {
            fVars[ivar] = other.fVars[ivar];
            fVarsAlt[ivar] = other.fVarsAlt[ivar];
            fQuantiles[ivar] = other.fQuantiles[ivar];
          }
          
          for (int itgt=0; itgt<ntargets; ++itgt) {
            fTargets[itgt] = other.fTargets[itgt];
            fTransTargets[itgt] = other.fTransTargets[itgt];
	    fTransTargets2[itgt] = other.fTransTargets2[itgt];	    
	    fSmoothedTargets[itgt] = other.fSmoothedTargets[itgt];
	    fCurrentNodes[itgt] = other.fCurrentNodes[itgt];
          }
          
          for (int iparm=0; iparm<nparms; ++iparm) {
	    fDerivatives[iparm] = other.fDerivatives[iparm];
	    fDerivatives2[iparm] = other.fDerivatives2[iparm];
	  }
          
        }
        
       ~HybridGBREvent();
       
       bool IsPrimary()    const { return true; }
       float InverseGenPdf() const { return fInverseGenPdf; }
       unsigned int CurrentNode() const { return fCurrentNode; }
       unsigned char Class() const { return fClass; }
       float Var(int i) const { return fVars[i]; }
       float VarAlt(int i) const { return fVarsAlt[i]; }       
       unsigned short Quantile(int i) const { return fQuantiles[i]; }
       float Target() const   { return Target(2);  }
       float TransTarget() const { return TransTarget(2);  }       
       float Weight() const   { return fWeight;  }
       float WeightedTransTarget() const { return fWeight*TransTarget(2); }
       float WeightedTransTarget2() const { return fWeight*TransTarget(2)*TransTarget(2); }
       
       void SetClass(unsigned char i) { fClass = i; }
       void SetVar(int i, float x) { fVars[i] = x; }
       void SetVarAlt(int i, float x) { fVarsAlt[i] = x; }       
       void SetQuantile(int i, int q) { fQuantiles[i] = q; }
       void SetTarget(float x)     { fTarget = x; }
       //cache computed qunatities needed for split-search
       void SetTransTarget(float x) { fTransTarget = x; fWeightedTransTarget = fWeight*fTransTarget; fWeightedTransTarget2 = fWeightedTransTarget*fTransTarget; }
       void SetWeight(float x)     { fWeight = x; }
       void SetCurrentNode(unsigned int i) { fCurrentNode = i; }
       void SetInverseGenPdf(float x) { fInverseGenPdf = x; }
       void SetIsPrimary(bool b) { fIsPrimary = b; }
       
       float Target(int i) const { return fTargets[i]; }
       float TransTarget(int i) const { return fTransTargets[i]; }
       float SmoothedTarget(int i) const { return fSmoothedTargets[i];  }       
       void SetTarget(int i, float x) { fTargets[i] = x; }
       void SetTransTarget(int i, float x) { fTransTargets[i] = x; }
       void SetSmoothedTarget(int i, float x) { fSmoothedTargets[i] = x; }
       
       double PdfVal() const { return fPdfVal; }
       void SetPdfVal(double x) { fPdfVal = x; }

       double Derivative(int i) const { return fDerivatives[i];  }       
       void SetDerivative(int i, double x) { fDerivatives[i] = x; }     
       
       double Derivative2(int i) const { return fDerivatives2[i];  }       
       void SetDerivative2(int i, double x) { fDerivatives2[i] = x; }        

       double TransTarget2(int i) const { return fTransTargets2[i];  }       
       void SetTransTarget2(int i, double x) { fTransTargets2[i] = x; }        
       
       unsigned int CurrentNode(int i) const { return fCurrentNodes[i];  }       
       void SetCurrentNode(int i, unsigned int x) { fCurrentNodes[i] = x; }         
       
       std::vector<std::pair<const HybridGBREvent*, float> > &NearestNeighbours() { return fNearestNeighbours; }
       
       double SumNNWeights() const { return fSumNNWeights; }
       void SetSumNNWeights(double x) { fSumNNWeights = x; }
       
       TVectorD &ValVector() { return fValVector; }
       TMatrixDSym &ParmMatrix() { return fParmMatrix; }
       
       
    protected:
      float                    *fVars;
      float                    *fVarsAlt;      
      float                    *fTargets;
      float                    *fTransTargets;
      float                    *fTransTargets2;      
      float                    *fSmoothedTargets;      
      int                      *fQuantiles;
      double                   *fDerivatives;
      double                   *fDerivatives2;
      unsigned int             *fCurrentNodes;
      TVectorD                  fValVector;
      TMatrixDSym               fParmMatrix;
      double                    fPdfVal;
      int                       fOldestIdx;
      float                     fTarget;
      float                     fTransTarget;
      float                     fWeight;
      float                     fWeightedTransTarget;
      float                     fWeightedTransTarget2;
      float                     fInverseGenPdf;
      unsigned char             fClass;
      unsigned int		fCurrentNode;
      bool                      fIsPrimary;
      double                    fSumNNWeights;
      std::vector<std::pair<const HybridGBREvent*, float> > fNearestNeighbours;
  };
  
  
  class GBRTargetCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRTargetCMP() {}
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return ev1->Target()<ev2->Target() ? true : false; }
  };
  
  class GBRAbsTargetCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRAbsTargetCMP() {}
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return fabs(ev1->Target())<fabs(ev2->Target()) ? true : false; }
  };  
  
  class GBRVarCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRVarCMP() {}
      GBRVarCMP(int idx) : fVarIdx(idx) {}      
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return ev1->Var(fVarIdx)<ev2->Var(fVarIdx) ? true : false; }
      
    protected:
      int fVarIdx;
  };
  
  class GBRVarAltCMP : public std::binary_function<HybridGBREvent*, HybridGBREvent*, bool> {
    public:
      GBRVarAltCMP() {}
      GBRVarAltCMP(int idx) : fVarIdx(idx) {}      
      bool operator() (const HybridGBREvent *ev1, const HybridGBREvent *ev2) const { return ev1->VarAlt(fVarIdx)<ev2->VarAlt(fVarIdx) ? true : false; }
      
    protected:
      int fVarIdx;
  };  
  
#endif
