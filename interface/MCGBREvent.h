#ifndef EGAMMAOBJECTS_MCGBREvent
#define EGAMMAOBJECTS_MCGBREvent

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MCGBREvent                                                             //
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
  
  class MCGBREvent {

    public:

       MCGBREvent(int nvars); 
       MCGBREvent(int nvars, const MCGBREvent &other) :
        fFuncVal(other.fFuncVal),
        fFuncValAlt(other.fFuncValAlt),
        fWeight(other.fWeight),
        fTarget(other.fTarget),
        fTargetMin(other.fTargetMin),
        fTarget3(other.fTarget3),
        fArg(other.fArg),
        fArgLog(other.fArgLog) {
	  
          fVars = new float[nvars];
          fQuantiles = new int[nvars];
        	  
          for (int ivar=0; ivar<nvars; ++ivar) {
            fVars[ivar] = other.fVars[ivar];
            fQuantiles[ivar] = other.fQuantiles[ivar];
          }
                    
        }
        
       ~MCGBREvent();
       
       float Var(int i) const { return fVars[i]; }
       int Quantile(int i) const { return fQuantiles[i]; }   
       double Weight() const   { return fWeight;  }
       
       void SetVar(int i, float x) { fVars[i] = x; }
       void SetQuantile(int i, int q) { fQuantiles[i] = q; }
       //cache computed qunatities needed for split-search
       void SetWeight(double x)     { fWeight = x; }
       
       double FuncVal() const { return fFuncVal; }
       void SetFuncVal(double x) { fFuncVal = x; }

       double FuncValAlt() const { return fFuncValAlt; }
       void SetFuncValAlt(double x) { fFuncValAlt = x; }       
       
       double Target() const { return fTarget; }
       void SetTarget(double x) { fTarget = x; }       
       
       double TargetMin() const { return fTargetMin; }
       void SetTargetMin(double x) { fTargetMin = x; }         
       
       double Target3() const { return fTarget3; }
       void SetTarget3(double x) { fTarget3 = x; }
       
       double Arg() const { return fArg; }
       void SetArg(double x) { fArg = x; }
       
       double ArgLog() const { return fArgLog; }
       void SetArgLog(double x) { fArgLog = x; }
       
       double ArgSig1() const { return fArgSig1; }
       void SetArgSig1(double x) { fArgSig1 = x; }

       double ArgSig2() const { return fArgSig2; }
       void SetArgSig2(double x) { fArgSig2 = x; }
       
       double ArgSig3() const { return fArgSig3; }
       void SetArgSig3(double x) { fArgSig3 = x; }
       
       double ArgSig4() const { return fArgSig4; }
       void SetArgSig4(double x) { fArgSig4 = x; }
       
    protected:
      float                    *fVars;
      int                      *fQuantiles;
      double                    fFuncVal;
      double                    fFuncValAlt;
      double                    fWeight;
      double                    fTarget;
      double                    fTargetMin;
      double                    fTarget3;
      double                    fArg;
      double                    fArgLog;
      double                    fArgSig1;
      double                    fArgSig2;
      double                    fArgSig3;
      double                    fArgSig4;
  };
  
  
/*  class GBRTargetCMP : public std::binary_function<MCGBREvent*, MCGBREvent*, bool> {
    public:
      GBRTargetCMP() {}
      bool operator() (const MCGBREvent *ev1, const MCGBREvent *ev2) const { return ev1->Target()<ev2->Target() ? true : false; }
  };
  
  class GBRAbsTargetCMP : public std::binary_function<MCGBREvent*, MCGBREvent*, bool> {
    public:
      GBRAbsTargetCMP() {}
      bool operator() (const MCGBREvent *ev1, const MCGBREvent *ev2) const { return fabs(ev1->Target())<fabs(ev2->Target()) ? true : false; }
  }; */ 
  
  class MCGBRVarCMP : public std::binary_function<MCGBREvent*, MCGBREvent*, bool> {
    public:
      MCGBRVarCMP() {}
      MCGBRVarCMP(int idx) : fVarIdx(idx) {}      
      bool operator() (const MCGBREvent *ev1, const MCGBREvent *ev2) const { return ev1->Var(fVarIdx)<ev2->Var(fVarIdx) ? true : false; }
      
    protected:
      int fVarIdx;
  };
  
  class MCGBRTgtCMP : public std::binary_function<MCGBREvent*, MCGBREvent*, bool> {
    public:
      MCGBRTgtCMP() {}
      bool operator() (const MCGBREvent *ev1, const MCGBREvent *ev2) const { return ev1->FuncVal()-ev1->Target() < ev2->FuncVal()-ev2->Target() ? true : false; }
//       bool operator() (const MCGBREvent *ev1, const MCGBREvent *ev2) const { return ev1->FuncVal()-ev1->TargetMin() < ev2->FuncVal()-ev2->TargetMin() ? true : false; }

      
  };  
  
  
#endif
