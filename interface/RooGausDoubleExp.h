#ifndef ROOGAUSDOUBLEEXP
#define ROOGAUSDOUBLEEXP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooGausDoubleExp : public RooAbsPdf {
public:
  RooGausDoubleExp();
  RooGausDoubleExp(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _alpha2
           );
  RooGausDoubleExp(const RooGausDoubleExp& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooGausDoubleExp(*this,newname); }
  inline virtual ~RooGausDoubleExp() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy alpha2;

  Double_t evaluate() const ;

private:

  ClassDef(RooGausDoubleExp,1)
};
#endif
