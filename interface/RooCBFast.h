#ifndef ROOCBFAST
#define ROOCBFAST

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooCBFast : public RooAbsPdf {
public:
  RooCBFast();
  RooCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1
           );
  RooCBFast(const RooCBFast& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCBFast(*this,newname); }
  inline virtual ~RooCBFast() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooCBFast,1)
};
#endif
