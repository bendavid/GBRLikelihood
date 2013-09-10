#ifndef ROOREVCBFAST
#define ROOREVCBFAST

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooRevCBFast : public RooAbsPdf {
public:
  RooRevCBFast();
  RooRevCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2
           );
  RooRevCBFast(const RooRevCBFast& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooRevCBFast(*this,newname); }
  inline virtual ~RooRevCBFast() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooRevCBFast,1)
};
#endif
