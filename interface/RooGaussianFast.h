#ifndef ROOGAUSSIANFAST
#define ROOGAUSSIANFAST

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooGaussianFast : public RooAbsPdf {
public:
  RooGaussianFast();
  RooGaussianFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width
           );
  RooGaussianFast(const RooGaussianFast& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooGaussianFast(*this,newname); }
  inline virtual ~RooGaussianFast() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooGaussianFast,1)
};
#endif
