#ifndef ROOCBEXP
#define ROOCBEXP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooCBExp : public RooAbsPdf {
public:
  RooCBExp();
  RooCBExp(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2
           );
  RooCBExp(const RooCBExp& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCBExp(*this,newname); }
  inline virtual ~RooCBExp() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooCBExp,1)
};
#endif
