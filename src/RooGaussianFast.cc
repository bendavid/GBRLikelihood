#include <iostream>
#include <math.h>
#include "TMath.h"

//#include "../interface/RooGaussianFast.h"
//#include "../interface/RooFermi.h"
//#include "../interface/RooRelBW.h"
#include "../interface/RooGaussianFast.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "../interface/GBRMath.h"

using namespace RooFit;

 ClassImp(RooGaussianFast) 

 RooGaussianFast::RooGaussianFast(){}

 RooGaussianFast::RooGaussianFast(const char *name, const char *title, 
                    RooAbsReal& _x,
                    RooAbsReal& _mean,
                    RooAbsReal& _width
                    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width)
 { 
 } 


 RooGaussianFast::RooGaussianFast(const RooGaussianFast& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width)

 { 
 } 
 
 double RooGaussianFast::evaluate() const 
 { 
   double t = (x-mean)*vdt::fast_inv(width);
   return vdt::fast_exp(-0.5*t*t);    
 } 

 Int_t RooGaussianFast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const 
 {
   if (matchArgs(allVars,analVars,x)) return 1;
   return 0;
 }

 Double_t RooGaussianFast::analyticalIntegral(Int_t code, const char* rangeName) const 
 {
   assert(code==1) ;
   
   static const double rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
   static const double root2Pi = 2.0*rootPiBy2;
   static const double invRoot2 = 1.0/sqrt(2);
 
   double xmin = x.min(rangeName);
   double xmax = x.max(rangeName);   
   double invwidth = vdt::fast_inv(width);
   
   double tmin = (xmin-mean)*invwidth;
   double tmax = (xmax-mean)*invwidth;
   
   bool isfullrange = (tmin<-1000. && tmax>1000.);

   if (isfullrange) {
     return width*root2Pi;
   }
   else {
     return rootPiBy2*width*(TMath::Erf(tmax*invRoot2)-TMath::Erf(tmin*invRoot2));
   }
 
 }
 