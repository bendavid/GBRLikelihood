#include <iostream>
#include <math.h>
#include "TMath.h"

//#include "../interface/RooFermi.h"
//#include "../interface/RooRelBW.h"
#include "../interface/RooGausDoubleExp.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "../interface/GBRMath.h"

using namespace RooFit;

 ClassImp(RooGausDoubleExp)

 RooGausDoubleExp::RooGausDoubleExp(){}

 RooGausDoubleExp::RooGausDoubleExp(const char *name, const char *title,
                    RooAbsReal& _x,
                    RooAbsReal& _mean,
                    RooAbsReal& _width,
                    RooAbsReal& _alpha1,
                    RooAbsReal& _alpha2
                    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   alpha1("alpha1","alpha1",this,_alpha1),
   alpha2("alpha2","alpha2",this,_alpha2)
 {
 } 


 RooGausDoubleExp::RooGausDoubleExp(const RooGausDoubleExp& other, const char* name) :
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   alpha1("alpha1",this,other.alpha1),
   alpha2("alpha2",this,other.alpha2)
 {
 }
 
 double RooGausDoubleExp::evaluate() const
 { 
   double t = (x-mean)*vdt::fast_inv(width);
   if(t>-alpha1 && t<alpha2){
     return vdt::fast_exp(-0.5*t*t);
   }else if(t<=-alpha1){
     double A1 = vdt::fast_exp(0.5*alpha1*alpha1);
     return A1*vdt::fast_exp(alpha1*t);
   }else if(t>=alpha2){
     double A2 = vdt::fast_exp(0.5*alpha2*alpha2);
     return A2*vdt::fast_exp(-alpha2*t);
   }//else{
     //cout << "ERROR evaluating range..." << endl;
   return -99.;
   //}
    
 } 

 Int_t RooGausDoubleExp::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const
 {
   if (matchArgs(allVars,analVars,x)) return 1;
   return 0;
 }

 Double_t RooGausDoubleExp::analyticalIntegral(Int_t code, const char* rangeName) const
 {
   assert(code==1) ;
 
   double central=0;
   double left=0;
   double right=0;
 
   static const double root2 = sqrt(2) ;
   static const double rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
   double xscale = root2*width;
 
   double invwidth = vdt::fast_inv(width);
   
   //compute gaussian contribution
   double central_low =std::max(x.min(rangeName),mean - alpha1*width );
   double central_high=std::min(x.max(rangeName),mean + alpha2*width );
   if(central_low < central_high) // is the gaussian part in range?
     central = rootPiBy2*width*(TMath::Erf((central_high-mean)/xscale)-TMath::Erf((central_low-mean)/xscale));
 
   //compute left tail;
   double A1 = vdt::fast_exp(0.5*alpha1*alpha1);

   double left_low=x.min(rangeName);
   double left_high=std::min(x.max(rangeName),mean - alpha1*width);
   if(left_low < left_high){ //is the left tail in range?
     left = A1*width*vdt::fast_inv(alpha1)*(vdt::fast_exp((left_high-mean)*invwidth*alpha1) - vdt::fast_exp((left_low-mean)*invwidth*alpha1));
   }
 
   //compute right tail;
   double A2 = vdt::fast_exp(0.5*alpha2*alpha2);

   double right_low=std::max(x.min(rangeName),mean + alpha2*width);
   double right_high=x.max(rangeName);
   if(right_low < right_high){ //is the right tail in range?
     right = -A2*width*vdt::fast_inv(alpha2)*(vdt::fast_exp(-(right_high-mean)*invwidth*alpha2) - vdt::fast_exp(-(right_low-mean)*invwidth*alpha2));
   }

   //printf("left = %5f, central = %5f, right = %5f\n",left,central,right);
     
   return left+central+right;
 
 }
 