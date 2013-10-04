#include <iostream>
#include <math.h>
#include "TMath.h"

//#include "../interface/RooRevCBFast.h"
//#include "../interface/RooFermi.h"
//#include "../interface/RooRelBW.h"
#include "../interface/RooRevCBFast.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "../interface/GBRMath.h"

using namespace RooFit;

 ClassImp(RooRevCBFast) 

 RooRevCBFast::RooRevCBFast(){}

 RooRevCBFast::RooRevCBFast(const char *name, const char *title, 
                    RooAbsReal& _x,
                    RooAbsReal& _mean,
                    RooAbsReal& _width,
                    RooAbsReal& _alpha2,
                    RooAbsReal& _n2
                    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   alpha2("alpha2","alpha2",this,_alpha2),
   n2("n2","n2",this,_n2)
 { 
 } 


 RooRevCBFast::RooRevCBFast(const RooRevCBFast& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   alpha2("alpha2",this,other.alpha2),
   n2("n2",this,other.n2)

 { 
 } 
 
 double RooRevCBFast::evaluate() const 
 { 
   double t = (x-mean)*vdt::fast_inv(width);
   double val = -99.;
   if(t<alpha2){
     val = vdt::fast_exp(-0.5*t*t);
   }
   else {
     double alpha2invn2 = alpha2*vdt::fast_inv(n2);
     val = vdt::fast_exp(-0.5*alpha2*alpha2)*gbrmath::fast_pow(1. - alpha2invn2*(alpha2-t), -n2);     
     
//      double n2invalpha2 = n2*vdt::fast_inv(fabs(alpha2));
//      double A2 = gbrmath::fast_pow(n2invalpha2,n2)*vdt::fast_exp(-alpha2*alpha2/2.);
//      double B2 = n2invalpha2-fabs(alpha2);
//      val = A2*gbrmath::fast_pow(B2+t,-n2);
   }//else{
     //cout << "ERROR evaluating range..." << endl;
     
//    if (!std::isnormal(val)) {
//      printf("bad val: x = %5f, t = %5f, mean = %5f, sigma = %5f, alpha2 = %5f, n2 = %5f\n",double(x), t, double(mean),double(width),double(alpha2), double(n2));
//      printf("val = %5f\n",val);
//    }
     
   return val;
   //return std::max(double(std::numeric_limits<float>::min()),val);
   //return std::max(1e-3,val);
   //}
    
 } 

 Int_t RooRevCBFast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const 
 {
   if (matchArgs(allVars,analVars,x)) return 1;
   return 0;
 }

 Double_t RooRevCBFast::analyticalIntegral(Int_t code, const char* rangeName) const 
 {
   assert(code==1) ;
 
   double central=0;
   double left=0;
   double right=0;
   
   double xmin = x.min(rangeName);
   double xmax = x.max(rangeName);
   
   //printf("xmin = %5e, xmax = %5e\n",xmin,xmax);
   
   //bool isfullrange = xmin<=-RooNumber::infinity() && xmax>=RooNumber::infinity();
 
   static const double rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
   static const double invRoot2 = 1.0/sqrt(2);   
 
   double invwidth = vdt::fast_inv(width);
   
   double tmin = (xmin-mean)*invwidth;
   double tmax = (xmax-mean)*invwidth;
   
   bool isfullrange = (tmin<-1000. && tmax>1000.);
   
   //compute gaussian contribution
   double central_low = xmin;
   double central_high=std::min(xmax,mean + alpha2*width );
   if(central_low < central_high)  {// is the gaussian part in range?
     double tcentral_low = tmin;
     double tcentral_high = (central_high-mean)*invwidth;     
     central = rootPiBy2*width*(TMath::Erf(tcentral_high*invRoot2)-TMath::Erf(tcentral_low*invRoot2));
   }
 
 
   //compute right tail;
   if (isfullrange && (n2-1.0)>1.e-5) {
     right = width*vdt::fast_exp(-0.5*alpha2*alpha2)*n2*vdt::fast_inv(alpha2*(n2-1.));
   }
   else {    
    double right_low=std::max(xmin,mean + alpha2*width);
    double right_high=xmax;
    double tlow = (right_low - mean)*invwidth;
    
    if(right_low < right_high){ //is the right tail in range?
      double n2invalpha2 = n2*vdt::fast_inv(fabs(alpha2)); 
      if(fabs(n2-1.0)>1.e-5) {
	double invn2m2 = vdt::fast_inv(n2-1.);
	double rightpow = gbrmath::fast_pow(n2invalpha2,-n2*invn2m2);
	double right0 = width*vdt::fast_exp(-0.5*alpha2*alpha2)*invn2m2;
	double right1, right2;
	
	if (xmin<(mean+alpha2*width)) right1 = n2invalpha2;
	else right1 = gbrmath::fast_pow( rightpow*(n2invalpha2 - alpha2 + tlow), 1.-n2);
	
	if (tmax>1000.) right2 = 0.;
	else right2 = gbrmath::fast_pow( rightpow*(n2invalpha2 - alpha2 + tmax), 1.-n2);
	
	right = right0*(right1-right2);	
	
	//right = A2*vdt::fast_inv(-n2+1.0)*width*(gbrmath::fast_pow(B2+(right_high-mean)*invwidth,-n2+1.)-gbrmath::fast_pow(B2+(right_low-mean)*invwidth,-n2+1.));
      }
      else {
	double A2 = gbrmath::fast_pow(n2invalpha2,n2)*vdt::fast_exp(-0.5*alpha2*alpha2);
	double B2 = n2invalpha2-fabs(alpha2);
	right = A2*width*(vdt::fast_log(B2+(right_high-mean)*invwidth) - vdt::fast_log(B2+(right_low-mean)*invwidth) );
      }
    }
   }
   
   double sum = left + central + right;
   
   //if (!std::isnormal(left) || !std::isnormal(central) || !std::isnormal(right)) {
   //if (!std::isnormal(sum)) {
   //  printf("bad int: mean = %5f, sigma = %5f, alpha2 = %5f, n2 = %5f\n",double(mean),double(width),double(alpha2), double(n2));
   //  //printf("left = %5f, central = %5f, right = %5f, A1 = %5f, B1 = %5f, A2 = %5f, B2 = %5f, integral = %5f\n",left,central,right,A1,B1,A2,B2,left+central+right);
   //  printf("left = %5f, central = %5f, right = %5f, integral = %5f\n",left,central,right,sum);
   //}
     
   return sum;
 
 }
 
