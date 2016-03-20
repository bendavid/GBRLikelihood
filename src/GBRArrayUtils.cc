#include "../interface/GBRArrayUtils.h" 
#include "../interface/GBRMath.h"
#include <limits>
#include <algorithm>    
#include "Eigen/Dense"
    
void GBRArrayUtils::InitArrays(int *__restrict__ ns, double *__restrict__ tgts, double *__restrict__ tgt2s, float *__restrict__ bsepgains, const int nbins) {
 
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  ns = (int*)__builtin_assume_aligned(ns,32);
  tgts = (double*)__builtin_assume_aligned(tgts,32);
  tgt2s = (double*)__builtin_assume_aligned(tgt2s,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    ns[ibin] = 0;
    tgts[ibin] = 0.;
    tgt2s[ibin] = 0.;     
    
    bsepgains[ibin] = -std::numeric_limits<float>::max();
  }
   
}

void GBRArrayUtils::InitArrays(int *__restrict__ ns, double *__restrict__ tgts, double *__restrict__ tgt2s, double *__restrict__ bsepgains, const int nbins) {
 
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  ns = (int*)__builtin_assume_aligned(ns,32);
  tgts = (double*)__builtin_assume_aligned(tgts,32);
  tgt2s = (double*)__builtin_assume_aligned(tgt2s,32);
  bsepgains = (double*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    ns[ibin] = 0;
    tgts[ibin] = 0.;
    tgt2s[ibin] = 0.;     
    
    bsepgains[ibin] = -std::numeric_limits<double>::max();
  }
   
}
    
void GBRArrayUtils::ZeroArray(double *__restrict__ wscls, const int nbins) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  wscls = (double*)__builtin_assume_aligned(wscls,32);
#endif
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    wscls[ibin] = 0.;
  }
}

void GBRArrayUtils::MaxArray(double *__restrict__ wscls, const int nbins) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  wscls = (double*)__builtin_assume_aligned(wscls,32);
#endif
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    wscls[ibin] = std::numeric_limits<double>::max();
  }
}

void GBRArrayUtils::MinArray(double *__restrict__ wscls, const int nbins) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  wscls = (double*)__builtin_assume_aligned(wscls,32);
#endif
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    wscls[ibin] = -std::numeric_limits<double>::max();
  }
}

void GBRArrayUtils::MaxArray(int *__restrict__ wscls, const int nbins) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  wscls = (int*)__builtin_assume_aligned(wscls,32);
#endif
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    wscls[ibin] = std::numeric_limits<int>::max();
  }
}

void GBRArrayUtils::MinArray(int *__restrict__ wscls, const int nbins) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)
  wscls = (int*)__builtin_assume_aligned(wscls,32);
#endif
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    wscls[ibin] = std::numeric_limits<int>::lowest();
  }
}
  
void GBRArrayUtils::MinMaxQuants(int &__restrict__ minquant, int &__restrict__ maxquant, const int *__restrict__ quants, const int nev) {
  
#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  quants = (const int*)__builtin_assume_aligned(quants,32);
#endif  
  
  minquant = std::numeric_limits<int>::max();
  maxquant = 0;
   
  for (int iev = 0; iev<nev; ++iev) {
    if (quants[iev]<minquant) minquant = quants[iev];
    if (quants[iev]>maxquant) maxquant = quants[iev];
  }      
  
} 
 
void GBRArrayUtils::FillBinQuants(int *__restrict__ binquants, const unsigned int offset, const unsigned int pscale, const unsigned int nquantiles, const unsigned int nbins) {

#if  __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)    
  binquants = (int*)__builtin_assume_aligned(binquants,32);
#endif
    
  for (unsigned int ibin=0; ibin<nbins; ++ibin) { 
    //int scaledbin
    //int quant = ((ibin+1)<<pscale) + offset - 1;
    unsigned int quant = ((ibin+1)<<pscale) + offset - 1;
    //unsigned short quant = (ibin<<pscale) + offset - 1;
    //int quant = ((ibin+1)<<pscale) + offset - 1;
    binquants[ibin] = std::min(quant, nquantiles-1);
    //binquants[ibin] = quant < nquantiles ? quant : nquantiles-1;
  }  
  
}
 
void GBRArrayUtils::FillSepGains(const double *__restrict__ sumtgts, const double *__restrict__ sumtgt2s, float *__restrict__ bsepgains, const double fulldiff, const double sumtgt, const double sumtgt2, const int nbins) {
  
#if __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  sumtgts = (const double*)__builtin_assume_aligned(sumtgts,32);
  sumtgt2s = (const double*)__builtin_assume_aligned(sumtgt2s,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
  for (int ibin=0; ibin<nbins; ++ibin) {     
        
    double leftdiff = std::min(0.,-0.5*sumtgts[ibin]*sumtgts[ibin]*vdt::fast_inv(sumtgt2s[ibin]));
    //double leftdiff = std::min(0.,-0.5*sumtgts[ibin]*sumtgts[ibin]/sumtgt2s[ibin]);
    //double leftdiff = -0.5*sumtgts[ibin]*sumtgts[ibin]/sumtgt2s[ibin];

    double righttgtsum = sumtgt - sumtgts[ibin];
    double righttgt2sum = sumtgt2 - sumtgt2s[ibin];
    
    double rightdiff = std::min(0.,-0.5*righttgtsum*righttgtsum*vdt::fast_inv(righttgt2sum));
    //double rightdiff = std::min(0.,-0.5*righttgtsum*righttgtsum/righttgt2sum);
    //double rightdiff = -0.5*righttgtsum*righttgtsum/righttgt2sum;

	  
    //weighted improvement in variance from this split     
    //bsepgains[ibin] = std::max(0.,fulldiff - leftdiff - rightdiff);
    bsepgains[ibin] = fulldiff - leftdiff - rightdiff;
    
    //float valid =  sumtgt2s[ibin]==0.;// || righttgt2sum==0. || leftdiff>0. || rightdiff>0.);

    

    
  }  
  
  
}

void GBRArrayUtils::FillSepGainsMC(const double *sumtgts, const double *sumtgt2s, const double *sumfmins, const double *sumfminsr, const double *sumws, float *bsepgains, const double curval, const double sumtgt, const double sumtgt2, const double sumw, const int nbins, const double shrinkage) {
  
#if __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  sumtgts = (const double*)__builtin_assume_aligned(sumtgts,32);
  sumws = (const double*)__builtin_assume_aligned(sumws,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
  for (int ibin=0; ibin<nbins; ++ibin) {     
        
//     double leftval = -sumtgts[ibin]*log(sumtgts[ibin]/(sumtgts[ibin]+sumws[ibin])) - sumws[ibin]*log(sumws[ibin]/(sumtgts[ibin]+sumws[ibin]));
// 
//     
//     double righttgtsum = sumtgt - sumtgts[ibin];
//     double rightwsum = sumw - sumws[ibin];
//     
//     
//     double rightval = -righttgtsum*log(righttgtsum/(righttgtsum+rightwsum)) - rightwsum*log(rightwsum/(righttgtsum+rightwsum));

    double righttgtsum = sumtgt - sumtgts[ibin];
    double rightwsum = sumw - sumws[ibin];
    double righttgt2sum = sumtgt2 - sumtgt2s[ibin];
    
//     double leftresponse = std::max(shrinkage*sumtgts[ibin]/sumws[ibin],-sumfmins[ibin]);
//     double rightresponse = std::max(shrinkage*righttgtsum/rightwsum,-sumfminsr[ibin]);
   
//     double leftval = sumtgt2s[ibin] + sumws[ibin]*leftresponse*leftresponse - 2.*leftresponse*sumtgts[ibin];
//     double rightval = righttgt2sum + rightwsum*rightresponse*rightresponse - 2.*rightresponse*righttgtsum;    
    
//     printf("leftresponse = %5f, rightresponse = %5f, curval = %5f, leftval = %5f, rightval = %5f, sepgain = %5f\n",leftresponse,rightresponse,curval,leftval,rightval,curval-leftval-rightval);
    
//     double righttgt2sum = sumtgt3s[ibin];
    
//     double leftval = sumtgt2s[ibin]*sumws[ibin]/sumtgts[ibin];
//     double rightval = righttgt2sum*rightwsum/righttgtsum;
    
//     double leftval = -sumtgts[ibin]/(sumws[ibin]*(sumtgt2s[ibin]+righttgt2sum));
//     double rightval = -righttgtsum/(rightwsum*(sumtgt2s[ibin]+righttgt2sum));
    
//     double leftval = sumtgt2s[ibin]*sumws[ibin] - sumtgts[ibin];
//     double rightval = righttgt2sum*rightwsum - righttgtsum;
    
//     double leftval = sumws[ibin]*sumtgt2s[ibin];
//     double rightval = rightwsum*righttgt2sum;

    double leftval = sumtgt2s[ibin] - sumtgts[ibin]*sumtgts[ibin]/sumws[ibin];
    double rightval = righttgt2sum - righttgtsum*righttgtsum/rightwsum;
    
	  
    //weighted improvement in variance from this split     
    //bsepgains[ibin] = std::max(0.,fulldiff - leftdiff - rightdiff);
    bsepgains[ibin] = curval - leftval - rightval;
    
//     printf("curval = %5f, leftval = %5f, rightval = %5f, gain = %5f\n",curval,leftval,rightval,curval-leftval-rightval);
    
    //float valid =  sumtgt2s[ibin]==0.;// || righttgt2sum==0. || leftdiff>0. || rightdiff>0.);

    

    
  }  
  
}

void GBRArrayUtils::FillSepGainsMCEnvelope(const double *sumtgts, const double *sumtgt2s, const double *sumtgt3s, const double *sumtgt4s, const double *sumtgtmaxs, const double *sumtgtmaxsr, const double *sumfmins, const double *sumfminsr, const double *sumws, float *bsepgains, const double curval, const double sumtgt, const double sumtgt2, const double sumtgt3, const double sumtgt4, const double sumtgtmax, const double sumtgtmin, const double sumw, const int nbins) {
  
#if __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  sumtgts = (const double*)__builtin_assume_aligned(sumtgts,32);
  sumws = (const double*)__builtin_assume_aligned(sumws,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif  
  
//   double dh = sumtgtmax;
//   double expdf = (sumtgt + sumw*dh)/sumtgt2;
  
//   double df = sumtgtmax;
//   double expdf = vdt::fast_exp(df);
//   double dh = (expdf*sumtgt2-sumtgt)/sumw;  
  
//   const double sigma = 1e-4;
//   const double k = 0.5/sigma/sigma;
  const double k = 10.;
  
  for (int ibin=0; ibin<nbins; ++ibin) {     
        

//     double lefttgtsum = sumtgts[ibin];
    double leftwsum = sumws[ibin];
//     double lefttgt2sum = sumtgt2s[ibin];    
//     double lefttgt3sum = sumtgt3s[ibin];    
//     double lefttgt4sum = sumtgt4s[ibin];    
    double lefttgtmax = sumtgtmaxs[ibin];
    double lefttgtmin = sumfmins[ibin];
    
//     double righttgtsum = sumtgt - sumtgts[ibin];
    double rightwsum = sumw - sumws[ibin];
//     double righttgt2sum = sumtgt2 - sumtgt2s[ibin];
//     double righttgt3sum = sumtgt3 - sumtgt3s[ibin];
//     double righttgt4sum = sumtgt4 - sumtgt4s[ibin];
    double righttgtmax = sumtgtmaxsr[ibin];
    double righttgtmin = sumfminsr[ibin];
    
//     double leftval = sumws[ibin]*sumtgtmaxs[ibin];
//     double rightval = rightwsum*sumtgtmaxsr[ibin];
    
//     double leftval = -sumtgts[ibin]/(sumws[ibin]*sumtgtmaxs[ibin]);
//     double rightval = -righttgtsum/(rightwsum*sumtgtmaxsr[ibin]);
   
//     double leftval = sumtgtmaxs[ibin]*sumtgt2s[ibin];
//     double rightval = sumtgtmaxsr[ibin]*righttgt2sum;
    
//     double leftval = (sumtgtmaxs[ibin] - sumfmins[ibin])*sumtgt2s[ibin];
//     double rightval = (sumtgtmaxsr[ibin] - sumfminsr[ibin])*righttgt2sum;

//     double leftval = sumws[ibin]*sumtgtmaxs[ibin];
//     double rightval = rightwsum*sumtgtmaxsr[ibin];
//     double leftval = sumtgtmaxs[ibin];
//     double rightval = sumtgtmaxsr[ibin];    
    
    
//     double dhleft = lefttgtmax;
//     double dsleft = (lefttgtsum - 8.*lefttgt2sum)/(8.*leftwsum);
//     double leftval = 8.*dsleft*lefttgt2sum + 4.*leftwsum*dsleft*dsleft - dsleft*lefttgtsum + dhleft*lefttgt2sum + leftwsum*dhleft*dsleft + 0.5*leftwsum*dhleft;
//     double leftval = 0.5*leftwsum*dhleft;

//     double dhright = righttgtmax;
//     double dsright = (righttgtsum - 8.*righttgt2sum)/(8.*rightwsum);
//     double rightval = 8.*dsright*righttgt2sum + 4.*rightwsum*dsright*dsright - dsright*righttgtsum + dhright*righttgt2sum + rightwsum*dhright*dsright + 0.5*rightwsum*dhright;
//     double rightval = 0.5*rightwsum*dhright;
    
//     double leftval = 2.*leftwsum*lefttgtmax - 0.5*vdt::fast_exp(-lefttgtmax)*lefttgt2sum;
//     double rightval = 2.*rightwsum*righttgtmax - 0.5*vdt::fast_exp(-righttgtmax)*righttgt2sum;
    
/*    double leftval = -vdt::fast_exp(-lefttgtmax)*lefttgt2sum;
    double rightval = -vdt::fast_exp(-righttgtmax)*righttgt2sum;  */  

//     double dhleft = lefttgtmax;
//     double dfleft = log( leftwsum / (dhleft*lefttgt2sum - lefttgtsum) );
//     double leftval = leftwsum*(1.-dfleft);
// 
//     double dhright = righttgtmax;
//     double dfright = log( rightwsum / (dhright*righttgt2sum - righttgtsum) );
//     double rightval = rightwsum*(1.-dfright);    
    
  
//     double dhleft = lefttgtmax;
//     double expndfleft = lefttgt2sum/(leftwsum*dhleft+lefttgt3sum);
//     double leftval = -expndfleft*lefttgtsum;    
//     
//     double dhright = righttgtmax;
//     double expndfright = righttgt2sum/(rightwsum*dhright+righttgt3sum);
//     double rightval = -expndfright*righttgtsum;      
    
//     -vdt::fast_exp(-tgtmax)*sumtgt2;
   
//     double leftval = -vdt::fast_exp(-lefttgtmax)*lefttgt2sum;
//     double rightval = -vdt::fast_exp(-righttgtmax)*righttgt2sum;


//     double dhleft = lefttgtmax;
//     double dlleft = -lefttgtmin + dhleft;
//     double dhmdlleft = dlleft;
// //     double leftval = dhmdlleft*lefttgtsum - 0.5*dhmdlleft*dhmdlleft*lefttgt2sum - dhleft*lefttgt3sum + 0.5*dhleft*dhleft*lefttgt4sum;   
// //     double leftval = dhmdlleft*lefttgtsum - 0.5*dhmdlleft*dhmdlleft*lefttgt2sum + dhmdlleft*dhmdlleft*dhmdlleft*lefttgt3sum/3. - dhmdlleft*dhmdlleft*dhmdlleft*dhmdlleft*lefttgt4sum/4.;
//     double leftval = dhmdlleft*lefttgtsum - 0.5*dhmdlleft*dhmdlleft*lefttgt2sum;   
//    
//     double dhright = righttgtmax;
//     double dlright = -righttgtmin + dhright;
//     double dhmdlright = dlright;
// //     double rightval = dhmdlright*righttgtsum - 0.5*dhmdlright*dhmdlright*righttgt2sum - dhright*righttgt3sum + 0.5*dhright*dhright*righttgt4sum;   
// //     double rightval = dhmdlright*righttgtsum - 0.5*dhmdlright*dhmdlright*righttgt2sum + dhmdlright*dhmdlright*dhmdlright*righttgt3sum/3. - dhmdlright*dhmdlright*dhmdlright*dhmdlright*righttgt4sum/4.;
//     double rightval = dhmdlright*righttgtsum - 0.5*dhmdlright*dhmdlright*righttgt2sum;  
   
/*    double dfleft = lefttgtmax;
    double dsleft = log( -leftwsum / (vdt::fast_exp(-dfleft)*lefttgtsum - lefttgt2sum) );
    double leftval = leftwsum*(1.-dsleft);   
   
    double dfright = righttgtmax;
    double dsright = log( -rightwsum / (vdt::fast_exp(-dfright)*righttgtsum - righttgt2sum) );
    double rightval = rightwsum*(1.-dsright);*/     
    
//     double leftval = -lefttgtsum*lefttgtsum;
//     double rightval = -righttgtsum*righttgtsum;

//     double leftlambda = lefttgtsum - expdf*lefttgt2sum + leftwsum*dh;
//     leftlambda = 0.;
//     double leftf = lefttgt3sum/expdf;
//     double dfleft = lefttgtmax;
//     double expdfleft = vdt::fast_exp(dfleft);
//     double leftfnow = lefttgt3sum/expdfleft;
//     double leftval = leftfnow*leftfnow - leftlambda*leftlambda - leftf*leftf;
// 
//     double rightlambda = righttgtsum - expdf*righttgt2sum + rightwsum*dh;
//     rightlambda = 0.;
//     double rightf = righttgt3sum/expdf;
//     double dfright = righttgtmax;
//     double expdfright = vdt::fast_exp(dfright);
//     double rightfnow = righttgt3sum/expdfright;
//     double rightval = rightfnow*rightfnow - rightlambda*rightlambda - rightf*rightf;
    
/*    double dh = sumtgtmax;
    
    double dhleft = lefttgtmax;
    double leftbefore = lefttgtsum - leftwsum*dh;
    double leftafter = lefttgtsum - leftwsum*dhleft;
    double leftval = leftafter*leftafter - leftbefore*leftbefore;

    double dhright = righttgtmax;
    double rightbefore = righttgtsum - rightwsum*dh;
    double rightafter = righttgtsum - rightwsum*dhright;
    double rightval = rightafter*rightafter - rightbefore*rightbefore;  */  
    
//     double leftval = pow(lefttgtsum - leftwsum*lefttgtmax

//     printf("ibin = %i, leftval = %5e, rightval = %5e, curval = %5e, diff = %5e\n",ibin,leftval,rightval,curval,curval-leftval-rightval);
    
//     printf("ibin = %i, lefttgtsum = %5f, lefttgt2sum = %5f\n",ibin,lefttgtsum,lefttgt2sum);
//     printf("ibin = %i, curval = %5e, leftsumw = %5f, lefttgtmax = %5f, rightsumw = %5f, righttgtmax = %5f, leftval = %5e, rightval = %5e, sepgain = %5f, lefttgtsum = %5f, lefttgt2sum = %5f, righttgtsum = %5f, righttgt2sum = %5f\n",ibin,curval,sumws[ibin],sumtgtmaxs[ibin],rightwsum,sumtgtmaxsr[ibin],leftval,rightval, curval-leftval-rightval,lefttgtsum,lefttgt2sum,righttgtsum,righttgtsum);
    
//     bsepgains[ibin] = curval - leftval - rightval;


//     double sepgain = (lefttgtmax>righttgtmax) ? (righttgtmax-sumtgtmax)*( : 2.*leftwsum*leftwsum;

//     double leftval = -(lefttgtmax-sumtgtmax)*(lefttgtmax-sumtgtmax);
//     double rightval = -(righttgtmax-sumtgtmax)*(righttgtmax-sumtgtmax);
//     double sepgain = -leftval -rightval;

//     double sepgain = (lefttgtmax > righttgtmax) ? 2.*righttgtsum*righttgtsum : 2.*lefttgtsum*lefttgtsum;

//     printf("ibin = %i, lefttgtmax = %5f, righttgtmax = %5f, sumtgtmax = %5f, lefttgtsum = %5f, righttgtsum = %5f, sepgain = %5f\n",ibin,lefttgtmax,righttgtmax,sumtgtmax,lefttgtsum,righttgtsum,sepgain);

/*    double dhleft = lefttgtsum/leftwsum;
    double leftval = -2.*dhleft*lefttgtsum + dhleft*dhleft - vdt::fast_exp(-lefttgtmin)*lefttgt2sum;

    double dhright = righttgtsum/rightwsum;
    double rightval = -2.*dhright*righttgtsum + dhright*dhright - vdt::fast_exp(-righttgtmin)*righttgt2sum; */   

//     double leftval = - vdt::fast_exp(-lefttgtmin)*lefttgt2sum + lefttgt4sum - lefttgtmax*lefttgt3sum;;
//     double rightval = - vdt::fast_exp(-righttgtmin)*righttgt2sum + righttgt4sum - righttgtmax*righttgt3sum;;

//     double leftval = - vdt::fast_exp(-lefttgtmin)*lefttgt2sum;
//     double rightval = - vdt::fast_exp(-righttgtmin)*righttgt2sum;
    
//     double leftval = -0.5*lefttgtsum*lefttgtsum/lefttgt2sum;
//     double rightval = -0.5*righttgtsum*righttgtsum/righttgt2sum;
    
//     double leftval = lefttgtmax - 0.*lefttgt2sum*exp(-lefttgtmax);
//     double rightval = righttgtmax - 0.*righttgt2sum*exp(-righttgtmax);
    
//     double dhleft = lefttgtmax;
//     double dwleft = lefttgtmin + dhleft;
//     double leftval = lefttgtsum*dwleft - 0.5*lefttgt2sum*dwleft*dwleft + (1./3.)*lefttgt3sum*dwleft*dwleft*dwleft;    
// 
//     double dhright = righttgtmax;
//     double dwright = righttgtmin + dhright;
//     double rightval = righttgtsum*dwright - 0.5*righttgt2sum*dwright*dwright + (1./3.)*righttgt3sum*dwright*dwright*dwright;  
    
//     const double k = 5.;
//     
//     double dhleft = lefttgtmax;
//     double leftval = lefttgtsum*dhleft - 0.5*lefttgt2sum*dhleft*dhleft + (1./3.)*lefttgt3sum*dhleft*dhleft*dhleft - 0.25*lefttgt4sum*dhleft*dhleft*dhleft*dhleft - k*(-lefttgt2sum*dhleft + lefttgt3sum*dhleft*dhleft - lefttgt4sum*dhleft*dhleft*dhleft);
//     
//     double dhright = righttgtmax;
//     double rightval = righttgtsum*dhright - 0.5*righttgt2sum*dhright*dhright + (1./3.)*righttgt3sum*dhright*dhright*dhright - 0.25*righttgt4sum*dhright*dhright*dhright*dhright - k*(-righttgt2sum*dhright + righttgt3sum*dhright*dhright - righttgt4sum*dhright*dhright*dhright); 
    
    
//     double leftval = k*(lefttgt2sum - lefttgtsum*lefttgtsum/leftwsum);
//     double rightval = k*(righttgt2sum - righttgtsum*righttgtsum/rightwsum);
    
//     double dhleft = std::max(0.,lefttgtsum/leftwsum);
//     double leftval = k*(lefttgt2sum - 2.*dhleft*lefttgtsum + leftwsum*dhleft*dhleft);    

//     double dhright = std::max(0.,righttgtsum/rightwsum);
//     double rightval = k*(righttgt2sum - 2.*dhright*righttgtsum + rightwsum*dhright*dhright);      

//     double leftval = -k*lefttgt2sum*lefttgt2sum/lefttgtsum;
//     double rightval = -k*righttgt2sum*righttgt2sum/righttgtsum;

//     double leftval = k*leftwsum*lefttgtmin + k*leftwsum*std::max(0.,lefttgtmax);
//     double rightval = k*rightwsum*righttgtmin + k*rightwsum*std::max(0.,righttgtmax);

   
//     double leftval = k*leftwsum*lefttgtmin;
//     double rightval = k*rightwsum*righttgtmin;

    double leftval = k*leftwsum*lefttgtmin + k*leftwsum*std::max(0.,lefttgtmax);
    double rightval = k*rightwsum*righttgtmin + k*rightwsum*std::max(0.,righttgtmax);

//     double leftval = k*leftwsum*lefttgtmin + k*leftwsum*lefttgtmax;
//     double rightval = k*rightwsum*righttgtmin + k*rightwsum*righttgtmax;


//     double dhleft = lefttgtsum/leftwsum;
//     double leftval = k*leftwsum*lefttgtmin + 1e4*(leftwsum*dhleft*dhleft - 2.*dhleft*lefttgtsum);
// 
//     double dhright = righttgtsum/rightwsum;
//     double rightval = k*rightwsum*righttgtmin + 1e4*(rightwsum*dhright*dhright - 2.*dhright*righttgtsum);   

/*    double dhleft = std::max(0.,lefttgtsum/leftwsum);
    double leftval = k*leftwsum*lefttgtmin + 1e4*(leftwsum*dhleft*dhleft - 2.*dhleft*lefttgtsum);

    double dhright = std::max(0.,righttgtsum/rightwsum);
    double rightval = k*rightwsum*righttgtmin + 1e4*(rightwsum*dhright*dhright - 2.*dhright*righttgtsum);  */  

//     double dhleft = std::max(0.,lefttgtmax);
//     double leftval = k*leftwsum*lefttgtmin + k*(lefttgtsum*dhleft - lefttgt2sum*dhleft*dhleft/2. + lefttgt3sum*dhleft*dhleft*dhleft/3. - lefttgt4sum*dhleft*dhleft*dhleft*dhleft/4.);
// 
//     double dhright = std::max(0.,righttgtmax);
//     double rightval = k*rightwsum*righttgtmin + k*(righttgtsum*dhright - righttgt2sum*dhright*dhright/2. + righttgt3sum*dhright*dhright*dhright/3. - righttgt4sum*dhright*dhright*dhright*dhright/4.);
    
    double sepgain = curval - leftval - rightval;

//     printf("ibin = %i, leftval = %5f, rightval = %5f, sepgain = %5f, curval = %5f\n",ibin,leftval,rightval,sepgain, curval);
    
    bsepgains[ibin] = sepgain;

    

    
  }  
  
}

void GBRArrayUtils::FillSepGainsMCMatrix(const double *sumtgts, const double *sumtgt2s, const double *sumtgt3s, const double *sumtgt4s, const double *sumtgt5s, const double *sumws, float *bsepgains, const double curval, const double sumtgt, const double sumtgt2, const double sumtgt3, const double sumtgt4, const double sumtgt5, const double sumw, const int nbins) {
  
#if __GNUC__>4 || (__GNUC__==4 && __GNUC_MINOR__>=7)  
  sumtgts = (const double*)__builtin_assume_aligned(sumtgts,32);
  sumtgt2s = (const double*)__builtin_assume_aligned(sumtgt2s,32);
  sumtgt3s = (const double*)__builtin_assume_aligned(sumtgt3s,32);
  sumtgt4s = (const double*)__builtin_assume_aligned(sumtgt4s,32);
  sumtgt5s = (const double*)__builtin_assume_aligned(sumtgt5s,32);
  sumws = (const double*)__builtin_assume_aligned(sumws,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
#endif    
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    
//     double leftwsum = sumws[ibin];
    double lefttgtsum = sumtgts[ibin];
//     double lefttgt2sum = sumtgt2s[ibin];
//     double lefttgt3sum = sumtgt3s[ibin];
//     double lefttgt4sum = sumtgt4s[ibin];
//     double lefttgt5sum = sumtgt5s[ibin]; 
    
//     double rightwsum = sumw - sumws[ibin];
    double righttgtsum = sumtgt - sumtgts[ibin];
//     double righttgt2sum = sumtgt2 - sumtgt2s[ibin];
//     double righttgt3sum = sumtgt3 - sumtgt3s[ibin];
//     double righttgt4sum = sumtgt4 - sumtgt4s[ibin];
//     double righttgt5sum = sumtgt5 - sumtgt5s[ibin];

    
    
    
/*    Eigen::Vector2d dfleft(sumtgts[ibin],sumtgt2s[ibin]);
    Eigen::Matrix2d d2fleft;
    d2fleft(0,0) = sumtgt3s[ibin];
    d2fleft(0,1) = sumtgt5s[ibin];
    d2fleft(1,0) = sumtgt5s[ibin];
    d2fleft(1,1) = sumtgt4s[ibin];
    
    Eigen::Vector2d dxleft = -d2fleft.ldlt().solve(dfleft);

    double leftval = -0.5*dxleft.transpose()*d2fleft*dxleft;    

    Eigen::Vector2d dfright(righttgtsum,righttgt2sum);
    Eigen::Matrix2d d2fright;
    d2fright(0,0) = righttgt3sum;
    d2fright(0,1) = righttgt5sum;
    d2fright(1,0) = righttgt5sum;
    d2fright(1,1) = righttgt4sum;
    
    Eigen::Vector2d dxright = -d2fright.ldlt().solve(dfright);

    double rightval = -0.5*dxright.transpose()*d2fright*dxright;  */  

//     double dhleft = lefttgtsum/lefttgt2sum;
//     double dsleft = 0.5*log(leftwsum/(dhleft*dhleft*lefttgt2sum-dhleft*lefttgtsum+lefttgt3sum));
//     double exp2dsleft = vdt::fast_exp(2.*dsleft);
//     double leftval = -leftwsum*dsleft - 0.5*dhleft*exp2dsleft*lefttgtsum + 0.5*dhleft*dhleft*exp2dsleft*lefttgt2sum + 0.5*exp2dsleft*lefttgt3sum;
// 
//     double dhright = righttgtsum/righttgt2sum;
//     double dsright = 0.5*log(rightwsum/(dhright*dhright*righttgt2sum-dhright*righttgtsum+righttgt3sum));
//     double exp2dsright = vdt::fast_exp(2.*dsright);
//     double rightval = -rightwsum*dsright - 0.5*dhright*exp2dsright*righttgtsum + 0.5*dhright*dhright*exp2dsright*righttgt2sum + 0.5*exp2dsright*righttgt3sum;
    
//     double leftval = lefttgt2sum - lefttgtsum*lefttgtsum/leftwsum;
//     double rightval = righttgt2sum - righttgtsum*righttgtsum/rightwsum;
    
//     double dsleft = 0.5*log(leftwsum/(dhleft*dhleft*lefttgt2sum-dhleft*lefttgtsum+lefttgt3sum));
//     double exp2dsleft = vdt::fast_exp(2.*dsleft);
//     double leftval = -leftwsum*dsleft - 0.5*dhleft*exp2dsleft*lefttgtsum + 0.5*dhleft*dhleft*exp2dsleft*lefttgt2sum + 0.5*exp2dsleft*lefttgt3sum;    
    
/*    double dsleft = -0.5*log(lefttgt3sum/leftwsum);
    double leftval = leftwsum*(0.5-dsleft);

    double dsright = -0.5*log(righttgt3sum/rightwsum);
    double rightval = rightwsum*(0.5-dsright); */   

//     double leftval = lefttgt2sum - lefttgtsum*lefttgtsum/leftwsum;
//     double rightval = righttgt2sum - righttgtsum*righttgtsum/rightwsum;

    double leftval = -lefttgtsum*lefttgtsum;
    double rightval = -righttgtsum*righttgtsum;
    
    bsepgains[ibin] = curval - leftval - rightval;
    
    printf("ibin = %i, curval = %5e, leftval = %5e, rightval = %5e, sepgain = %5f\n",ibin,curval,leftval,rightval, curval-leftval-rightval);

    
//     printf("ibin = %i, curval = %5e, leftval = %5e, rightval = %5e, sepgain = %5f, lefth = %5f, lefts = %5f, righth = %5f, rights = %5f\n",ibin,curval,leftval,rightval, curval-leftval-rightval,dxleft[0],dxleft[1],dxright[0],dxright[1]);
    
  }
  
}
