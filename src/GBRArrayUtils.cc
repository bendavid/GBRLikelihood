#include "../interface/GBRArrayUtils.h" 
#include <limits>
    
void GBRArrayUtils::InitArrays(int *__restrict__ ns, double *__restrict__ tgts, double *__restrict__ tgt2s, float *__restrict__ bsepgains, const int nbins) {
 
  ns = (int*)__builtin_assume_aligned(ns,32);
  tgts = (double*)__builtin_assume_aligned(tgts,32);
  tgt2s = (double*)__builtin_assume_aligned(tgt2s,32);
  bsepgains = (float*)__builtin_assume_aligned(bsepgains,32);
  
  for (int ibin=0; ibin<nbins; ++ibin) {
    ns[ibin] = 0;
    tgts[ibin] = 0.;
    tgt2s[ibin] = 0.;     
    
    bsepgains[ibin] = -std::numeric_limits<float>::max();
  }
   
}
  