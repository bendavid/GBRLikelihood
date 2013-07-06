#ifndef GBRARRAYUTILS
#define GBRARRAYUTILS

class GBRArrayUtils {
  
  friend class RooHybridBDTAutoPdf;
  
protected:
 static void InitArrays(int *ns, double *tgts, double *tgt2s, float *bsepgains, const int nbins);
  
};
 
#endif
