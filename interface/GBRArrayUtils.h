#ifndef GBRARRAYUTILS
#define GBRARRAYUTILS

class GBRArrayUtils {
  
  friend class RooHybridBDTAutoPdf;
  friend class MCGBRIntegrator;
  
protected:
 static void InitArrays(int *ns, double *tgts, double *tgt2s, float *bsepgains, const int nbins);
 static void InitArrays(int *ns, double *tgts, double *tgt2s, double *bsepgains, const int nbins);
 static void ZeroArray(double *wscls, const int nbins);
 static void MaxArray(double *wscls, const int nbins);
 static void MinArray(double *wscls, const int nbins);
 static void MaxArray(int *wscls, const int nbins);
 static void MinArray(int *wscls, const int nbins); 
 static void MinMaxQuants(int &minquant, int &maxquant, const int *quants, const int nev);
 static void FillBinQuants(int *binquants, const unsigned int offset, const unsigned int pscale, const unsigned int nquantiles, const unsigned int nbins);
 static void FillSepGains(const double *sumtgts, const double *sumtgt2s, float *bsepgains, const double fulldiff, const double sumtgt, const double sumtgt2, const int nbins);
 static void FillSepGainsMC(const double *sumtgts, const double *sumtgt2s, const double *sumfmins, const double *sumfminsr, const double *sumws, float *bsepgains, const double curval, const double sumtgt, const double sumtgt2, const double sumw, const int nbins, const double shrinkage);
 static void FillSepGainsMCEnvelope(const double *sumtgts, const double *sumtgt2s, const double *sumtgt3s, const double *sumtgt4s, const double *sumtgtmaxs, const double *sumtgtmaxsr, const double *sumfmins, const double *sumfminsr, const double *sumws, float *bsepgains, const double curval, const double sumtgt, const double sumtgt2, const double sumtgt3, const double sumtgt4, const double sumtgtmax, const double sumtgtmin, const double sumw, const int nbins);
 static void FillSepGainsMCMatrix(const double *sumtgts, const double *sumtgt2s, const double *sumtgt3s, const double *sumtgt4s, const double *sumtgt5s, const double *sumws, float *bsepgains, const double curval, const double sumtgt, const double sumtgt2, const double sumtgt3, const double sumtgt4, const double sumtgt5, const double sumw, const int nbins);


};
 
#endif
