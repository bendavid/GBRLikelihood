#ifndef GBRLikelihood_WMassFitter_h
#define GBRLikelihood_WMassFitter_h

#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <RooDataSet.h>
#include "Math/Minimizer.h"


class WMassFitter {

    public:
      static constexpr unsigned int npdferr_ = 26;
      
      class WMassMinimFunc : public ROOT::Math::IBaseFunctionMultiDim {
        public:
          WMassMinimFunc(WMassFitter *fitter) : fitter_(fitter) {}
          unsigned int NDim() const { return 3 + fitter_->npdferr_; }
          IBaseFunctionMultiDim * Clone() const { return new WMassMinimFunc(*this); }
       
        private:
          double DoEval(const double * x) const;
          
          inline double LogKappa(double theta, double logkappaHigh, double logkappaLow) const;
          
          
          WMassFitter *fitter_;
       
      };
      
      struct WMassDataset {
        
        WMassDataset(unsigned int nevents) :
          prmwup_(nevents),
          prmwdown_(nevents),
          pdfkappaup_(nevents),
          pdfkappadown_(nevents),
          weight_(nevents),
          origweight_(nevents)
          { }
        
        
        std::vector<double> prmwup_;
        std::vector<double> prmwdown_;
        std::vector<std::array<double,WMassFitter::npdferr_> > pdfkappaup_;
        std::vector<std::array<double,WMassFitter::npdferr_> > pdfkappadown_;
        std::vector<double> weight_;      
        std::vector<double> origweight_;
        

      };      
     
      WMassFitter(std::vector<RooDataSet*> &indata, const char *pdfname = "expfalt");

      void Fit();

    protected:
      
      std::vector<WMassDataset> data_;
      
      std::array<double,npdferr_> kpdf_;
      std::array<double,npdferr_> kpdfmean_;
      
      std::string pdfname_;



};

inline double WMassFitter::WMassMinimFunc::LogKappa(double theta, double logkappaHigh, double logkappaLow) const {
  
  if (std::abs(theta) >= 0.5) {
    return (theta >= 0 ? logkappaHigh : - logkappaLow);
  }
  else {
    double logKhi =  logkappaHigh;
    double logKlo = -logkappaLow;
    double avg = 0.5*(logKhi + logKlo), halfdiff = 0.5*(logKhi - logKlo);
    double twotheta = theta+theta, twotheta2 = twotheta*twotheta;
    double alpha = 0.125 * twotheta * (twotheta2 * (3*twotheta2 - 10.) + 15.);
    return avg + alpha*halfdiff;
  }
  return 0.;
  
}


#endif