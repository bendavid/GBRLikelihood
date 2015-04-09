#include "../interface/WMassFitter.h"

#include "Math/Factory.h"
#include "../interface/GBRMath.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "omp.h"


WMassFitter::WMassFitter(std::vector<RooDataSet*> &indata, const char *pdfname) : pdfname_(pdfname) {
 
  kpdfmean_.fill(0.);
    
  const unsigned int imwup = 176;
  const unsigned int imwdown = 26;
  const unsigned int ipdfnom = 309;      
  
  
  for (unsigned int idata=0; idata<indata.size(); ++idata) {
    RooDataSet *data = indata[idata];
    
    printf("idata = %i, nev = %i, sumw = %5f\n",int(idata),data->numEntries(),data->sumEntries());
    
    data_.emplace_back(data->numEntries());
    WMassDataset &outdset = data_.back();
    
    TString sample;
    if (idata==0 || idata==1) {
      sample = "wplus";
    }
    else if (idata==2 || idata==3) {
      sample = "wminus";
    }
        
    const RooArgSet *dset = data->get();
    
    RooAbsReal *mwupreal = static_cast<RooAbsReal*>(dset->find(TString::Format("%s_%s_%i",pdfname_.c_str(),sample.Data(),imwup)));
    RooAbsReal *mwdownreal = static_cast<RooAbsReal*>(dset->find(TString::Format("%s_%s_%i",pdfname_.c_str(),sample.Data(),imwdown)));
    
    RooAbsReal *pdfnomreal = static_cast<RooAbsReal*>(dset->find(TString::Format("%s_%s_%i",pdfname_.c_str(),sample.Data(),ipdfnom)));
    
    std::vector<RooAbsReal*> pdfupreal(npdferr_);
    std::vector<RooAbsReal*> pdfdownreal(npdferr_);
    
    for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
      unsigned int ipdfup = ipdfnom + 2*ipdferr + 1;
      unsigned int ipdfdown = ipdfnom + 2*ipdferr + 2;      
      
      pdfupreal[ipdferr] = static_cast<RooAbsReal*>(dset->find(TString::Format("%s_%s_%i",pdfname_.c_str(),sample.Data(),ipdfup)));
      pdfdownreal[ipdferr] = static_cast<RooAbsReal*>(dset->find(TString::Format("%s_%s_%i",pdfname_.c_str(),sample.Data(),ipdfdown)));
    }

    
    for (int ievent=0; ievent<data->numEntries(); ++ievent) {
      data->get(ievent);
      double weight = data->weight();
      outdset.weight_[ievent] = weight;
      outdset.origweight_[ievent] = weight;
      outdset.prmwup_[ievent] = mwupreal->getVal();
      outdset.prmwdown_[ievent] = mwdownreal->getVal();
      for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
        outdset.pdfkappaup_[ievent][ipdferr] = vdt::fast_log(pdfupreal[ipdferr]->getVal()/pdfnomreal->getVal());
        outdset.pdfkappadown_[ievent][ipdferr] = vdt::fast_log(pdfdownreal[ipdferr]->getVal()/pdfnomreal->getVal());        
      }
      
    }
  }
  
}

void WMassFitter::Fit() {
  
  WMassMinimFunc *minimfunc = new WMassMinimFunc(this);

  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
//   ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  
  
//   ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
  min->SetFunction(*minimfunc);
  min->SetStrategy(2);
//   min->SetErrorDef(0.5);
//   min->SetTolerance(50.);
//   min->SetTolerance(1.);
//   min->SetPrecision(1e-21);
  min->SetPrecision(1e-28);
//   min->SetTolerance(1.);
  min->SetTolerance(1.);
//   min->SetTolerance(1.);
  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(10000);
  min->SetPrintLevel(9);
  
//   min->SetLimitedVariable(0,"mw",80.398,0.03,80.200,80.596);
//   min->SetLimitedVariable(1,"kplus",1.0,0.01,0.5,2.0);
//   min->SetLimitedVariable(2,"kminus",1.0,0.01,0.5,2.0);
  min->SetVariable(0,"mw",80.398,0.03);
  min->SetVariable(1,"kplus",1.0,0.01);
  min->SetVariable(2,"kminus",1.0,0.01);  
  for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
    min->SetVariable(3+ipdferr,TString::Format("kpdf_%i",ipdferr).Data(),0.,1.);
//     min->SetLimitedVariable(3+ipdferr,TString::Format("kpdf_%i",ipdferr).Data(),0.,1.,-5.,5.);
  }
  
//   for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
//     min->FixVariable(3+ipdferr);
//   }  
//   min->FixVariable(2);
  
/*  RooRealVar *mwfit = new RooRealVar("mwfit","",80.398);
  RooDataSet *fitres = new RooDataSet("fitres","",*mwfit);
  
  RooRealVar *mwmean = new RooRealVar("mwmean","",80.398,80.2,80.596);
  RooRealVar *mwsigma = new RooRealVar("mwsigma","",0.01,1e-6,1.0);
  
  RooGaussian *mwgaus = new RooGaussian("mwgaus","",*mwfit,*mwmean,*mwsigma);
  
  const unsigned int ntoy = 1000;
//   const unsigned int ntoy = 1;
  for (unsigned int itoy=0; itoy<ntoy; ++itoy) {
    //randomize nuisances
    for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
      double krand = gRandom->Gaus();
      min->SetVariableValue(3+ipdferr,krand);
      kpdfmean_[ipdferr] = krand;
//       min->ReleaseVariable(3+ipdferr);
//       min->FixVariable(3+ipdferr);
    }
    min->SetVariableValue(0,80.398);
    min->SetVariableValue(1,1.0);
    min->SetVariableValue(2,1.0);
    min->SetStrategy(0);
    
    min->Minimize();
    double mwfitfull = min->X()[0];

    mwfit->setVal(mwfitfull);
    fitres->add(*mwfit);
  }
  
  mwgaus->fitTo(*fitres);
  
  new TCanvas;
  RooPlot *plot = mwfit->frame(80.2,80.596,100);
  fitres->plotOn(plot);
  mwgaus->plotOn(plot);
  plot->Draw();
  return;*/    
  
  
  RooRealVar *mwfit = new RooRealVar("mwfit","",80.398);
  RooDataSet *fitres = new RooDataSet("fitres","",*mwfit);
  
  RooRealVar *mwmean = new RooRealVar("mwmean","",80.398,80.2,80.596);
  RooRealVar *mwsigma = new RooRealVar("mwsigma","",0.01,1e-6,1.0);
  
  RooGaussian *mwgaus = new RooGaussian("mwgaus","",*mwfit,*mwmean,*mwsigma);
  
  const unsigned int ntoy = 100;
//   const unsigned int ntoy = 1;
  for (unsigned int itoy=0; itoy<ntoy; ++itoy) {
    //randomize nuisances
    for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
      double krand = gRandom->Gaus();
      min->SetVariableValue(3+ipdferr,krand);
      kpdfmean_[ipdferr] = krand;
//       min->ReleaseVariable(3+ipdferr);
//       min->FixVariable(3+ipdferr);
    }
    min->SetVariableValue(0,80.398);
    min->SetVariableValue(1,1.0);
    min->SetVariableValue(2,1.0);
    min->SetStrategy(0.);
    
    //randomize data (bootstrapping)
    const unsigned int ndata = data_.size();
    for (unsigned int idata=0; idata<ndata; ++idata) {
      if (idata%2 != 0) continue;
      WMassDataset &data = data_[idata];
      unsigned int nev = data.weight_.size();
      for (unsigned int iev=0; iev<nev; ++iev) {
        data.weight_[iev] = gRandom->Poisson(data.origweight_[iev]);
      } 
      
    }
    
    min->Minimize();
    double mwfitfull = min->X()[0];
    
//     for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
//       double krand = 0.;
//       min->SetVariableValue(3+ipdferr,krand);
//       kpdfmean_[ipdferr] = krand;
//       min->FixVariable(3+ipdferr);
//     }
//     min->SetVariableValue(0,80.398);
//     min->SetVariableValue(1,1.0);
//     min->SetVariableValue(2,1.0);    
//     min->Minimize();
//     double mwfitnosyst = min->X()[0];

    
//     mwfit->setVal(min->X()[0]);
//     mwfit->setVal(mwfitfull - mwfitnosyst + 80.398);
    mwfit->setVal(mwfitfull);
    fitres->add(*mwfit);
  }
  
  mwgaus->fitTo(*fitres);
//   
  new TCanvas;
  RooPlot *plot = mwfit->frame(80.2,80.596,100);
  fitres->plotOn(plot);
  mwgaus->plotOn(plot);
  plot->Draw();
  return;  

//   min->FixVariable(4);  
//   min->FixVariable(0);  
  
//   min->SetStrategy(0.);  
//   min->Minimize();
//   min->SetStrategy(1.);
//   min->Minimize();
//   min->SetStrategy(2.);

  min->Minimize();

  return;
  
  min->SetStrategy(0);
  min->SetTolerance(1.);
  
//   return;
  
/*  double uperr;
  double downerr;
  min->GetMinosError(0,uperr,downerr);*/  
  
//   return;
// //   
//   min->FixVariable(0);
//   min->Minimize();
//   
//   return;
  
  double mw0 = min->X()[0];
  
//   double errtgt = 0.1;
  double mwerr = min->Errors()[0];
  mwerr = 0.022;
  
  min->FixVariable(0);
//   for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
//     min->FixVariable(3+ipdferr);
//   }
  
  double errscale = 1.2;
  double mwstart = mw0-errscale*mwerr;
  const unsigned int nstep = 19;
  double stepsize = 2.*errscale*mwerr/double(nstep-1);  
  
  
  std::vector<double> nllvals(nstep);
  std::vector<double> mwvals(nstep);
  double nllmin = 0.;
  for (unsigned int istep=0; istep<nstep; ++istep) {
//   for (unsigned int istep=0; istep<1; ++istep) {    
    double mwval = mwstart + istep*stepsize;
    mwvals[istep] = mwval;
    min->SetVariableValue(0,mwval);
    
//     min->SetStrategy(0.);
//     min->Minimize();
//     min->SetStrategy(1.);
//     min->Minimize();
//     min->SetStrategy(2.);
    min->Minimize();    
    
    nllvals[istep] = min->MinValue();
    
    if (istep == (nstep/2)) {
      nllmin = min->MinValue();
    }
  }
  
  const unsigned int nhalfstep = nstep/2;
  
  TGraph *hnll = new TGraph(nstep);
  TGraph *hlow = new TGraph(nhalfstep);
  TGraph *hhigh = new TGraph(nhalfstep);
  for (unsigned int istep=0; istep<nstep; ++istep) {
    hnll->SetPoint(istep,mwvals[istep],nllvals[istep]-nllmin);
  }
  
  for (unsigned int ipoint=0; ipoint<nhalfstep; ++ipoint) {
    int istep = nstep/2-1-ipoint;
    hlow->SetPoint(ipoint,nllvals[istep]-nllmin,mwvals[istep]);
  }
  for (unsigned int ipoint=0; ipoint<nhalfstep; ++ipoint) {
    int istep = nstep/2+1+ipoint;
    hhigh->SetPoint(ipoint,nllvals[istep]-nllmin,mwvals[istep]);
  }
  
  double mwlow = hlow->Eval(1.0,0,"S");
  double mwhigh = hhigh->Eval(1.0,0,"S");
  
  double mwerrlow = mw0-mwlow;
  double mwerrhigh = mwhigh-mw0;
  
  printf("mw0 = %5f, mwlow = %5f, mwhigh = %5f\n",mw0,mwlow,mwhigh);
  printf("mw = %5f +%5f -%5f\n",mw0,mwerrhigh,mwerrlow);
  
  new TCanvas;
  hnll->Draw("ALP");

  new TCanvas;
  hhigh->Draw("ALP");
  
  new TCanvas;
  hlow->Draw("ALP");  
  
  return;
  
  
  
//   double uperr;
//   double downerr;
//   min->GetMinosError(0,uperr,downerr);  
//   return;
  
//   for (unsigned int ipdferr=0; ipdferr<npdferr_; ++ipdferr) {
//     min->FixVariable(3+ipdferr);
//   }
  

  
  

}

double WMassFitter::WMassMinimFunc::DoEval(const double * x) const {
  
  double mw = x[0];
  double kplus = x[1];
  double kminus = x[2];
//   kminus = kplus;
  
  double frachigh = (mw-80.200)/0.396;
  double fraclow = 1.-frachigh;
  
  constexpr unsigned int npdferr = fitter_->npdferr_;
  
//   const double scale = 10.;
  
  double minval = 0.;
  
  const double *kpdfmeans = fitter_->kpdfmean_.data();
  for (unsigned int ipdferr=0; ipdferr<npdferr; ++ipdferr) {  
    double kpdf = x[3+ipdferr];
    double kmean = kpdfmeans[ipdferr];
    minval += 0.5*(kpdf-kmean)*(kpdf-kmean);
  }
  
  const unsigned int ndata = fitter_->data_.size();
  for (unsigned int idata=0; idata<ndata; ++idata) {
    const WMassDataset &data = fitter_->data_[idata];
    unsigned int nev = data.weight_.size();
    
    double k = (idata<=1) ? kplus : kminus;
  
    const double *weights = data.weight_.data();
    const double *prmwup = data.prmwup_.data();
    const double *prmwdown = data.prmwdown_.data();
  
    if (idata%2==0) {      
    //data
      #pragma omp parallel for simd reduction(+:minval)
      for (unsigned int iev=0; iev<nev; ++iev) {
        const double *pdfkappaup = data.pdfkappaup_[iev].data();
        const double *pdfkappadown = data.pdfkappadown_[iev].data();
        
        double weight = weights[iev];
        minval += -weight*vdt::fast_log(k*(frachigh*prmwup[iev] + fraclow*prmwdown[iev])); 
        double sumthetalogkappa = 0;
        for (unsigned int ipdferr=0; ipdferr<npdferr; ++ipdferr) {
          double logkappaHigh = pdfkappaup[ipdferr];
          double logkappaLow = pdfkappadown[ipdferr];
          double theta = (1./1.645)*x[3+ipdferr];
          double logKappa = LogKappa(theta,logkappaHigh,logkappaLow);
          sumthetalogkappa += theta*logKappa;
        }
        minval += -weight*sumthetalogkappa;
      }
    }
    else {
      //refdata
      #pragma omp parallel for simd reduction(+:minval)
      for (unsigned int iev=0; iev<nev; ++iev) {
        const double *pdfkappaup = data.pdfkappaup_[iev].data();
        const double *pdfkappadown = data.pdfkappadown_[iev].data();
        
        double weight = weights[iev];
        double evtprod = weight*k*(frachigh*prmwup[iev] + fraclow*prmwdown[iev]); 
        double sumexp = 0;
        for (unsigned int ipdferr=0; ipdferr<npdferr; ++ipdferr) {
          double logkappaHigh = pdfkappaup[ipdferr];
          double logkappaLow = pdfkappadown[ipdferr];
          double theta = (1./1.645)*x[3+ipdferr];
          double logKappa = LogKappa(theta,logkappaHigh,logkappaLow);
          sumexp += theta*logKappa;
        }
        evtprod *= vdt::fast_exp(sumexp);        
        minval += evtprod;
      }
    }
  }
  
  
  return 2.0*minval;
}