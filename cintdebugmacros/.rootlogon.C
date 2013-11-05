{

  {

  TString libstr(Form("%s/lib/%s/%s",
                      gSystem->Getenv("CMSSW_BASE"),
                      gSystem->Getenv("SCRAM_ARCH"),
                      "libHiggsAnalysisCombinedLimit.so"));

  gSystem->Load(libstr);

  }

  {

  TString libstr(Form("%s/lib/%s/%s",
                      gSystem->Getenv("CMSSW_BASE"),
                      gSystem->Getenv("SCRAM_ARCH"),
                      "libHiggsAnalysisGBRLikelihood.so"));

  gSystem->Load(libstr);

  }

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/interface");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/interface");



}