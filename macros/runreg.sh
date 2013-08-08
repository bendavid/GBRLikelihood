#!/bin/bash
root -l -q $CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/macros/eregtraining.C+\(true,false\)
root -l -q $CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/macros/eregtraining.C+\(false,false\)
root -l -q $CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/macros/eregtraining.C+\(true,true\)
root -l -q $CMSSW_BASE/src/HiggsAnalysis/GBRLikelihood/macros/eregtraining.C+\(false,true\)
