#include "../interface/HybridGBREvent.h"

//_______________________________________________________________________
HybridGBREvent::HybridGBREvent(int nvars, int ntargets, int nparms) : 
  fVars(new float[nvars]),
  fVarsAlt(new float[nvars]),
  fTargets(new float[ntargets]),
  fTransTargets(new float[ntargets]),
  fTransTargets2(new float[ntargets]),
  fSmoothedTargets(new float[ntargets]),
  fQuantiles(new int[nvars]),
  fDerivatives(new double[nparms]),
  fDerivatives2(new double[nparms]),
  fCurrentNodes(new unsigned int[ntargets]),
  fValVector(ntargets),
  fParmMatrix(nparms),  
  fPdfVal(0.),
  fTarget(0.0),
  fTransTarget(0.0),
  fWeight(1.0),
  fWeightedTransTarget(0.),
  fWeightedTransTarget2(0.),
  fInverseGenPdf(1.0),
  fClass(0),
  fCurrentNode(0),
  fIsPrimary(true)
{

  for (int itgt=0; itgt<ntargets; ++itgt) {
    fTargets[itgt] = 0.;
    fTransTargets[itgt] = 0.;
    fTransTargets2[itgt] = 0.;
    fSmoothedTargets[itgt] = 0.;
    fCurrentNodes[itgt] = 0;
  }  
  
  for (int ivar=0; ivar<nparms; ++ivar) {
    fDerivatives[ivar] = 0.;
    fDerivatives2[ivar] = 0.;
  }   
  
}

//_______________________________________________________________________
HybridGBREvent::~HybridGBREvent() 
{
  if (fVars) delete[] fVars;
  if (fVarsAlt) delete[] fVarsAlt;  
  if (fTargets) delete[] fTargets;
  if (fTransTargets) delete[] fTransTargets;
  if (fSmoothedTargets) delete[] fTargets;  
  if (fQuantiles) delete [] fQuantiles;
  if (fDerivatives) delete [] fDerivatives;
  if (fDerivatives2) delete [] fDerivatives2;
  if (fCurrentNodes) delete [] fCurrentNodes;
  
}
