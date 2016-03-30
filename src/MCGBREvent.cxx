#include "../interface/MCGBREvent.h"

//_______________________________________________________________________
MCGBREvent::MCGBREvent(int nvars) : 
  fVars(new float[nvars]),
  fQuantiles(new int[nvars]),
  fFuncVal(0.),
  fFuncValAlt(0.),
  fWeight(1.0),
  fTarget(0.),
  fTargetMin(0.),
  fTarget3(0.),
  fArg(0.),
  fArgLog(0.)
{
  
}

//_______________________________________________________________________
MCGBREvent::~MCGBREvent() 
{
  if (fVars) delete[] fVars;
  if (fQuantiles) delete [] fQuantiles;
  
}
