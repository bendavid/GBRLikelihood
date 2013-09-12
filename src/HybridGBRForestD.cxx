#include "../interface/HybridGBRForestD.h"




//_______________________________________________________________________
HybridGBRForestD::HybridGBRForestD(int ntargets)
{
  fInitialResponse.resize(ntargets,0.);
//  fResponses.resize(ntargets);
  fTrees.resize(ntargets);
}

//_______________________________________________________________________
HybridGBRForestD::~HybridGBRForestD() 
{
}
