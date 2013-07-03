#include "../interface/HybridGBRForest.h"




//_______________________________________________________________________
HybridGBRForest::HybridGBRForest(int ntargets)
{
  fInitialResponse.resize(ntargets,0.);
  fResponses.resize(ntargets);
  fTrees.resize(ntargets);
}

//_______________________________________________________________________
HybridGBRForest::~HybridGBRForest() 
{
}
