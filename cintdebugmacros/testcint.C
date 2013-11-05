#include "TClass.h"

void testcint() {
 
  TClass::GetClass("HybridGBRForestFlex")->Print();
  TClass::GetClass("RooSpline1D")->Print();
  
}