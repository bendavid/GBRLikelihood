#ifndef GBRMATH
#define GBRMATH

#include "DataFormats/Math/interface/VDTMath.h"

namespace gbrmath {
 
  inline double fast_pow(double base, double exponent) {
    return vdt::fast_exp(exponent*vdt::fast_log(base));
  }

//   inline float fast_powf(float base, float exponent) {
//     return vdt::fast_expf(exponent*vdt::fast_logf(base));
//   }  
  
}

#endif