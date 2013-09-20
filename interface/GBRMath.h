#ifndef GBRMATH
#define GBRMATH

#include "vdt/vdtMath.h"
#include <cmath>

namespace gbrmath {
 
  inline double fast_pow(double base, double exponent) {
    if (base==0. && exponent>0.) return 0.;
    else if (base>0.) return vdt::fast_exp(exponent*vdt::fast_log(base));
    else if (base<0.) return -vdt::fast_exp(exponent*vdt::fast_log(-base));
    else return std::nan("");
  }

  inline float fast_powf(float base, float exponent) {
    if (base==0. && exponent>0.) return 0.;
    else if (base>0.) return vdt::fast_expf(exponent*vdt::fast_logf(base));
    else if (base<0.) return -vdt::fast_expf(exponent*vdt::fast_logf(-base));
    else return std::nanf(""); 
  }
}

#endif