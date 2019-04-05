/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

namespace sva
{

/**
 * Compute 1/sinc(x).
 * This code is inspired by boost/math/special_functions/sinc.hpp.
 */
template<typename T>
T sinc_inv(const T x)
{
  const T taylor_0_bound = std::numeric_limits<T>::epsilon();
  const T taylor_2_bound = std::sqrt(taylor_0_bound);
  const T taylor_n_bound = std::sqrt(taylor_2_bound);

  // We use the 4th order taylor series around 0 of x/sin(x) to compute
  // this function:
  //      2      4
  //     x    7⋅x     ⎛ 6⎞
  // 1 + ── + ──── + O⎝x ⎠
  //     6    360
  // this approximation is valid around 0.
  // if x is far from 0, our approximation is not valid
  // since x^6 becomes non negligable we use the normal computation of the function
  // (i.e. taylor_2_bound^6 + taylor_0_bound == taylor_0_bound but
  //       taylor_n_bound^6 + taylor_0_bound != taylor_0).

  if(std::abs(x) >= taylor_n_bound)
  {
    return (x / std::sin(x));
  }
  else
  {
    // x is bellow taylor_n_bound so we don't care of the 6th order term of
    // the taylor series.
    // We set the 0 order term.
    T result = static_cast<T>(1);

    if(std::abs(x) >= taylor_0_bound)
    {
      // x is above the machine epsilon so x^2 is meaningful.
      T x2 = x * x;
      result += x2 / static_cast<T>(6);

      if(std::abs(x) >= taylor_2_bound)
      {
        // x is upper the machine sqrt(epsilon) so x^4 is meaningful.
        result += static_cast<T>(7) * (x2 * x2) / static_cast<T>(360);
      }
    }

    return (result);
  }
}

} // namespace sva
