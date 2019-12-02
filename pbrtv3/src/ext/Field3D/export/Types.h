//----------------------------------------------------------------------------//

/*
 * Copyright (c) 2009 Sony Pictures Imageworks Inc
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the
 * distribution.  Neither the name of Sony Pictures Imageworks nor the
 * names of its contributors may be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//----------------------------------------------------------------------------//

/*! \file Types.h
  \brief Contains typedefs for the commonly used types in Field3D.
  \ingroup field
*/

//----------------------------------------------------------------------------//

#ifndef _INCLUDED_Field3D_Types_H_
#define _INCLUDED_Field3D_Types_H_

//----------------------------------------------------------------------------//

#include <vector>
#include <limits>

#ifdef FIELD3D_CUSTOM_MATH_LIB
#  include FIELD3D_MATH_LIB_INCLUDE
#else
#  include "StdMathLib.h"
#endif

//----------------------------------------------------------------------------//
// Interval
//----------------------------------------------------------------------------//

//! Represents a single integration interval. 
//! The interval is assumed to be inclusive, i.e. [t0,t1].
struct Interval
{
  // Constructor ---------------------------------------------------------------

  //! Default constructor
  Interval(double start, double end, double step)
    : t0(start), t1(end), stepLength(step) 
  { }

  // Public data members -------------------------------------------------------

  //! The start of the interval (inclusive)
  double t0;
  //! The end of the interval (inclusive)
  double t1;
  //! The world space step length that is reasonable to use for the given 
  //! interval.
  double stepLength;
};

//----------------------------------------------------------------------------//

typedef std::vector<Interval> IntervalVec;

//----------------------------------------------------------------------------//
// Utilities
//----------------------------------------------------------------------------//

template <typename From_T, typename To_T>
To_T clampForType(const From_T v)
{
  // Different behavior for integer vs fp types
  To_T lowestTo;
  From_T lowest;
  if (std::numeric_limits<To_T>::is_integer) {
    lowestTo = std::numeric_limits<To_T>::min();
    lowest   = static_cast<From_T>(lowestTo);
  } else {
    lowestTo  = -std::numeric_limits<To_T>::max();
    lowest    = static_cast<From_T>(lowestTo);
  }
  const To_T   highestTo = std::numeric_limits<To_T>::max();
  const From_T highest   = static_cast<From_T>(highestTo);
  // Perform check
  if (v < lowest) {
    return lowest;
  } else if (v > highest) {
    return highest;
  }
  return v;
}

//----------------------------------------------------------------------------//

#endif // Include guard

