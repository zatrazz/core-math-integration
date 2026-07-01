//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef REFIMPLS_MODES_H
#define REFIMPLS_MODES_H

// Rounding mode indices shared between the C reference implementations and the
// C++ driver.  The reference functions compute the multi-precision result once
// and round it to every mode selected in a bitmask (1u << REF_RND*), storing
// each rounded result in out[REF_RND*].  The order matches the C99 rounding
// modes used by the driver (FE_TONEAREST, FE_UPWARD, FE_DOWNWARD,
// FE_TOWARDZERO).

enum ref_round_mode
{
  REF_RNDN,
  REF_RNDU,
  REF_RNDD,
  REF_RNDZ,
  REF_NRND
};

#endif
