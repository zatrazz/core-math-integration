//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

// MPFR reference implementations for binary{32,64}.
//
// Each function computes the transcendental once at INTERNAL_PRECISION (two
// bits more than the target format) in round-toward-zero, then rounds that
// single high-precision result to the target format for every rounding mode
// selected in MASK.  MPFR's round-toward-zero result plus its exact ternary is
// folded into a sticky bit (round-to-odd) by set_sticky(); rounding a
// round-to-odd value to the narrower target in any mode is double-rounding
// safe, so the per-mode mpfr_get_*() calls yield the same values as computing
// each mode independently (and mpfr_get_*() applies the target exponent range,
// including gradual underflow to subnormals and overflow to infinity).
//
// Each ref_NAME/ref_NAMEf pair is a C++ wrapper (declared in refimpls_mpfr.h)
// around one template instantiation.

#include <bit>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>

#include <mpfr.h>

#include "refimpls_modes.h"
#include "refimpls_mpfr.h"

// Exception side-channel shared with the driver (see refimpls.cc): when
// refimpls_compute_exc is non-zero each evaluation publishes, per rounding
// mode, the FE_* exceptions its rounded result is expected to raise.
extern "C"
{
  extern int refimpls_compute_exc;
  extern __thread unsigned refimpls_last_exc[REF_NRND];
}

namespace
{

static const mpfr_rnd_t ref_rnd_modes[REF_NRND]
    = { MPFR_RNDN, MPFR_RNDU, MPFR_RNDD, MPFR_RNDZ };

// Per-format constants and MPFR conversions.
template <typename F> struct Fmt;
template <> struct Fmt<double>
{
  using U = std::uint64_t;
  static constexpr int prec = DBL_MANT_DIG + 2;
  static constexpr int mant_dig = DBL_MANT_DIG;
  static constexpr U exp_mask = 0x7ffULL;
  static int set (mpfr_ptr y, double x) { return mpfr_set_d (y, x, MPFR_RNDN); }
  static double get (mpfr_srcptr y, mpfr_rnd_t r) { return mpfr_get_d (y, r); }
};
template <> struct Fmt<float>
{
  using U = std::uint32_t;
  static constexpr int prec = FLT_MANT_DIG + 2;
  static constexpr int mant_dig = FLT_MANT_DIG;
  static constexpr U exp_mask = 0xffU;
  static int set (mpfr_ptr y, float x) { return mpfr_set_flt (y, x, MPFR_RNDN); }
  static float get (mpfr_srcptr y, mpfr_rnd_t r) { return mpfr_get_flt (y, r); }
};

// Per-thread scratch reused across calls to avoid a malloc/free of the MPFR
// limbs (and the sticky-bit mpz) on every evaluation.  Initialized lazily on
// first use and intentionally never freed (the process owns it until exit).
template <typename F> struct Scratch
{
  bool inited;
  mpfr_t a, b, c;
  mpz_t sticky;
};

template <typename F>
Scratch<F> &
scratch ()
{
  static thread_local Scratch<F> s;
  if (!s.inited)
    {
      mpfr_init2 (s.a, Fmt<F>::prec);
      mpfr_init2 (s.b, Fmt<F>::prec);
      mpfr_init2 (s.c, Fmt<F>::prec);
      mpz_init (s.sticky);
      s.inited = true;
    }
  return s;
}

// Fold the round-toward-zero ternary INEX into a sticky bit (round-to-odd):
// if the result was inexact, force the least significant bit of the mantissa.
template <typename F>
void
set_sticky (mpfr_ptr r, int inex)
{
  if (inex == 0 || !mpfr_number_p (r))
    return;
  auto &s = scratch<F> ();
  mpfr_exp_t e = mpfr_get_z_2exp (s.sticky, r);
  if (mpz_sgn (s.sticky) < 0)
    {
      mpz_neg (s.sticky, s.sticky);
      mpz_setbit (s.sticky, 0);
      mpz_neg (s.sticky, s.sticky);
    }
  else
    mpz_setbit (s.sticky, 0);
  mpfr_set_z_2exp (r, s.sticky, e, MPFR_RNDN);
}

// Round the high-precision value HI (computed round-toward-zero, INEX its
// ternary) to the target format for every rounding mode selected in MASK.
// When exception computation is enabled, also determine the IEEE exceptions
// each rounded result raises, from the true value HI (still in a wide exponent
// range) versus the rounded result: NaN result -> invalid, infinite result
// from a finite argument -> divide-by-zero (pole), a finite true value that
// does not fit -> overflow, a tiny inexact result -> underflow, and any inexact
// rounding -> inexact.
template <typename F>
void
round_all (mpfr_ptr hi, int inex, unsigned mask, F out[REF_NRND])
{
  set_sticky<F> (hi, inex);

  int base = 0;
  if (refimpls_compute_exc)
    {
      if (mpfr_nan_p (hi))
	base = FE_INVALID;
      else if (mpfr_inf_p (hi))
	base = FE_DIVBYZERO;
    }

  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      {
	F v = Fmt<F>::get (hi, ref_rnd_modes[i]);
	out[i] = v;
	if (!refimpls_compute_exc)
	  continue;

	int e = base;
	if (base == 0) // HI is finite
	  {
	    int inexact = (inex != 0) || (mpfr_cmp_d (hi, v) != 0);
	    // Overflow: the result rounds to infinity, or a directed mode
	    // clamps it to the max finite value while the true value reaches
	    // the first binade above it.  A value merely between the max finite
	    // value and that binade rounds down without overflowing.
	    if (std::isinf (v)
		|| (std::fabs (v) == std::numeric_limits<F>::max () && inexact
		    && mpfr_get_exp (hi)
			   >= std::numeric_limits<F>::max_exponent + 1))
	      e |= FE_OVERFLOW | FE_INEXACT;
	    else if (inexact)
	      {
		e |= FE_INEXACT;
		if (std::fabs (v) < std::numeric_limits<F>::min ()
		    && !mpfr_zero_p (hi))
		  e |= FE_UNDERFLOW;
	      }
	  }
	refimpls_last_exc[i] = (unsigned) e;
      }
}

// Store the mode-independent value V (raising exceptions EXC) into every slot
// selected in MASK.  Used for special cases (NaN/Inf/zero) whose result does
// not depend on rounding.
template <typename F>
void
broadcast (F v, int exc, unsigned mask, F out[REF_NRND])
{
  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      {
	out[i] = v;
	if (refimpls_compute_exc)
	  refimpls_last_exc[i] = (unsigned) exc;
      }
}

template <typename F>
int
arg_is_snan (F x)
{
  using U = typename Fmt<F>::U;
  U a = (U) (std::bit_cast<U> (x) << 1);
  return a > (Fmt<F>::exp_mask << Fmt<F>::mant_dig)
	 && a < ((U) ((Fmt<F>::exp_mask << 1) | 1) << (Fmt<F>::mant_dig - 1));
}

// round_all() for a single-argument function.  round_all() infers the raised
// exceptions from the magnitude of the *result*, which is wrong when the
// *argument* is non-finite: an infinite argument yielding an infinite result
// is not a pole (no divide-by-zero, e.g. exp(inf)=inf), and a quiet-NaN
// argument raises nothing.  Handle those two cases here from the argument, and
// defer every finite argument to round_all() (poles, overflow, underflow and
// finite-domain errors are already correct there).
template <typename F>
void
round_all1 (F x, mpfr_ptr hi, int inex, unsigned mask, F out[REF_NRND])
{
  if (std::isfinite (x))
    {
      round_all<F> (hi, inex, mask, out);
      return;
    }

  int exc = 0;
  if (refimpls_compute_exc)
    exc = std::isnan (x)
	      ? (arg_is_snan<F> (x) ? FE_INVALID : 0)
	      // Infinite argument: signals only when the function is undefined
	      // there and the reference yields NaN, e.g. sin(inf), log(-inf).
	      : (mpfr_nan_p (hi) ? FE_INVALID : 0);

  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      {
	// A NaN argument propagates (glibc returns the input, quieted); an
	// infinite argument takes the reference value (exp(inf)=inf,
	// atan(inf)=pi/2, ...).
	out[i] = std::isnan (x) ? x : Fmt<F>::get (hi, ref_rnd_modes[i]);
	if (refimpls_compute_exc)
	  refimpls_last_exc[i] = (unsigned) exc;
      }
}

// round_all() for a two-operand function whose second floating argument may be
// non-finite (pow, powr) or an integer promoted to floating (pown, rootn,
// compoundn, always finite).  round_all() derives the exceptions from the
// result, which is wrong when an argument is non-finite: the value is always
// the reference's, and the exception follows the argument -- a NaN argument
// propagates and signals only when signaling, while an infinite argument
// signals invalid exactly when the operation is undefined there and the
// reference yields NaN.  With every floating argument finite round_all() is
// already correct.
template <typename F>
void
round_all2 (F x, F y, mpfr_ptr hi, int inex, unsigned mask, F out[REF_NRND])
{
  if (std::isfinite (x) && std::isfinite (y))
    {
      round_all<F> (hi, inex, mask, out);
      return;
    }
  int exc = 0;
  if (refimpls_compute_exc)
    exc = (std::isnan (x) || std::isnan (y))
	      ? ((arg_is_snan<F> (x) || arg_is_snan<F> (y)) ? FE_INVALID : 0)
	      : (mpfr_nan_p (hi) ? FE_INVALID : 0);
  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      {
	out[i] = Fmt<F>::get (hi, ref_rnd_modes[i]);
	if (refimpls_compute_exc)
	  refimpls_last_exc[i] = (unsigned) exc;
      }
}

// Generic bodies, parameterized by the MPFR function.

using MpfrF = int (*) (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
using MpfrFF = int (*) (mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);

template <typename F>
void
ref1 (F x, unsigned mask, F out[REF_NRND], MpfrF fn)
{
  auto &s = scratch<F> ();
  Fmt<F>::set (s.a, x);
  int inex = fn (s.a, s.a, MPFR_RNDZ);
  round_all1<F> (x, s.a, inex, mask, out);
}

template <typename F>
void
ref2 (F x, F y, unsigned mask, F out[REF_NRND], MpfrFF fn)
{
  auto &s = scratch<F> ();
  Fmt<F>::set (s.b, x);
  Fmt<F>::set (s.c, y);
  int inex = fn (s.a, s.b, s.c, MPFR_RNDZ);
  round_all2<F> (x, y, s.a, inex, mask, out);
}

template <typename F, typename IntFn>
void
refi (F x, long long int y, unsigned mask, F out[REF_NRND], IntFn fn)
{
  auto &s = scratch<F> ();
  Fmt<F>::set (s.b, x);
  int inex = fn (s.a, s.b, y, MPFR_RNDZ);
  round_all2<F> (x, (F) y, s.a, inex, mask, out);
}

// Special cases needing bespoke handling.

template <typename F>
void
ref_rsqrt_impl (F x, unsigned mask, F out[REF_NRND])
{
  // mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says rsqrt(-0)
  // should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf.
  if (x == F (0) && F (1) / x < F (0))
    {
      broadcast<F> (F (1) / x, FE_DIVBYZERO, mask, out);
      return;
    }
  auto &s = scratch<F> ();
  Fmt<F>::set (s.a, x);
  int inex = mpfr_rec_sqrt (s.a, s.a, MPFR_RNDZ);
  round_all1<F> (x, s.a, inex, mask, out);
}

template <typename F>
void
ref_atan2_impl (F a0, F a1, unsigned mask, F out[REF_NRND])
{
  // atan2 is finite everywhere except at a NaN argument; a quiet NaN
  // propagates without signalling (round_all would wrongly tag it invalid).
  if (std::isnan (a0) || std::isnan (a1))
    {
      int snan = arg_is_snan<F> (a0) || arg_is_snan<F> (a1);
      broadcast<F> (a0 + a1, snan ? FE_INVALID : 0, mask, out);
      return;
    }
  auto &s = scratch<F> ();
  Fmt<F>::set (s.b, a0);
  Fmt<F>::set (s.c, a1);
  int inex = mpfr_atan2 (s.a, s.b, s.c, MPFR_RNDZ);
  round_all<F> (s.a, inex, mask, out);
}

template <typename F>
void
ref_hypot_impl (F x, F y, unsigned mask, F out[REF_NRND])
{
  // hypot(+-inf, y) = +inf for any y, even a NaN; other non-finite arguments
  // follow NaN propagation.  round_all cannot derive these from the result.
  if (std::isinf (x) || std::isinf (y))
    {
      int snan = arg_is_snan<F> (x) || arg_is_snan<F> (y);
      broadcast<F> ((F) INFINITY, snan ? FE_INVALID : 0, mask, out);
      return;
    }
  if (std::isnan (x) || std::isnan (y))
    {
      int snan = arg_is_snan<F> (x) || arg_is_snan<F> (y);
      broadcast<F> (x + y, snan ? FE_INVALID : 0, mask, out);
      return;
    }
  auto &s = scratch<F> ();
  Fmt<F>::set (s.b, x);
  Fmt<F>::set (s.c, y);
  int inex = mpfr_hypot (s.a, s.b, s.c, MPFR_RNDZ);
  round_all<F> (s.a, inex, mask, out);
}

template <typename F>
void
ref_lgamma_impl (F x, unsigned mask, F out[REF_NRND])
{
  int sign;
  auto &s = scratch<F> ();
  Fmt<F>::set (s.a, x);
  int inex = mpfr_lgamma (s.a, &sign, s.a, MPFR_RNDZ);
  round_all1<F> (x, s.a, inex, mask, out);
}

template <typename F>
void
ref_sincos_impl (F x, unsigned mask, F sinp[REF_NRND], F cosp[REF_NRND])
{
  auto &s = scratch<F> ();
  Fmt<F>::set (s.a, x);
  int inex_sin = mpfr_sin (s.a, s.a, MPFR_RNDZ);
  round_all1<F> (x, s.a, inex_sin, mask, sinp);

  // The single sincos call raises the union of the exceptions of both results;
  // stash sin's before cos overwrites the side-channel, then combine.
  unsigned sinexc[REF_NRND];
  if (refimpls_compute_exc)
    for (int i = 0; i < REF_NRND; i++)
      if (mask & (1u << i))
	sinexc[i] = refimpls_last_exc[i];

  Fmt<F>::set (s.b, x);
  int inex_cos = mpfr_cos (s.b, s.b, MPFR_RNDZ);
  round_all1<F> (x, s.b, inex_cos, mask, cosp);

  if (refimpls_compute_exc)
    for (int i = 0; i < REF_NRND; i++)
      if (mask & (1u << i))
	refimpls_last_exc[i] |= sinexc[i];
}

} // namespace

// extern "C" wrappers: one ref_NAME (double) and ref_NAMEf (float) per function,
// matching the symbols the driver links against.

#define DEF1(cname, mpfrfn)                                                   \
  void ref_##cname (double x, unsigned m, double o[REF_NRND])      \
  {                                                                           \
    ref1<double> (x, m, o, mpfrfn);                                           \
  }                                                                           \
  void ref_##cname##f (float x, unsigned m, float o[REF_NRND])     \
  {                                                                           \
    ref1<float> (x, m, o, mpfrfn);                                            \
  }

#define DEF2(cname, mpfrfn)                                                   \
  void ref_##cname (double x, double y, unsigned m,               \
			       double o[REF_NRND])                            \
  {                                                                           \
    ref2<double> (x, y, m, o, mpfrfn);                                        \
  }                                                                           \
  void ref_##cname##f (float x, float y, unsigned m,              \
				  float o[REF_NRND])                           \
  {                                                                           \
    ref2<float> (x, y, m, o, mpfrfn);                                         \
  }

#define DEFI(cname, mpfrfn)                                                   \
  void ref_##cname (double x, long long int y, unsigned m,        \
			       double o[REF_NRND])                            \
  {                                                                           \
    refi<double> (x, y, m, o, mpfrfn);                                        \
  }                                                                           \
  void ref_##cname##f (float x, long long int y, unsigned m,      \
				  float o[REF_NRND])                           \
  {                                                                           \
    refi<float> (x, y, m, o, mpfrfn);                                         \
  }

DEF1 (acos, mpfr_acos)
DEF1 (acosh, mpfr_acosh)
DEF1 (acospi, mpfr_acospi)
DEF1 (asin, mpfr_asin)
DEF1 (asinh, mpfr_asinh)
DEF1 (asinpi, mpfr_asinpi)
DEF1 (atan, mpfr_atan)
DEF1 (atanh, mpfr_atanh)
DEF1 (atanpi, mpfr_atanpi)
DEF1 (cbrt, mpfr_cbrt)
DEF1 (cos, mpfr_cos)
DEF1 (cosh, mpfr_cosh)
DEF1 (cospi, mpfr_cospi)
DEF1 (erf, mpfr_erf)
DEF1 (erfc, mpfr_erfc)
DEF1 (exp, mpfr_exp)
DEF1 (exp10, mpfr_exp10)
DEF1 (exp10m1, mpfr_exp10m1)
DEF1 (exp2, mpfr_exp2)
DEF1 (exp2m1, mpfr_exp2m1)
DEF1 (expm1, mpfr_expm1)
DEF1 (log, mpfr_log)
DEF1 (log1p, mpfr_log1p)
DEF1 (log2, mpfr_log2)
DEF1 (log2p1, mpfr_log2p1)
DEF1 (log10, mpfr_log10)
DEF1 (log10p1, mpfr_log10p1)
DEF1 (sin, mpfr_sin)
DEF1 (sinh, mpfr_sinh)
DEF1 (sinpi, mpfr_sinpi)
DEF1 (tan, mpfr_tan)
DEF1 (tanh, mpfr_tanh)
DEF1 (tanpi, mpfr_tanpi)
DEF1 (tgamma, mpfr_gamma)

DEF2 (pow, mpfr_pow)
DEF2 (powr, mpfr_powr)

DEFI (compoundn, mpfr_compound_si)
DEFI (pown, mpfr_pown)
DEFI (rootn, mpfr_rootn_si)

void
ref_rsqrt (double x, unsigned m, double o[REF_NRND])
{
  ref_rsqrt_impl<double> (x, m, o);
}
void
ref_rsqrtf (float x, unsigned m, float o[REF_NRND])
{
  ref_rsqrt_impl<float> (x, m, o);
}

void
ref_atan2 (double y, double x, unsigned m, double o[REF_NRND])
{
  ref_atan2_impl<double> (y, x, m, o);
}
void
ref_atan2f (float x, float y, unsigned m, float o[REF_NRND])
{
  ref_atan2_impl<float> (x, y, m, o);
}

void
ref_hypot (double x, double y, unsigned m, double o[REF_NRND])
{
  ref_hypot_impl<double> (x, y, m, o);
}
void
ref_hypotf (float x, float y, unsigned m, float o[REF_NRND])
{
  ref_hypot_impl<float> (x, y, m, o);
}

void
ref_lgamma (double x, unsigned m, double o[REF_NRND])
{
  ref_lgamma_impl<double> (x, m, o);
}
void
ref_lgammaf (float x, unsigned m, float o[REF_NRND])
{
  ref_lgamma_impl<float> (x, m, o);
}

void
ref_sincos (double x, unsigned m, double s[REF_NRND], double c[REF_NRND])
{
  ref_sincos_impl<double> (x, m, s, c);
}
void
ref_sincosf (float x, unsigned m, float s[REF_NRND], float c[REF_NRND])
{
  ref_sincos_impl<float> (x, m, s, c);
}
