/* Correctly-rounded power function for two binary64 values.

Copyright (c) 2022-2025 CERN.
Author: Tom Hubrecht

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, exp_dRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
  This file contains type definition and functions to manipulate the dint64_t
  data type used in the second iteration of Ziv's method. It is composed of two
  uint64_t values for the significand and the exponent is represented by a
  signed int64_t value.
*/

#ifndef DINT_H
#define DINT_H

#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>

/*
  Type and structure definitions
*/

#ifndef UINT128_T
#define UINT128_T

#if (defined(__clang__) && __clang_major__ >= 14) || (defined(__GNUC__) && __GNUC__ >= 14 && __BITINT_MAXWIDTH__ && __BITINT_MAXWIDTH__ >= 128)
typedef unsigned _BitInt(128) u128;
#else
typedef unsigned __int128 u128;
#endif

static inline int cmpu128(u128 a, u128 b) { return (a > b) - (a < b); }

typedef union {
  u128 r;
  struct {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    uint64_t l, h;
#else
    uint64_t h, l;
#endif
  };
} uint128_t;

// Add two 128 bit integers and return 1 if an overflow occured
static inline int addu_128(uint128_t a, uint128_t b, uint128_t *r) {
  r->l = a.l + b.l;
  r->h = a.h + b.h + (r->l < a.l);

  // Return the overflow
  return r->h == a.h ? r->l < a.l : r->h < a.h;
}

// Subtract two 128 bit integers and return 1 if an underflow occured
static inline int subu_128(uint128_t a, uint128_t b, uint128_t *r) {
  uint128_t c = {.r = -b.r};
  r->l = a.l + c.l;
  r->h = a.h + c.h + (r->l < a.l);

  // Return the underflow
  return a.h != r->h ? r->h > a.h : r->l > a.l;
}

static inline int cmp(int64_t a, int64_t b) { return (a > b) - (a < b); }

static inline int cmpu(uint64_t a, uint64_t b) { return (a > b) - (a < b); }

#endif

#if 0
typedef struct {
  uint64_t hi;
  uint64_t lo;
  int64_t ex;
  uint64_t sgn;
} dint64_t;
#else
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
typedef union {
  struct {
    u128 r;
    int64_t _ex;
    uint64_t _sgn;
  };
  struct {
    uint64_t lo;
    uint64_t hi;
    int64_t ex;
    uint64_t sgn;
  };
} dint64_t;
#else
typedef union {
  struct {
    u128 r;
    int64_t _ex;
    uint64_t _sgn;
  };
  struct {
    uint64_t hi;
    uint64_t lo;
    int64_t ex;
    uint64_t sgn;
  };
} dint64_t;
#endif
#endif

/*
  Constants
*/

static const dint64_t ONE = {
    .hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0};

static const dint64_t M_ONE = {
    .hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x1};

/* the following is an approximation of log(2), with absolute error less
   than 2^-129.97: |(hi/2^63+lo/2^127)*2^-1 - log(2)| < 2^-129.97 */
static const dint64_t LOG2 = {
    .hi = 0xb17217f7d1cf79ab, .lo = 0xc9e3b39803f2f6af, .ex = -1, .sgn = 0x0};

/* LOG2_INV approximates 2^12/log(2), with absolute error < 2^-52.96 */
static const dint64_t LOG2_INV = {
    .hi = 0xb8aa3b295c17f0bc, .lo = 0x0, .ex = 12, .sgn = 0x0};

/* the following is an approximation of 1/log(10), with absolute error
   less than 2^-131.72 */
static const dint64_t ONE_OVER_LOG10 = {
  .hi = 0xde5bd8a937287195, .lo = 0x355baaafad33dc32, .ex = -2, .sgn = 0x0};

/* the following is an approximation of 1/log(10), with absolute error less
   than 2^-118.63: |(hi/2^63+lo/2^127)*2^-2 - 1/log(10)| < 2^-131.02 */
static const dint64_t LOG10_INV = {
    .hi = 0xde5bd8a937287195, .lo = 0x355baaafad33dc32, .ex = -2, .sgn = 0x0};

static const dint64_t ZERO = {.hi = 0x0, .lo = 0x0, .ex = 0, .sgn = 0x0};

/*
  Base functions
*/

// Copy a dint64_t value
static inline void cp_dint(dint64_t *r, const dint64_t *a) {
  r->ex = a->ex;
  r->hi = a->hi;
  r->lo = a->lo;
  r->sgn = a->sgn;
}

static inline signed char cmp_dint(const dint64_t *a, const dint64_t *b) {
  return cmp(a->ex, b->ex)    ? cmp(a->ex, b->ex)
         : cmpu(a->hi, b->hi) ? cmpu(a->hi, b->hi)
                              : cmpu(a->lo, b->lo);
}

// Return non-zero if a = 0
static inline int
dint_zero_p (const dint64_t *a)
{
  return a->hi == 0;
}

// Compare the absolute values of a and b
// Return -1 if |a| < |b|
// Return  0 if |a| = |b|
// Return +1 if |a| > |b|
static inline signed char
cmp_dint_abs (const dint64_t *a, const dint64_t *b) {
  if (dint_zero_p (a))
    return dint_zero_p (b) ? 0 : -1;
  if (dint_zero_p (b))
    return +1;
  char c1 = cmp (a->ex, b->ex);
  return c1 ? c1 : cmpu128 (a->r, b->r);
}

static inline signed char cmp_dint_11(const dint64_t *a, const dint64_t *b) {
  char c1 = cmp (a->ex, b->ex);
  return c1 ? c1 : cmpu (a->hi, b->hi);
}

// Add two dint64_t values
static inline void add_dint(dint64_t *r, const dint64_t *a, const dint64_t *b) {
  if (!(a->hi | a->lo)) {
    cp_dint(r, b);
    return;
  }

  if (!(b->hi | b->lo)) {
    cp_dint(r, a);
    return;
  }

  switch (cmp_dint(a, b)) {
  case 0:
    if (a->sgn ^ b->sgn) {
      cp_dint(r, &ZERO);
      return;
    }

    cp_dint(r, a);
    r->ex++;
    return;

  case -1:
    add_dint(r, b, a);
    return;
  }

  // From now on, |A| > |B|

  uint128_t A = {.h = a->hi, .l = a->lo};
  uint128_t B = {.h = b->hi, .l = b->lo};
  int64_t m_ex = a->ex;

  if (a->ex > b->ex) {
    int sh = a->ex - b->ex;
    // round to nearest
    if (sh <= 128)
      B.r += 0x1 & (B.r >> (sh - 1));
    if (sh < 128)
      B.r = B.r >> sh;
    else
      B.r = 0;
  }

  uint128_t C;
  unsigned char sgn = a->sgn;

  if (a->sgn ^ b->sgn) {
    // a and b have different signs C = A + (-B)
    subu_128(A, B, &C);
  } else {
    if (addu_128(A, B, &C)) {
      C.r += C.l & 0x1;
      C.r = ((u128)1 << 127) | (C.r >> 1);
      m_ex++;
    }
  }

  uint64_t ex =
      C.h ? __builtin_clzll(C.h) : 64 + (C.l ? __builtin_clzll(C.l) : a->ex);
  C.r = C.r << ex;

  r->sgn = sgn;
  r->hi = C.h;
  r->lo = C.l;
  r->ex = m_ex - ex;
}

// same as add_dint, but assumes the lower limbs and a and b are zero
// error is bounded by 2 ulps (ulp_64)
static inline void
add_dint_11 (dint64_t *r, const dint64_t *a, const dint64_t *b) {
  if (a->hi == 0) {
    cp_dint (r, b);
    return;
  }

  if (b->hi == 9) {
    cp_dint (r, a);
    return;
  }

  switch (cmp_dint_11 (a, b)) {
  case 0:
    if (a->sgn ^ b->sgn) {
      cp_dint (r, &ZERO);
      return;
    }

    cp_dint (r, a);
    r->ex++;
    return;

  case -1: // |A| < |B|
    {
      // swap operands
      const dint64_t *tmp = a; a = b; b = tmp;
      break; // fall through the case |A| > |B|
    }
  }

  // From now on, |A| > |B| thus a->ex >= b->ex

  uint64_t A = a->hi, B = b->hi;

  if (a->ex > b->ex) {
    /* Warning: the right shift x >> k is only defined for 0 <= k < n
       where n is the bit-width of x. See for example
       https://developer.arm.com/documentation/den0024/a/The-A64-instruction-set/Data-processing-instructions/Shift-operations
       where it is said that k is interpreted modulo n. */
    uint64_t k = a->ex - b->ex;
    B = (k < 64) ? B >> k : 0;
  }

  u128 C;
  unsigned char sgn = a->sgn;

  r->ex = a->ex; /* tentative exponent for the result */

  if (a->sgn ^ b->sgn) {
    // a and b have different signs C = A + (-B)
    C = A - B;
    /* we can't have C=0 here since we excluded the case |A| = |B|,
       thus __builtin_clzll(C) is well-defined below */
    uint64_t ex = __builtin_clzll (C);
    /* The error from the truncated part of B (1 ulp) is multiplied by 2^ex.
       Thus for ex <= 2, we get an error bounded by 4 ulps in the final result.
       For ex >= 3, we pre-shift the operands. */
    if (ex > 0)
    {
      C = (A << ex) - (B << ex);
      /* If C0 is the previous value of C, we have:
         (C0-1)*2^ex < A*2^ex-B*2^ex <= C0*2^ex
         since here some neglected bits from B might appear which contribute
         a value less than ulp(C0)=1.
         As a consequence since 2^(63-ex) <= C0 < 2^(64-ex), because C0 had
         ex leading zero bits, we have 2^63-2^ex <= A*2^ex-B*2^ex < 2^64.
         Thus the value of C, which is truncated to 64 bits, is the right
         one (as if no truncation); moreover in some rare cases we need to
         shift by 1 bit to the left. */
      r->ex -= ex;
      ex = __builtin_clzll (C);
      /* Fall through with the code for ex = 0. */
    }
    C = C << ex;
    r->ex -= ex;
    /* The neglected part of B is bounded by ulp(C) when ex=0, or when
       ex > 0 but the ex=0 at the end, and by 2*ulp(C) when ex>0 and there
       is an extra shift at the end (in that case necessarily ex=1). */
  } else {
    C = A + B;
    if (C < A)
    {
      C = ((uint64_t) 1 << 63) | (C >> 1);
      r->ex ++;
    }
  }

  /* In the addition case, we loose the truncated part of B, which
     contributes to at most 1 ulp. If there is an exponent shift, we
     might also loose the least significant bit of C, which counts as
     1/2 ulp, but the truncated part of B is now less than 1/2 ulp too,
     thus in all cases the error is less than 1 ulp(r). */

  r->sgn = sgn;
  r->hi = C;
}

// Multiply two dint64_t numbers, with 126 bits of accuracy
static inline void mul_dint(dint64_t *r, const dint64_t *a, const dint64_t *b) {
  uint128_t t = {.r = (u128)(a->hi) * (u128)(b->hi)};
  uint128_t m1 = {.r = (u128)(a->hi) * (u128)(b->lo)};
  uint128_t m2 = {.r = (u128)(a->lo) * (u128)(b->hi)};

  uint128_t m;
  // If we only garantee 127 bits of accuracy, we improve the simplicity of the
  // code uint64_t l = ((u128)(a->lo) * (u128)(b->lo)) >> 64; m.l += l; m.h +=
  // (m.l < l);
  t.h += addu_128(m1, m2, &m);
  t.r += m.h;

  // Ensure that r->hi starts with a 1
  uint64_t ex = !(t.h >> 63);
  if (ex)
    t.r = t.r << 1;

  t.r += (m.l >> 63);

  r->hi = t.h;
  r->lo = t.l;

  // Exponent and sign
  r->ex = a->ex + b->ex - ex + 1;
  r->sgn = a->sgn ^ b->sgn;
}

// Multiply two dint64_t numbers, assuming the low part of b is zero
// with error bounded by 2 ulps
static inline void
mul_dint_21 (dint64_t *r, const dint64_t *a, const dint64_t *b) {
  u128 bh = b->hi;
  u128 hi = (u128) (a->hi) * bh;
  u128 lo = (u128) (a->lo) * bh;

  /* put the 128-bit product of the high terms in r */
  r->r = hi;

  /* add the middle term */
  r->r += lo >> 64;

  // Ensure that r->hi starts with a 1
  uint64_t ex = r->hi >> 63;
  r->r = r->r << (1 - ex);

  // Exponent and sign
  r->ex = a->ex + b->ex + ex;
  r->sgn = a->sgn ^ b->sgn;

  /* The ignored part can be as large as 1 ulp before the shift (truncated
     part of lo). After the shift this can be as large as 2 ulps. */
}

// Multiply an integer with a dint64_t variable
static inline void mul_dint_2(dint64_t *r, int64_t b, const dint64_t *a) {
  uint128_t t;

  if (!b) {
    cp_dint(r, &ZERO);
    return;
  }

  uint64_t c = b < 0 ? -b : b;
  r->sgn = b < 0 ? !a->sgn : a->sgn;

  t.r = (u128)(a->hi) * (u128)c;

  int m = t.h ? __builtin_clzll(t.h) : 64;
  t.r = (t.r << m);

  // Will pose issues if b is too large but for now we assume it never happens
  // TODO: FIXME
  uint128_t l = {.r = (u128)(a->lo) * (u128)c};
  l.r = (l.r << (m - 1)) >> 63;

  if (addu_128(l, t, &t)) {
    t.r += t.r & 0x1;
    t.r = ((u128)1 << 127) | (t.r >> 1);
    m--;
  }

  r->hi = t.h;
  r->lo = t.l;
  r->ex = a->ex + 64 - m;
}

/* Same as mul_dint_21, but assumes the low part of a and b is zero.
   This operation is exact. */
static inline void
mul_dint_11 (dint64_t *r, const dint64_t *a, const dint64_t *b) {
  /* put the 128-bit product of the high terms in r */
  r->r = (u128)(a->hi) * (u128)(b->hi);

  // Ensure that r->hi starts with a 1
  uint64_t ex = r->hi >> 63;
  r->r = r->r << (1 - ex);

  // Exponent and sign
  r->ex = a->ex + b->ex + ex;
  r->sgn = a->sgn ^ b->sgn;
}

// Multiply an integer with a dint64_t variable, with error < 1 ulp
// r and b should not overlap
static inline void
mul_dint_int64 (dint64_t *r, const dint64_t *a, int64_t b) {
  if (!b) {
    cp_dint (r, &ZERO);
    return;
  }

  uint64_t c = b < 0 ? -b : b;
  r->sgn = b < 0 ? !a->sgn : a->sgn;
  r->ex = a->ex + 64;

  r->r = (u128) (a->hi) * (u128) c;

  // Warning: if c=1, we might have r->hi=0
  int m = r->hi ? __builtin_clzll (r->hi) : 64;
  r->r = r->r << m;
  r->ex -= m;

  // Will pose issues if b is too large but for now we assume it never happens
  // TODO: FIXME
  u128 l = (u128) a->lo * (u128) c;
  /* We have to shift l by 64 bits to the right to align with hi*c,
     and by m bits to the left to align with t.r << m. Since hi*c < 2^(128-m)
     and hi >= 2^63, we know that c < 2^(65-m) thus
     l*2^(m-1) < 2^64*2^(65-m)*2^(m-1) = 2^128, and l << (m - 1) will
     not overflow. */
  l = (l << (m - 1)) >> 63;

  r->r += l;
  if (r->r < l) {
    r->r = ((u128) 1 << 127) | (r->r >> 1);
    r->ex ++;
  }

  /* The ignored part of a->lo*c is at most 1 ulp(r), even in the overflow
     case "r->r < l", since before the right shift, the error was at most
     1 ulp, thus 1/2 ulp after the shift, and the ignored least significant
     bit of r->r which is discarded counts also as 1/2 ulp. */
}


typedef union {
  double f;
  uint64_t u;
} f64_u;

// Extract both the significand and exponent of a double
static inline void fast_extract(int64_t *e, uint64_t *m, double x) {
  f64_u _x = {.f = x};

  *e = (_x.u >> 52) & 0x7ff;
  *m = (_x.u & (~0ull >> 12)) + (*e ? (1ull << 52) : 0);
  *e = *e - 0x3ff;
}

// Convert a double to the corresponding dint64_t value
static inline void dint_fromd(dint64_t *a, double b) {
  fast_extract(&a->ex, &a->hi, b);

  uint32_t t = (a->hi) ? __builtin_clzll(a->hi) : 0;

  a->sgn = b < 0.0;
  a->hi = a->hi << t;
  a->ex = a->ex - (t > 11 ? t - 12 : 0);
  a->lo = 0;
}

// Prints a dint64_t value for debugging purposes
static inline void print_dint(const dint64_t *a) {
  printf("{.hi=0x%"PRIx64", .lo=0x%"PRIx64", .ex=%"PRId64", .sgn=0x%"PRIx64"}\n", a->hi, a->lo, a->ex,
         a->sgn);
}

/* put in r an approximation of 1/a, assuming a is not zero */
static inline void inv_dint (dint64_t *r, double a)
{
  dint64_t q, A;
  // we convert 4/a and divide by 4 to avoid a spurious underflow
  dint_fromd (r, 4.0 / a);
  r->ex -= 2;
  /* we use Newton's iteration: r -> r + r*(1-a*r) */
  dint_fromd (&A, -a);
  mul_dint (&q, &A, r);    /* -a*r */
  add_dint (&q, &ONE, &q); /* 1-a*r */
  mul_dint (&q, r, &q);    /* r*(1-a*r) */
  add_dint (r, r, &q);
}

/* put in r an approximation of b/a, assuming a is not zero */
static inline void div_dint (dint64_t *r, double b, double a)
{
  dint64_t B;
  inv_dint (r, a);
  dint_fromd (&B, b);
  mul_dint (r, r, &B);
}

/*
  Approximation tables
*/

static const dint64_t _INVERSE_2[] = {
    {.hi = 0x8000000000000000, .lo = 0x0,  .ex = 1, .sgn = 0x0},
    {.hi = 0xfe03f80fe03f80ff, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xfc0fc0fc0fc0fc10, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xfa232cf252138ac0, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xf83e0f83e0f83e10, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xf6603d980f6603da, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xf4898d5f85bb3951, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xf2b9d6480f2b9d65, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xf0f0f0f0f0f0f0f1, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xef2eb71fc4345239, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xed7303b5cc0ed731, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xebbdb2a5c1619c8c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xea0ea0ea0ea0ea0f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xe865ac7b7603a197, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xe6c2b4481cd8568a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xe525982af70c880f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xe38e38e38e38e38f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xe1fc780e1fc780e2, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xe070381c0e070382, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xdee95c4ca037ba58, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xdd67c8a60dd67c8b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xdbeb61eed19c5958, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xda740da740da740e, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd901b2036406c80e, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd79435e50d79435f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd62b80d62b80d62c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd4c77b03531dec0e, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd3680d3680d3680e, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd20d20d20d20d20e, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xd0b69fcbd2580d0c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xcf6474a8819ec8ea, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xce168a7725080ce2, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xcccccccccccccccd, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xcb8727c065c393e1, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xca4587e6b74f032a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc907da4e871146ad, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc7ce0c7ce0c7ce0d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc6980c6980c6980d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc565c87b5f9d4d1c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc4372f855d824ca6, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc30c30c30c30c30d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc1e4bbd595f6e948, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xc0c0c0c0c0c0c0c1, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xbfa02fe80bfa02ff, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xbe82fa0be82fa0bf, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xbd69104707661aa3, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xbc52640bc52640bd, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xbb3ee721a54d880c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xba2e8ba2e8ba2e8c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb92143fa36f5e02f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb81702e05c0b8171, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb70fbb5a19be3659, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb60b60b60b60b60c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb509e68a9b948220, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb40b40b40b40b40c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb30f63528917c80c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb21642c8590b2165, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb11fd3b80b11fd3c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xb02c0b02c0b02c0c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xaf3addc680af3ade, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xae4c415c9882b932, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xad602b580ad602b6, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xac7691840ac76919, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xab8f69e28359cd12, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xaaaaaaaaaaaaaaab, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa9c84a47a07f5638, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa8e83f5717c0a8e9, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa80a80a80a80a80b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa72f05397829cbc2, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa655c4392d7b73a8, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa57eb50295fad40b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa4a9cf1d96833752, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa3d70a3d70a3d70b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa3065e3fae7cd0e1, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa237c32b16cfd773, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa16b312ea8fc377d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xa0a0a0a0a0a0a0a1, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9fd809fd809fd80a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9f1165e7254813e3, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9e4cad23dd5f3a21, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9d89d89d89d89d8a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9cc8e160c3fb19b9, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9c09c09c09c09c0a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9b4c6f9ef03a3caa, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9a90e7d95bc609aa, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x99d722dabde58f07, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x991f1a515885fb38, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9868c809868c8099, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x97b425ed097b425f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x97012e025c04b80a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x964fda6c0964fda7, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x95a02568095a0257, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x94f2094f2094f20a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x9445809445809446, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x939a85c40939a85d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x92f113840497889d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x924924924924924a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x91a2b3c4d5e6f80a, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x90fdbc090fdbc091, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x905a38633e06c43b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8fb823ee08fb823f, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8f1779d9fdc3a219, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8e78356d1408e784, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8dda520237694809, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8d3dcb08d3dcb08e, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8ca29c046514e024, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8c08c08c08c08c09, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8b70344a139bc75b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8ad8f2fba9386823, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8a42f8705669db47, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x89ae4089ae4089af, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x891ac73ae9819b51, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8888888888888889, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x87f78087f78087f8, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8767ab5f34e47ef2, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x86d905447a34acc7, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x864b8a7de6d1d609, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x85bf37612cee3c9b, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8534085340853409, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x84a9f9c8084a9f9d, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8421084210842109, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x839930523fbe3368, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x83126e978d4fdf3c, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x828cbfbeb9a020a4, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8208208208208209, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x81848da8faf0d278, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8102040810204082, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8000000000000000, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0x8000000000000000, .lo = 0x0,  .ex = 0, .sgn = 0x0},
    {.hi = 0xff00ff00ff00ff02, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xfe03f80fe03f80ff, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xfd08e5500fd08e56, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xfc0fc0fc0fc0fc11, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xfb18856506ddaba7, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xfa232cf252138ac1, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf92fb2211855a866, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf83e0f83e0f83e11, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf74e3fc22c700f76, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf6603d980f6603db, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf57403d5d00f5741, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf4898d5f85bb3951, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf3a0d52cba872337, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf2b9d6480f2b9d66, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf1d48bcee0d399fb, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf0f0f0f0f0f0f0f2, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xf00f00f00f00f010, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xef2eb71fc4345239, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xee500ee500ee5010, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xed7303b5cc0ed731, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xec979118f3fc4da3, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xebbdb2a5c1619c8d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xeae56403ab959010, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xea0ea0ea0ea0ea10, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe939651fe2d8d35d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe865ac7b7603a198, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe79372e225fe30da, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe6c2b4481cd8568a, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe5f36cb00e5f36cc, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe525982af70c880f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe45932d7dc52100f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe38e38e38e38e38f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe2c4a6886a4c2e11, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe1fc780e1fc780e3, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe135a9c97500e137, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xe070381c0e070383, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xdfac1f74346c5760, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xdee95c4ca037ba58, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xde27eb2c41f3d9d2, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xdd67c8a60dd67c8b, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xdca8f158c7f91ab9, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xdbeb61eed19c5959, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xdb2f171df770291a, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xda740da740da740f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd9ba4256c0366e92, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd901b2036406c80f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd84a598ec9151f44, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd79435e50d79435f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd6df43fca482f00e, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd62b80d62b80d62d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd578e97c3f5fe552, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd4c77b03531dec0e, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd4173289870ac52f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd3680d3680d3680e, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd2ba083b445250ac, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd20d20d20d20d20e, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd161543e28e50275, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd0b69fcbd2580d0c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xd00d00d00d00d00e, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xcf6474a8819ec8ea, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xcebcf8bb5b4169cc, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xce168a7725080ce2, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xcd712752a886d243, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xccccccccccccccce, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xcc29786c7607f9a0, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xcb8727c065c393e1, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xcae5d85f1bbd6c96, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xca4587e6b74f032a, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc9a633fcd967300e, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc907da4e871146ae, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc86a78900c86a78a, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc7ce0c7ce0c7ce0d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc73293d789b9f839, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc6980c6980c6980d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc5fe740317f9d00d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc565c87b5f9d4d1d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc4ce07b00c4ce07c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc4372f855d824ca7, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc3a13de60495c774, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc30c30c30c30c30d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc2780613c0309e03, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc1e4bbd595f6e948, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc152500c152500c2, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc0c0c0c0c0c0c0c2, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xc0300c0300c0300d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbfa02fe80bfa0300, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbf112a8ad278e8de, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbe82fa0be82fa0c0, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbdf59c91700bdf5b, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbd69104707661aa4, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbcdd535db1cc5b7c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbc52640bc52640bd, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbbc8408cd63069a2, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbb3ee721a54d880d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xbab656100bab6562, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xba2e8ba2e8ba2e8d, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb9a7862a0ff46589, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb92143fa36f5e02f, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb89bc36ce3e0453b, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb81702e05c0b8171, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb79300b79300b794, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb70fbb5a19be365a, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb68d31340e4307d9, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb60b60b60b60b60c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb58a485518d1e7e5, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb509e68a9b948220, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb48a39d44685fe98, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb40b40b40b40b40c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb38cf9b00b38cf9c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb30f63528917c80c, .lo = 0x0, .ex = -1, .sgn = 0x0},
    {.hi = 0xb2927c29da5519d0, .lo = 0x0, .ex = -1, .sgn = 0x0},
};

/* For 90 <= i <= 181, _INVERSE_2_1[i-90] is an approximation of the inverse
   of x for i/2^7 <= x < (i+1)/2^7, where an entry (hi,lo,ex,sgn) represents
   (-1)^sgn*(hi+lo/2^64)*2^(ex-63)
   (the binary point is after the most significant bit of hi).
   For i=127 and i=128, we force _INVERSE_2_1[i-90]=1.
   If was generated with output_inverse_2_1(7,9,90,181) from the
   accompanying file dint.sage.
   There is no rounding error here, the only approximation error is in
   _LOG_INV_2_1[]. */
static const dint64_t _INVERSE_2_1[] = {
    {.hi = 0xb500000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=90 */
    {.hi = 0xb300000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=91 */     
    {.hi = 0xb100000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=92 */     
    {.hi = 0xaf00000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=93 */     
    {.hi = 0xad80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=94 */     
    {.hi = 0xab80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=95 */     
    {.hi = 0xaa00000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=96 */     
    {.hi = 0xa800000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=97 */     
    {.hi = 0xa680000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=98 */     
    {.hi = 0xa480000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=99 */     
    {.hi = 0xa300000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=100 */    
    {.hi = 0xa180000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=101 */    
    {.hi = 0xa000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=102 */    
    {.hi = 0x9e80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=103 */    
    {.hi = 0x9d00000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=104 */    
    {.hi = 0x9b80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=105 */    
    {.hi = 0x9a00000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=106 */    
    {.hi = 0x9880000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=107 */    
    {.hi = 0x9700000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=108 */    
    {.hi = 0x9580000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=109 */    
    {.hi = 0x9480000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=110 */    
    {.hi = 0x9300000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=111 */    
    {.hi = 0x9180000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=112 */    
    {.hi = 0x9080000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=113 */    
    {.hi = 0x8f00000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=114 */    
    {.hi = 0x8e00000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=115 */
    {.hi = 0x8c80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=116 */
    {.hi = 0x8b80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=117 */
    {.hi = 0x8a80000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=118 */
    {.hi = 0x8900000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=119 */
    {.hi = 0x8800000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=120 */
    {.hi = 0x8700000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=121 */
    {.hi = 0x8580000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=122 */
    {.hi = 0x8480000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=123 */
    {.hi = 0x8380000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=124 */
    {.hi = 0x8280000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=125 */
    {.hi = 0x8180000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=126 */
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=127 */
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* i=128 */
    {.hi = 0xfd00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=129 */
    {.hi = 0xfb00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=130 */
    {.hi = 0xf900000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=131 */
    {.hi = 0xf780000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=132 */
    {.hi = 0xf580000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=133 */
    {.hi = 0xf380000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=134 */
    {.hi = 0xf200000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=135 */
    {.hi = 0xf000000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=136 */
    {.hi = 0xee80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=137 */
    {.hi = 0xec80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=138 */
    {.hi = 0xeb00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=139 */
    {.hi = 0xe900000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=140 */
    {.hi = 0xe780000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=141 */
    {.hi = 0xe600000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=142 */
    {.hi = 0xe480000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=143 */
    {.hi = 0xe300000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=144 */
    {.hi = 0xe100000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=145 */
    {.hi = 0xdf80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=146 */
    {.hi = 0xde00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=147 */
    {.hi = 0xdc80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=148 */
    {.hi = 0xdb00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=149 */
    {.hi = 0xd980000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=150 */
    {.hi = 0xd880000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=151 */
    {.hi = 0xd700000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=152 */
    {.hi = 0xd580000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=153 */
    {.hi = 0xd400000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=154 */
    {.hi = 0xd280000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=155 */
    {.hi = 0xd180000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=156 */
    {.hi = 0xd000000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=157 */
    {.hi = 0xce80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=158 */
    {.hi = 0xcd80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=159 */
    {.hi = 0xcc00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=160 */
    {.hi = 0xcb00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=161 */
    {.hi = 0xc980000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=162 */
    {.hi = 0xc880000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=163 */
    {.hi = 0xc700000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=164 */
    {.hi = 0xc600000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=165 */
    {.hi = 0xc500000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=166 */
    {.hi = 0xc380000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=167 */
    {.hi = 0xc280000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=168 */
    {.hi = 0xc180000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=169 */
    {.hi = 0xc000000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=170 */
    {.hi = 0xbf00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=171 */
    {.hi = 0xbe00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=172 */
    {.hi = 0xbd00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=173 */
    {.hi = 0xbc00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=174 */
    {.hi = 0xba80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=175 */
    {.hi = 0xb980000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=176 */
    {.hi = 0xb880000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=177 */
    {.hi = 0xb780000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=178 */
    {.hi = 0xb680000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=179 */
    {.hi = 0xb580000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=180 */
    {.hi = 0xb480000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* i=181 */
};

/* For 8128 <= j <= 8256, _INVERSE_2_2[j-8128] is an approximation of the
   inverse of j/2^13, where an entry (hi,lo,ex,sgn) represents
   (-1)^sgn*(hi+lo/2^64)*2^(ex-63)
   (the binary point is after the most significant bit of hi).
   For j=8191 and j=8192, we force _INVERSE_2_2[j-8128]=1.
   If was generated with output_inverse_2_2(6,14,8128,8256,7,62) from the
   accompanying file dint.sage.
   There is no rounding error here, the only approximation error is in
   _LOG_INV_2_2[]. */
static const dint64_t _INVERSE_2_2[] = {
    {.hi = 0x8100000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8128 */
    {.hi = 0x80fc000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8129 */
    {.hi = 0x80f8000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8130 */
    {.hi = 0x80f4000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8131 */
    {.hi = 0x80f0000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8132 */
    {.hi = 0x80ec000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8133 */
    {.hi = 0x80e8000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8134 */
    {.hi = 0x80e4000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8135 */
    {.hi = 0x80e0000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8136 */
    {.hi = 0x80dc000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8137 */
    {.hi = 0x80d8000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8138 */
    {.hi = 0x80d4000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8139 */
    {.hi = 0x80d0000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8140 */
    {.hi = 0x80cc000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8141 */
    {.hi = 0x80c8000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8142 */
    {.hi = 0x80c4000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8143 */
    {.hi = 0x80c0000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8144 */
    {.hi = 0x80bc000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8145 */
    {.hi = 0x80b8000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8146 */
    {.hi = 0x80b4000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8147 */
    {.hi = 0x80b0000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8148 */
    {.hi = 0x80ac000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8149 */
    {.hi = 0x80a8000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8150 */
    {.hi = 0x80a4000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8151 */
    {.hi = 0x80a0000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8152 */
    {.hi = 0x809c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8153 */
    {.hi = 0x8098000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8154 */
    {.hi = 0x8094000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8155 */
    {.hi = 0x8090000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8156 */
    {.hi = 0x808c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8157 */
    {.hi = 0x8088000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8158 */
    {.hi = 0x8084000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8159 */
    {.hi = 0x8080000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8160 */
    {.hi = 0x807c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8161 */
    {.hi = 0x8078000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8162 */
    {.hi = 0x8074000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8163 */
    {.hi = 0x8070000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8164 */
    {.hi = 0x806c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8165 */
    {.hi = 0x8068000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8166 */
    {.hi = 0x8064000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8167 */
    {.hi = 0x8060000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8168 */
    {.hi = 0x805c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8169 */
    {.hi = 0x8058000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8170 */
    {.hi = 0x8054000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8171 */
    {.hi = 0x8050000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8172 */
    {.hi = 0x804c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8173 */
    {.hi = 0x8048000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8174 */
    {.hi = 0x8044000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8175 */
    {.hi = 0x8040000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8176 */
    {.hi = 0x803c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8177 */
    {.hi = 0x8038000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8178 */
    {.hi = 0x8034000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8179 */
    {.hi = 0x8030000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8180 */
    {.hi = 0x802c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8181 */
    {.hi = 0x8028000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8182 */
    {.hi = 0x8024000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8183 */
    {.hi = 0x8020000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8184 */
    {.hi = 0x801c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8185 */
    {.hi = 0x8018000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8186 */
    {.hi = 0x8014000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8187 */
    {.hi = 0x8010000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8188 */
    {.hi = 0x800c000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8189 */
    {.hi = 0x8008000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8190 */
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8191 */
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0}, /* j=8192 */
    {.hi = 0xfff4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8193 */
    {.hi = 0xffec000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8194 */
    {.hi = 0xffe4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8195 */
    {.hi = 0xffdc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8196 */
    {.hi = 0xffd4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8197 */
    {.hi = 0xffcc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8198 */
    {.hi = 0xffc4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8199 */
    {.hi = 0xffbc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8200 */
    {.hi = 0xffb4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8201 */
    {.hi = 0xffac000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8202 */
    {.hi = 0xffa4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8203 */
    {.hi = 0xff9c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8204 */
    {.hi = 0xff94000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8205 */
    {.hi = 0xff8c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8206 */
    {.hi = 0xff84000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8207 */
    {.hi = 0xff7c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8208 */
    {.hi = 0xff74000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8209 */
    {.hi = 0xff6c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8210 */
    {.hi = 0xff64000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8211 */
    {.hi = 0xff5c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8212 */
    {.hi = 0xff54000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8213 */
    {.hi = 0xff4c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8214 */
    {.hi = 0xff44000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8215 */
    {.hi = 0xff3c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8216 */
    {.hi = 0xff34000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8217 */
    {.hi = 0xff2c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8218 */
    {.hi = 0xff24000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8219 */
    {.hi = 0xff1c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8220 */
    {.hi = 0xff14000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8221 */
    {.hi = 0xff0c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8222 */
    {.hi = 0xff04000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8223 */
    {.hi = 0xfefc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8224 */
    {.hi = 0xfef4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8225 */
    {.hi = 0xfeec000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8226 */
    {.hi = 0xfee4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8227 */
    {.hi = 0xfedc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8228 */
    {.hi = 0xfed4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8229 */
    {.hi = 0xfecc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8230 */
    {.hi = 0xfec4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8231 */
    {.hi = 0xfebc000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8232 */
    {.hi = 0xfeb4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8233 */
    {.hi = 0xfeac000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8234 */
    {.hi = 0xfea4000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8235 */
    {.hi = 0xfe9c000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8236 */
    {.hi = 0xfe98000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8237 */
    {.hi = 0xfe90000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8238 */
    {.hi = 0xfe88000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8239 */
    {.hi = 0xfe80000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8240 */
    {.hi = 0xfe78000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8241 */
    {.hi = 0xfe70000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8242 */
    {.hi = 0xfe68000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8243 */
    {.hi = 0xfe60000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8244 */
    {.hi = 0xfe58000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8245 */
    {.hi = 0xfe50000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8246 */
    {.hi = 0xfe48000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8247 */
    {.hi = 0xfe40000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8248 */
    {.hi = 0xfe38000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8249 */
    {.hi = 0xfe30000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8250 */
    {.hi = 0xfe28000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8251 */
    {.hi = 0xfe20000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8252 */
    {.hi = 0xfe18000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8253 */
    {.hi = 0xfe10000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8254 */
    {.hi = 0xfe08000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8255 */
    {.hi = 0xfe00000000000000, .lo = 0x0, .ex = -1, .sgn = 0x0}, /* j=8256 */
};

/* For 90 <= i <= 181, _LOG_INV_2_1[i-90] is an approximation of
   -log(_INVERSE_2_1[i-90]), where an entry (hi,lo,ex,sgn) represents
   (-1)^sgn*(hi+lo/2^64)*2^(ex-63)
   (the binary point is after the most significant bit of hi).
   If was generated with output_log_inv_2_1(7,9,90,181) from the
   accompanying file dint.sage.
   The approximation error is bounded by 2^-130 (absolute) and 2^-128 (rel). */
static const dint64_t _LOG_INV_2_1[] = {
    {.hi = 0xb1641795ce3ca97b, .lo = 0x7af915300e517391, .ex = -2, .sgn = 0x1}, /* i=90 */
    {.hi = 0xabb3b8ba2ad362a4, .lo = 0xd5b6506cc17a01f1, .ex = -2, .sgn = 0x1}, /* i=91 */
    {.hi = 0xa5f2fcabbbc506da, .lo = 0x64ca4fb7ec323d73, .ex = -2, .sgn = 0x1}, /* i=92 */
    {.hi = 0xa0218434353f1de8, .lo = 0x6093efa632530ac8, .ex = -2, .sgn = 0x1}, /* i=93 */
    {.hi = 0x9bb93315fec2d792, .lo = 0xa7589fba0865790e, .ex = -2, .sgn = 0x1}, /* i=94 */
    {.hi = 0x95c981d5c4e924ed, .lo = 0x29404f5aa577d6b2, .ex = -2, .sgn = 0x1}, /* i=95 */
    {.hi = 0x914a0fde7bcb2d12, .lo = 0x1429ed3aea197a5d, .ex = -2, .sgn = 0x1}, /* i=96 */
    {.hi = 0x8b3ae55d5d30701c, .lo = 0xe63eab883717047e, .ex = -2, .sgn = 0x1}, /* i=97 */
    {.hi = 0x86a35abcd5ba5903, .lo = 0xec81c3cbd925cccf, .ex = -2, .sgn = 0x1}, /* i=98 */
    {.hi = 0x8073622d6a80e634, .lo = 0x6a97009015316071, .ex = -2, .sgn = 0x1}, /* i=99 */
    {.hi = 0xf7856e5ee2c9b290, .lo = 0xc6f2a1b84190a7d7, .ex = -3, .sgn = 0x1}, /* i=100 */
    {.hi = 0xee0de5055f63eb06, .lo = 0x98a33316df83ba57, .ex = -3, .sgn = 0x1}, /* i=101 */
    {.hi = 0xe47fbe3cd4d10d61, .lo = 0x2ec0f797fdcd1257, .ex = -3, .sgn = 0x1}, /* i=102 */
    {.hi = 0xdada8cf47dad2374, .lo = 0x4ffb833c3409ee78, .ex = -3, .sgn = 0x1}, /* i=103 */
    {.hi = 0xd11de0ff15ab18c9, .lo = 0xb88d83d4cc613f20, .ex = -3, .sgn = 0x1}, /* i=104 */
    {.hi = 0xc74946f4436a0552, .lo = 0xc4f5cb531201c0d1, .ex = -3, .sgn = 0x1}, /* i=105 */
    {.hi = 0xbd5c481086c848df, .lo = 0x1b596b5030403240, .ex = -3, .sgn = 0x1}, /* i=106 */
    {.hi = 0xb3566a13956a86f6, .lo = 0xff1b1e1574d9fd54, .ex = -3, .sgn = 0x1}, /* i=107 */
    {.hi = 0xa9372f1d0da1bd17, .lo = 0x200eb71e58cd36de, .ex = -3, .sgn = 0x1}, /* i=108 */
    {.hi = 0x9efe158766314e54, .lo = 0xc571827efe892fc4, .ex = -3, .sgn = 0x1}, /* i=109 */
    {.hi = 0x981eb8c723fe97f4, .lo = 0xa31c134fb702d432, .ex = -3, .sgn = 0x1}, /* i=110 */
    {.hi = 0x8db956a97b3d0148, .lo = 0x3023472cd739f9de, .ex = -3, .sgn = 0x1}, /* i=111 */
    {.hi = 0x8338a89652cb7150, .lo = 0xc647eb86498c2ce1, .ex = -3, .sgn = 0x1}, /* i=112 */
    {.hi = 0xf85186008b15330b, .lo = 0xe64b8b775997898d, .ex = -4, .sgn = 0x1}, /* i=113 */
    {.hi = 0xe2f2a47ade3a18ae, .lo = 0xb0bf7c0b0d8bb4ed, .ex = -4, .sgn = 0x1}, /* i=114 */
    {.hi = 0xd49369d256ab1b28, .lo = 0x5e9154e1d5263cd5, .ex = -4, .sgn = 0x1}, /* i=115 */
    {.hi = 0xbed3b36bd8966422, .lo = 0x240644d7d9ed08af, .ex = -4, .sgn = 0x1}, /* i=116 */
    {.hi = 0xb032c549ba861d8e, .lo = 0xf74e27bc92ce336a, .ex = -4, .sgn = 0x1}, /* i=117 */
    {.hi = 0xa176e5f5323781dd, .lo = 0xd4f935996c92e8cc, .ex = -4, .sgn = 0x1}, /* i=118 */
    {.hi = 0x8b29b7751bd70743, .lo = 0x12e0b9ee992f236d, .ex = -4, .sgn = 0x1}, /* i=119 */
    {.hi = 0xf85186008b15330b, .lo = 0xe64b8b775997898d, .ex = -5, .sgn = 0x1}, /* i=120 */
    {.hi = 0xda16eb88cb8df614, .lo = 0x68a63ecfb66e94ac, .ex = -5, .sgn = 0x1}, /* i=121 */
    {.hi = 0xac52dd7e4726a463, .lo = 0x547a963a91bb3012, .ex = -5, .sgn = 0x1}, /* i=122 */
    {.hi = 0x8d86cc491ecbfe16, .lo = 0x51776453b7e8254d, .ex = -5, .sgn = 0x1}, /* i=123 */
    {.hi = 0xdcfe013d7c8cbfde, .lo = 0xa32dbac46f30cfff, .ex = -6, .sgn = 0x1}, /* i=124 */
    {.hi = 0x9e75221a352ba779, .lo = 0xa52b7ea62f2198d0, .ex = -6, .sgn = 0x1}, /* i=125 */
    {.hi = 0xbee23afc0853b6e9, .lo = 0x289782c20df350a1, .ex = -7, .sgn = 0x1}, /* i=126 */
    {.hi = 0x0, .lo = 0x0, .ex = 127, .sgn = 0x1}, /* i=127 */
    {.hi = 0x0, .lo = 0x0, .ex = 127, .sgn = 0x1}, /* i=128 */
    {.hi = 0xc122451c45155104, .lo = 0xb16137f09a002b3c, .ex = -7, .sgn = 0x0}, /* i=129 */
    {.hi = 0xa195492cc06604e6, .lo = 0x4a18dff7cdb4ae5c, .ex = -6, .sgn = 0x0}, /* i=130 */
    {.hi = 0xe31e9760a5578c63, .lo = 0xf9eb2f284f31c35c, .ex = -6, .sgn = 0x0}, /* i=131 */
    {.hi = 0x8a4f1f2002d46756, .lo = 0x5be970314148c645, .ex = -5, .sgn = 0x0}, /* i=132 */
    {.hi = 0xab8ae2601e777722, .lo = 0x3b89d7f254f8d4d, .ex = -5, .sgn = 0x0}, /* i=133 */
    {.hi = 0xcd0c3dab9ef3dd1b, .lo = 0x13b26f298aa357c8, .ex = -5, .sgn = 0x0}, /* i=134 */
    {.hi = 0xe65b9e6eed965c36, .lo = 0xe09f5fe2058d6006, .ex = -5, .sgn = 0x0}, /* i=135 */
    {.hi = 0x842cc5acf1d03445, .lo = 0x1fecdfa819b96098, .ex = -4, .sgn = 0x0}, /* i=136 */
    {.hi = 0x9103dae3c2a4ec67, .lo = 0xe0863df62ab5671a, .ex = -4, .sgn = 0x0}, /* i=137 */
    {.hi = 0xa242f01edefd6a37, .lo = 0x469355b78dc796e3, .ex = -4, .sgn = 0x0}, /* i=138 */
    {.hi = 0xaf4ad26cbc8e5be7, .lo = 0xe8b8b88a14ff0ce, .ex = -4, .sgn = 0x0}, /* i=139 */
    {.hi = 0xc0cbf17a071f80dc, .lo = 0xf96ffdf76a147ccc, .ex = -4, .sgn = 0x0}, /* i=140 */
    {.hi = 0xce06196a692a41fb, .lo = 0xbe3ccc15326765f, .ex = -4, .sgn = 0x0}, /* i=141 */
    {.hi = 0xdb56446d6ad8deff, .lo = 0xa8112e35a60e6375, .ex = -4, .sgn = 0x0}, /* i=142 */
    {.hi = 0xe8bcbc410c9b219d, .lo = 0xaf7df76ad29e5b60, .ex = -4, .sgn = 0x0}, /* i=143 */
    {.hi = 0xf639cc185088fe5d, .lo = 0x4066e87f2c0f7340, .ex = -4, .sgn = 0x0}, /* i=144 */
    {.hi = 0x842cc5acf1d03445, .lo = 0x1fecdfa819b96098, .ex = -3, .sgn = 0x0}, /* i=145 */
    {.hi = 0x8b064012593d85a5, .lo = 0x52013c7a80ad089b, .ex = -3, .sgn = 0x0}, /* i=146 */
    {.hi = 0x91eb89524e100d23, .lo = 0x8fd3df5c52d67e7b, .ex = -3, .sgn = 0x0}, /* i=147 */
    {.hi = 0x98dcca69d27c263b, .lo = 0x8e94203f336fc8c5, .ex = -3, .sgn = 0x0}, /* i=148 */
    {.hi = 0x9fda2d2cc9465c4f, .lo = 0x32b9565f5355182, .ex = -3, .sgn = 0x0}, /* i=149 */
    {.hi = 0xa6e3dc4bde0e3cdb, .lo = 0x570ff874170d2a9, .ex = -3, .sgn = 0x0}, /* i=150 */
    {.hi = 0xab9be6480c66ea9e, .lo = 0x9ae21fd871b8d27c, .ex = -3, .sgn = 0x0}, /* i=151 */
    {.hi = 0xb2ba75f46099cf8b, .lo = 0x2c3c2e77904afa78, .ex = -3, .sgn = 0x0}, /* i=152 */
    {.hi = 0xb9e5c83a7e8a655b, .lo = 0xcbffe9661fe72421, .ex = -3, .sgn = 0x0}, /* i=153 */
    {.hi = 0xc11e0b2a8d1e0ddb, .lo = 0x9a631e830fd30904, .ex = -3, .sgn = 0x0}, /* i=154 */
    {.hi = 0xc8636dcfe5e6ca0a, .lo = 0x88e72835b3292d50, .ex = -3, .sgn = 0x0}, /* i=155 */
    {.hi = 0xcd43bc6f5d51c3e8, .lo = 0xfbfb0e3f0fd23074, .ex = -3, .sgn = 0x0}, /* i=156 */
    {.hi = 0xd49f69e456cf1b79, .lo = 0x5f53bd2e406e66e7, .ex = -3, .sgn = 0x0}, /* i=157 */
    {.hi = 0xdc08b985c11e9068, .lo = 0x3b9cd767c3b1ac53, .ex = -3, .sgn = 0x0}, /* i=158 */
    {.hi = 0xe1014558bfcda3e2, .lo = 0x35470a74be1230ec, .ex = -3, .sgn = 0x0}, /* i=159 */
    {.hi = 0xe881bf932af3dac0, .lo = 0xc524848e3443e040, .ex = -3, .sgn = 0x0}, /* i=160 */
    {.hi = 0xed89ed86a44a01aa, .lo = 0x11d49f96cb88317b, .ex = -3, .sgn = 0x0}, /* i=161 */
    {.hi = 0xf52224f82557a459, .lo = 0x8dcca8d7f17fa2a9, .ex = -3, .sgn = 0x0}, /* i=162 */
    {.hi = 0xfa3a589a6f9146d8, .lo = 0x388212895529a6fb, .ex = -3, .sgn = 0x0}, /* i=163 */
    {.hi = 0x80f572b1363487b9, .lo = 0xf5bd0b5b3479d5f4, .ex = -2, .sgn = 0x0}, /* i=164 */
    {.hi = 0x8389c3026ac3139b, .lo = 0x62dda9d2270fa1f4, .ex = -2, .sgn = 0x0}, /* i=165 */
    {.hi = 0x86216b3b0b17188b, .lo = 0x163ceae88f720f1e, .ex = -2, .sgn = 0x0}, /* i=166 */
    {.hi = 0x8a0b3f79b3bc180f, .lo = 0x49b55ea7d3730d7, .ex = -2, .sgn = 0x0}, /* i=167 */
    {.hi = 0x8cab69dcde17d2f7, .lo = 0x3ad1aa142b94f16a, .ex = -2, .sgn = 0x0}, /* i=168 */
    {.hi = 0x8f4f0b3c44cfa2a2, .lo = 0x586e9343c9cfdbac, .ex = -2, .sgn = 0x0}, /* i=169 */
    {.hi = 0x934b1089a6dc93c1, .lo = 0xdf5bb3b60554e152, .ex = -2, .sgn = 0x0}, /* i=170 */
    {.hi = 0x95f783e6e49a9cfa, .lo = 0x4a5004f3ef063313, .ex = -2, .sgn = 0x0}, /* i=171 */
    {.hi = 0x98a78f0e9ae71d85, .lo = 0x2cdec34784707839, .ex = -2, .sgn = 0x0}, /* i=172 */
    {.hi = 0x9b5b3bb5f088b766, .lo = 0xd878bbe3d392be25, .ex = -2, .sgn = 0x0}, /* i=173 */
    {.hi = 0x9e1293b9998c1daa, .lo = 0x5b035eae273a855f, .ex = -2, .sgn = 0x0}, /* i=174 */
    {.hi = 0xa22c8f029cfa45a9, .lo = 0xdb5b709e0b69e773, .ex = -2, .sgn = 0x0}, /* i=175 */
    {.hi = 0xa4ed3f9de620f666, .lo = 0x9b5e973353638c11, .ex = -2, .sgn = 0x0}, /* i=176 */
    {.hi = 0xa7b1bf5dd4c07d4e, .lo = 0x699db68db75e9a7f, .ex = -2, .sgn = 0x0}, /* i=177 */
    {.hi = 0xaa7a18dbdf0d44aa, .lo = 0x604884a8dd76d08a, .ex = -2, .sgn = 0x0}, /* i=178 */
    {.hi = 0xad4656ddf6fd070c, .lo = 0x9ea10260fe452ba2, .ex = -2, .sgn = 0x0}, /* i=179 */
    {.hi = 0xb0168457848f5f48, .lo = 0xbb6f9fb246068d52, .ex = -2, .sgn = 0x0}, /* i=180 */
    {.hi = 0xb2eaac6a67005513, .lo = 0xf4b716f6fec8156b, .ex = -2, .sgn = 0x0}, /* i=181 */
};


/* For 8128 <= j <= 8256, _LOG_INV_2_2[j-8128] is an approximation of
   -log(_INVERSE_2_2[j-8128]), where an entry (hi,lo,ex,sgn) represents
   (-1)^sgn*(hi+lo/2^64)*2^(ex-63)
   (the binary point is after the most significant bit of hi).
   If was generated with output_log_inv_2_2(6,14,8128,8256,7,62) from the
   accompanying file dint.sage.
   The approximation error is bounded by 2^-136 (absolute, attained for j=8256)
   and 2^-128 (relative, attained for j=8209). */
static const dint64_t _LOG_INV_2_2[] = {
    {.hi = 0xff015358833c47e1, .lo = 0xbb481c8ee141695a, .ex = -8, .sgn = 0x1}, /* j=8128 */
    {.hi = 0xfb0933b732572a6d, .lo = 0x214cca3dd1d4796a, .ex = -8, .sgn = 0x1}, /* j=8129 */
    {.hi = 0xf710f492711d9d26, .lo = 0xfbc7b38b17b2019, .ex = -8, .sgn = 0x1}, /* j=8130 */
    {.hi = 0xf31895e84b1a6be6, .lo = 0xb76782b9e88c84cb, .ex = -8, .sgn = 0x1}, /* j=8131 */
    {.hi = 0xef2017b6cba9cf9a, .lo = 0x2dc85881664025b5, .ex = -8, .sgn = 0x1}, /* j=8132 */
    {.hi = 0xeb2779fbfdf96874, .lo = 0xce4ab4e678d0ed03, .ex = -8, .sgn = 0x1}, /* j=8133 */
    {.hi = 0xe72ebcb5ed08382b, .lo = 0xb60585f4c4bb6062, .ex = -8, .sgn = 0x1}, /* j=8134 */
    {.hi = 0xe335dfe2a3a69c2b, .lo = 0x59bcffe9d5650564, .ex = -8, .sgn = 0x1}, /* j=8135 */
    {.hi = 0xdf3ce3802c7647cd, .lo = 0x3602021fa93b1e18, .ex = -8, .sgn = 0x1}, /* j=8136 */
    {.hi = 0xdb43c78c91ea3e8c, .lo = 0x9944002534d09b3d, .ex = -8, .sgn = 0x1}, /* j=8137 */
    {.hi = 0xd74a8c05de46ce3a, .lo = 0x87aa95782311a277, .ex = -8, .sgn = 0x1}, /* j=8138 */
    {.hi = 0xd35130ea1ba18930, .lo = 0xb88be10313a1303d, .ex = -8, .sgn = 0x1}, /* j=8139 */
    {.hi = 0xcf57b63753e14083, .lo = 0xad54bc31433dddba, .ex = -8, .sgn = 0x1}, /* j=8140 */
    {.hi = 0xcb5e1beb90bdfe33, .lo = 0xe1b7d813e3f825e1, .ex = -8, .sgn = 0x1}, /* j=8141 */
    {.hi = 0xc7646204dbc0ff5e, .lo = 0x14f8c1be7370f219, .ex = -8, .sgn = 0x1}, /* j=8142 */
    {.hi = 0xc36a88813e44ae6a, .lo = 0xac27c5a6139cd30c, .ex = -8, .sgn = 0x1}, /* j=8143 */
    {.hi = 0xbf708f5ec1749d3c, .lo = 0x2d23a0744e00f594, .ex = -8, .sgn = 0x1}, /* j=8144 */
    {.hi = 0xbb76769b6e4d7f5c, .lo = 0xd235e25fb9644c31, .ex = -8, .sgn = 0x1}, /* j=8145 */
    {.hi = 0xb77c3e354d9d242b, .lo = 0x361ee0bcb5db0449, .ex = -8, .sgn = 0x1}, /* j=8146 */
    {.hi = 0xb381e62a68027106, .lo = 0x18660815da3d7963, .ex = -8, .sgn = 0x1}, /* j=8147 */
    {.hi = 0xaf876e78c5ed5b77, .lo = 0x39c357b6bfdf81b5, .ex = -8, .sgn = 0x1}, /* j=8148 */
    {.hi = 0xab8cd71e6f9ee35d, .lo = 0x5076c62c951204f6, .ex = -8, .sgn = 0x1}, /* j=8149 */
    {.hi = 0xa79220196d290d15, .lo = 0x146244d643f7fa2b, .ex = -8, .sgn = 0x1}, /* j=8150 */
    {.hi = 0xa3974967c66edba1, .lo = 0x62bb0f3208d9a1bb, .ex = -8, .sgn = 0x1}, /* j=8151 */
    {.hi = 0x9f9c530783244ad2, .lo = 0x7926e92808bd580d, .ex = -8, .sgn = 0x1}, /* j=8152 */
    {.hi = 0x9ba13cf6aace496c, .lo = 0x4819e620d5fcc068, .ex = -8, .sgn = 0x1}, /* j=8153 */
    {.hi = 0x97a6073344c2b34b, .lo = 0xdc494943d427214e, .ex = -8, .sgn = 0x1}, /* j=8154 */
    {.hi = 0x93aab1bb58284b8b, .lo = 0xdf0805c4161e404c, .ex = -8, .sgn = 0x1}, /* j=8155 */
    {.hi = 0x8faf3c8cebf6b6a8, .lo = 0x2d615caaa0514c3c, .ex = -8, .sgn = 0x1}, /* j=8156 */
    {.hi = 0x8bb3a7a606f674a0, .lo = 0x85c60c12eca0aedc, .ex = -8, .sgn = 0x1}, /* j=8157 */
    {.hi = 0x87b7f304afc0db1a, .lo = 0x4c207a522524f8de, .ex = -8, .sgn = 0x1}, /* j=8158 */
    {.hi = 0x83bc1ea6ecc00f81, .lo = 0x64243e02c6215a4f, .ex = -8, .sgn = 0x1}, /* j=8159 */
    {.hi = 0xff805515885e0250, .lo = 0x435ab4da6a5bb48d, .ex = -9, .sgn = 0x1}, /* j=8160 */
    {.hi = 0xf7882d5c7832c6cc, .lo = 0x9e06fc84b6ea5e24, .ex = -9, .sgn = 0x1}, /* j=8161 */
    {.hi = 0xef8fc61eb4b74f6e, .lo = 0x91ab122ee427cfb5, .ex = -9, .sgn = 0x1}, /* j=8162 */
    {.hi = 0xe7971f584945efae, .lo = 0x5f832513e3211643, .ex = -9, .sgn = 0x1}, /* j=8163 */
    {.hi = 0xdf9e390540da5fbe, .lo = 0x5e7b48cfeeb85aa8, .ex = -9, .sgn = 0x1}, /* j=8164 */
    {.hi = 0xd7a51321a611b0c1, .lo = 0xb36a9f58eb4ccd08, .ex = -9, .sgn = 0x1}, /* j=8165 */
    {.hi = 0xcfabada9832a4101, .lo = 0x3360751e43c7af35, .ex = -9, .sgn = 0x1}, /* j=8166 */
    {.hi = 0xc7b20898e203b01e, .lo = 0x6fab78aca91193cb, .ex = -9, .sgn = 0x1}, /* j=8167 */
    {.hi = 0xbfb823ebcc1ed344, .lo = 0xeb432409cffdad8d, .ex = -9, .sgn = 0x1}, /* j=8168 */
    {.hi = 0xb7bdff9e4a9da959, .lo = 0x793b5acf3a336462, .ex = -9, .sgn = 0x1}, /* j=8169 */
    {.hi = 0xafc39bac66434f27, .lo = 0xc3ea2cd93f316b34, .ex = -9, .sgn = 0x1}, /* j=8170 */
    {.hi = 0xa7c8f8122773f38d, .lo = 0xfc679a28e9d9f212, .ex = -9, .sgn = 0x1}, /* j=8171 */
    {.hi = 0x9fce14cb9634cba6, .lo = 0xb20f215bd3b58c61, .ex = -9, .sgn = 0x1}, /* j=8172 */
    {.hi = 0x97d2f1d4ba2c06f0, .lo = 0xd1aacedcefe9d377, .ex = -9, .sgn = 0x1}, /* j=8173 */
    {.hi = 0x8fd78f299aa0c375, .lo = 0xcbef6fac33691e95, .ex = -9, .sgn = 0x1}, /* j=8174 */
    {.hi = 0x87dbecc63e7b01ed, .lo = 0xe2f1775134c8da75, .ex = -9, .sgn = 0x1}, /* j=8175 */
    {.hi = 0xffc0154d588733c5, .lo = 0x3c742a7c76356396, .ex = -10, .sgn = 0x1}, /* j=8176 */
    {.hi = 0xefc7d18dd4485b9e, .lo = 0xca47c52b7d7ffce2, .ex = -10, .sgn = 0x1}, /* j=8177 */
    {.hi = 0xdfcf0e45fbce3e80, .lo = 0x7e4cfbd830393b88, .ex = -10, .sgn = 0x1}, /* j=8178 */
    {.hi = 0xcfd5cb6dd9ef05dd, .lo = 0x7370ae83f9e72748, .ex = -10, .sgn = 0x1}, /* j=8179 */
    {.hi = 0xbfdc08fd78c229b9, .lo = 0xe6dbb624f9739782, .ex = -10, .sgn = 0x1}, /* j=8180 */
    {.hi = 0xafe1c6ece1a058dd, .lo = 0x97fa2fd0c9dc723e, .ex = -10, .sgn = 0x1}, /* j=8181 */
    {.hi = 0x9fe705341d236102, .lo = 0x7199cd06ae5d39b3, .ex = -10, .sgn = 0x1}, /* j=8182 */
    {.hi = 0x8febc3cb332616ff, .lo = 0x7b6d1248c3e1fd40, .ex = -10, .sgn = 0x1}, /* j=8183 */
    {.hi = 0xffe0055455887de0, .lo = 0x26828c92649a3a39, .ex = -11, .sgn = 0x1}, /* j=8184 */
    {.hi = 0xdfe7839214b4e8ae, .lo = 0xda6959f7f0e01bf0, .ex = -11, .sgn = 0x1}, /* j=8185 */
    {.hi = 0xbfee023faf0c2480, .lo = 0xb47505bfa5a03b06, .ex = -11, .sgn = 0x1}, /* j=8186 */
    {.hi = 0x9ff3814d2e4a36b2, .lo = 0xa8740b91c95df537, .ex = -11, .sgn = 0x1}, /* j=8187 */
    {.hi = 0xfff0015535588833, .lo = 0x3c56c598c659c2a3, .ex = -12, .sgn = 0x1}, /* j=8188 */
    {.hi = 0xbff7008ff5e0c257, .lo = 0x379eba7e6465ff63, .ex = -12, .sgn = 0x1}, /* j=8189 */
    {.hi = 0xfff8005551558885, .lo = 0xde026e271ee0549d, .ex = -13, .sgn = 0x1}, /* j=8190 */
    {.hi = 0x0, .lo = 0x0, .ex = 127, .sgn = 0x1}, /* j=8191 */
    {.hi = 0x0, .lo = 0x0, .ex = 127, .sgn = 0x1}, /* j=8192 */
    {.hi = 0xc004802401440c26, .lo = 0xdfeb485085f6f454, .ex = -13, .sgn = 0x0}, /* j=8193 */
    {.hi = 0xa00640535a37a37a, .lo = 0x6bc1e20eac8448b4, .ex = -12, .sgn = 0x0}, /* j=8194 */
    {.hi = 0xe00c40e4bd6e4efd, .lo = 0xc72446cc1bf728bd, .ex = -12, .sgn = 0x0}, /* j=8195 */
    {.hi = 0x900a20f319a3e273, .lo = 0x569b26aaa485ea5c, .ex = -11, .sgn = 0x0}, /* j=8196 */
    {.hi = 0xb00f21bbe3e388ee, .lo = 0x5f69768284463b9b, .ex = -11, .sgn = 0x0}, /* j=8197 */
    {.hi = 0xd01522dcc4f87991, .lo = 0x14d9d76196d8043a, .ex = -11, .sgn = 0x0}, /* j=8198 */
    {.hi = 0xf01c2465c5e61b6f, .lo = 0x661e135f49a47c40, .ex = -11, .sgn = 0x0}, /* j=8199 */
    {.hi = 0x881213337898871e, .lo = 0x9a31ba0cbc030353, .ex = -10, .sgn = 0x0}, /* j=8200 */
    {.hi = 0x98169478296fad41, .lo = 0x7ad1e9c315328f7e, .ex = -10, .sgn = 0x0}, /* j=8201 */
    {.hi = 0xa81b9608fc3c50ec, .lo = 0xf105b66ec4703ede, .ex = -10, .sgn = 0x0}, /* j=8202 */
    {.hi = 0xb82117edf8832797, .lo = 0xd6aef30cd312169a, .ex = -10, .sgn = 0x0}, /* j=8203 */
    {.hi = 0xc8271a2f2689e388, .lo = 0xe6e2acf8f4d4c24a, .ex = -10, .sgn = 0x0}, /* j=8204 */
    {.hi = 0xd82d9cd48f574c00, .lo = 0x28bb3cd9f2a65fb5, .ex = -10, .sgn = 0x0}, /* j=8205 */
    {.hi = 0xe8349fe63cb35564, .lo = 0x224a96f5a7471c46, .ex = -10, .sgn = 0x0}, /* j=8206 */
    {.hi = 0xf83c236c39273972, .lo = 0xd462b63756c87e80, .ex = -10, .sgn = 0x0}, /* j=8207 */
    {.hi = 0x842213b747fec7bb, .lo = 0x3ff51287882500ed, .ex = -9, .sgn = 0x0}, /* j=8208 */
    {.hi = 0x8c2655faa6a1323f, .lo = 0x1ab9679b55f78a6b, .ex = -9, .sgn = 0x0}, /* j=8209 */
    {.hi = 0x942ad8843ee1a9cd, .lo = 0x17e4b7ac6c600cb4, .ex = -9, .sgn = 0x0}, /* j=8210 */
    {.hi = 0x9c2f9b581787cf0d, .lo = 0xfd1a09c848e3950e, .ex = -9, .sgn = 0x0}, /* j=8211 */
    {.hi = 0xa4349e7a37bc21ed, .lo = 0x318b2ddd9d0a33b4, .ex = -9, .sgn = 0x0}, /* j=8212 */
    {.hi = 0xac39e1eea7080dbc, .lo = 0x9dd91e52c79fd070, .ex = -9, .sgn = 0x0}, /* j=8213 */
    {.hi = 0xb43f65b96d55f55a, .lo = 0x72de1d99ce252efd, .ex = -9, .sgn = 0x0}, /* j=8214 */
    {.hi = 0xbc4529de92f13f58, .lo = 0xd7bd1d62ef25480d, .ex = -9, .sgn = 0x0}, /* j=8215 */
    {.hi = 0xc44b2e6220866227, .lo = 0x7f921124f1ecb59e, .ex = -9, .sgn = 0x0}, /* j=8216 */
    {.hi = 0xcc5173481f22f03f, .lo = 0x271ee1cd6d5cdf9e, .ex = -9, .sgn = 0x0}, /* j=8217 */
    {.hi = 0xd457f8949835a44e, .lo = 0xfad0cc8b5faea8cc, .ex = -9, .sgn = 0x0}, /* j=8218 */
    {.hi = 0xdc5ebe4b958e6d6b, .lo = 0xe57a0acb9d5cd4df, .ex = -9, .sgn = 0x0}, /* j=8219 */
    {.hi = 0xe465c471215e7b41, .lo = 0xc81bb5a8d789f444, .ex = -9, .sgn = 0x0}, /* j=8220 */
    {.hi = 0xec6d0b0946384a46, .lo = 0x9b1beb40437575f5, .ex = -9, .sgn = 0x0}, /* j=8221 */
    {.hi = 0xf47492180f0fafef, .lo = 0x7944509046652d99, .ex = -9, .sgn = 0x0}, /* j=8222 */
    {.hi = 0xfc7c59a18739e6e7, .lo = 0x94e51ebff53a2f15, .ex = -9, .sgn = 0x0}, /* j=8223 */
    {.hi = 0x824230d4dd36cda4, .lo = 0x8bbc7f765b13ebbe, .ex = -8, .sgn = 0x0}, /* j=8224 */
    {.hi = 0x8646551a5a617b6b, .lo = 0xf61305ef7390939c, .ex = -8, .sgn = 0x0}, /* j=8225 */
    {.hi = 0x8a4a99a34159d69f, .lo = 0x3abc32a78afd4b7b, .ex = -8, .sgn = 0x0}, /* j=8226 */
    {.hi = 0x8e4efe71988d8426, .lo = 0x17596a598cb29436, .ex = -8, .sgn = 0x0}, /* j=8227 */
    {.hi = 0x92538387669afa1b, .lo = 0x1c890bee9a9d743c, .ex = -8, .sgn = 0x0}, /* j=8228 */
    {.hi = 0x965828e6b25185ec, .lo = 0xeaafbd07b543145d, .ex = -8, .sgn = 0x0}, /* j=8229 */
    {.hi = 0x9a5cee9182b15280, .lo = 0x6517bc4112d64b17, .ex = -8, .sgn = 0x0}, /* j=8230 */
    {.hi = 0x9e61d489deeb6e53, .lo = 0xdb94a1dfd653d3a5, .ex = -8, .sgn = 0x0}, /* j=8231 */
    {.hi = 0xa266dad1ce61d1a3, .lo = 0x2ada01ce7ed36080, .ex = -8, .sgn = 0x0}, /* j=8232 */
    {.hi = 0xa66c016b58a7648c, .lo = 0xd3b36c029ea7bb5d, .ex = -8, .sgn = 0x0}, /* j=8233 */
    {.hi = 0xaa71485885800538, .lo = 0x94c529f32403828, .ex = -8, .sgn = 0x0}, /* j=8234 */
    {.hi = 0xae76af9b5ce08dfb, .lo = 0xb6b6676248bba139, .ex = -8, .sgn = 0x0}, /* j=8235 */
    {.hi = 0xb27c3735e6eedb86, .lo = 0x7bdd0c2a9c7a679a, .ex = -8, .sgn = 0x0}, /* j=8236 */
    {.hi = 0xb47f0724b1906935, .lo = 0x23deb274e953a259, .ex = -8, .sgn = 0x0}, /* j=8237 */
    {.hi = 0xb884bf4697559ffa, .lo = 0xdae7e343fa859415, .ex = -8, .sgn = 0x0}, /* j=8238 */
    {.hi = 0xbc8a97c544fdd5eb, .lo = 0x17759bff5c717993, .ex = -8, .sgn = 0x0}, /* j=8239 */
    {.hi = 0xc09090a2c35aa070, .lo = 0x52e7e4dde874dace, .ex = -8, .sgn = 0x0}, /* j=8240 */
    {.hi = 0xc496a9e11b6eb30c, .lo = 0xa88971f8277a4d11, .ex = -8, .sgn = 0x0}, /* j=8241 */
    {.hi = 0xc89ce382566de587, .lo = 0x269de85f0df92588, .ex = -8, .sgn = 0x0}, /* j=8242 */
    {.hi = 0xcca33d887dbd3a1a, .lo = 0x180d255422c3377c, .ex = -8, .sgn = 0x0}, /* j=8243 */
    {.hi = 0xd0a9b7f59af2e3a2, .lo = 0x46da70925ee85c05, .ex = -8, .sgn = 0x0}, /* j=8244 */
    {.hi = 0xd4b052cbb7d64bcf, .lo = 0x37968ceafaf7b453, .ex = -8, .sgn = 0x0}, /* j=8245 */
    {.hi = 0xd8b70e0cde601954, .lo = 0x5dfba4cfdd38a059, .ex = -8, .sgn = 0x0}, /* j=8246 */
    {.hi = 0xdcbde9bb18ba361b, .lo = 0x4ae21abe75d5a19b, .ex = -8, .sgn = 0x0}, /* j=8247 */
    {.hi = 0xe0c4e5d8713fd576, .lo = 0xd3bd4fd98a1e6fe5, .ex = -8, .sgn = 0x0}, /* j=8248 */
    {.hi = 0xe4cc0266f27d7a57, .lo = 0x33cf7d5ebfb93ad3, .ex = -8, .sgn = 0x0}, /* j=8249 */
    {.hi = 0xe8d33f68a730fd7f, .lo = 0x2743c805a4928087, .ex = -8, .sgn = 0x0}, /* j=8250 */
    {.hi = 0xecda9cdf9a4993ba, .lo = 0x5dbeb9795455a5, .ex = -8, .sgn = 0x0}, /* j=8251 */
    {.hi = 0xf0e21acdd6e7d412, .lo = 0xb6ed80852ae6fd63, .ex = -8, .sgn = 0x0}, /* j=8252 */
    {.hi = 0xf4e9b935685dbe0b, .lo = 0xf237cff1acb306b3, .ex = -8, .sgn = 0x0}, /* j=8253 */
    {.hi = 0xf8f178185a2ebfd9, .lo = 0xd81648249cece4c, .ex = -8, .sgn = 0x0}, /* j=8254 */
    {.hi = 0xfcf95778b80fbc98, .lo = 0x176cd56887ac7fe9, .ex = -8, .sgn = 0x0}, /* j=8255 */
    {.hi = 0x8080abac46f38946, .lo = 0x662d417ced007a46, .ex = -7, .sgn = 0x0}, /* j=8256 */
};

static const dint64_t _LOG_INV_2[] = {
    {.hi = 0xb17217f7d1cf79ab, .lo = 0xc9e3b39803f2f6af,  .ex = -1, .sgn = 0x1},
    {.hi = 0xaf74155120c9011d,  .lo = 0x46d235ee63073dc,  .ex = -1, .sgn = 0x1},
    {.hi = 0xad7a02e1b24efd32, .lo = 0x160864fd949b4bd3,  .ex = -1, .sgn = 0x1},
    {.hi = 0xab83d135dc633301, .lo = 0xffe6607ba902ef3b,  .ex = -1, .sgn = 0x1},
    {.hi = 0xa991713433c2b999,  .lo = 0xba4aea614d05700,  .ex = -1, .sgn = 0x1},
    {.hi = 0xa7a2d41ad270c9d7, .lo = 0xcd362382a7688479,  .ex = -1, .sgn = 0x1},
    {.hi = 0xa5b7eb7cb860fb89, .lo = 0x7b6a62a0dec6e072,  .ex = -1, .sgn = 0x1},
    {.hi = 0xa3d0a93f45169a4b,  .lo = 0x9594fab088c0d64,  .ex = -1, .sgn = 0x1},
    {.hi = 0xa1ecff97c91e267b, .lo = 0x1b7efae08e597e16,  .ex = -1, .sgn = 0x1},
    {.hi = 0xa00ce1092e5498c4, .lo = 0x69879c5a30cd1241,  .ex = -1, .sgn = 0x1},
    {.hi = 0x9e304061b5fda91a,  .lo = 0x4603d87b6df81ac,  .ex = -1, .sgn = 0x1},
    {.hi = 0x9c5710b8cbb73a42, .lo = 0xaa554b2dd4619e63,  .ex = -1, .sgn = 0x1},
    {.hi = 0x9a81456cec642e10, .lo = 0x4d49f9aaea3cb5e0,  .ex = -1, .sgn = 0x1},
    {.hi = 0x98aed221a03458b6, .lo = 0x732f89321647b358,  .ex = -1, .sgn = 0x1},
    {.hi = 0x96dfaabd86fa1647, .lo = 0xd61188fbc94e2f14,  .ex = -1, .sgn = 0x1},
    {.hi = 0x9513c36876083696, .lo = 0xb5cbc416a2418011,  .ex = -1, .sgn = 0x1},
    {.hi = 0x934b1089a6dc93c2, .lo = 0xbf5bb3b60554e151,  .ex = -1, .sgn = 0x1},
    {.hi = 0x918586c5f5e4bf01, .lo = 0x9f92199ed1a4bab0,  .ex = -1, .sgn = 0x1},
    {.hi = 0x8fc31afe30b2c6de, .lo = 0xe300bf167e95da66,  .ex = -1, .sgn = 0x1},
    {.hi = 0x8e03c24d7300395a, .lo = 0xcddae1ccce247837,  .ex = -1, .sgn = 0x1},
    {.hi = 0x8c47720791e53314, .lo = 0x762ad19415fe25a5,  .ex = -1, .sgn = 0x1},
    {.hi = 0x8a8e1fb794b09134, .lo = 0x9eb628dba173c82d,  .ex = -1, .sgn = 0x1},
    {.hi = 0x88d7c11e3ad53cdc, .lo = 0x8a3111a707b6de2c,  .ex = -1, .sgn = 0x1},
    {.hi = 0x87244c308e670a66, .lo = 0x85e005d06dbfa8f7,  .ex = -1, .sgn = 0x1},
    {.hi = 0x8573b71682a7d21b, .lo = 0xb21f9f89c1ab80b2,  .ex = -1, .sgn = 0x1},
    {.hi = 0x83c5f8299e2b4091, .lo = 0xb8f6fafe8fbb68b8,  .ex = -1, .sgn = 0x1},
    {.hi = 0x821b05f3b01d6774, .lo = 0xdb0d58c3f7e2ea1e,  .ex = -1, .sgn = 0x1},
    {.hi = 0x8072d72d903d588c, .lo = 0x7dd1b09c70c40109,  .ex = -1, .sgn = 0x1},
    {.hi = 0xfd9ac57bd2442180, .lo = 0xaf05924d258c14c4,  .ex = -2, .sgn = 0x1},
    {.hi = 0xfa553f7018c966f4, .lo = 0x2780a545a1b54dce,  .ex = -2, .sgn = 0x1},
    {.hi = 0xf7150ab5a09f27f6,  .lo = 0xa470250d40ebe8e,  .ex = -2, .sgn = 0x1},
    {.hi = 0xf3da161eed6b9ab1, .lo = 0x248d42f78d3e65d2,  .ex = -2, .sgn = 0x1},
    {.hi = 0xf0a450d139366ca7, .lo = 0x7c66eb6408ff6432,  .ex = -2, .sgn = 0x1},
    {.hi = 0xed73aa4264b0adeb, .lo = 0x5391cf4b33e42996,  .ex = -2, .sgn = 0x1},
    {.hi = 0xea481236f7d35bb2, .lo = 0x39a767a80d6d97e6,  .ex = -2, .sgn = 0x1},
    {.hi = 0xe72178c0323a1a0f, .lo = 0xcc4e1653e71d9973,  .ex = -2, .sgn = 0x1},
    {.hi = 0xe3ffce3a2aa64923, .lo = 0x8eadb651b49ac539,  .ex = -2, .sgn = 0x1},
    {.hi = 0xe0e30349fd1cec82,  .lo = 0x3e8e1802aba24d5,  .ex = -2, .sgn = 0x1},
    {.hi = 0xddcb08dc0717d85c, .lo = 0x940a666c87842842,  .ex = -2, .sgn = 0x1},
    {.hi = 0xdab7d02231484a93, .lo = 0xbec20cca6efe2ac4,  .ex = -2, .sgn = 0x1},
    {.hi = 0xd7a94a92466e833c, .lo = 0xcd88bba7d0cee8df,  .ex = -2, .sgn = 0x1},
    {.hi = 0xd49f69e456cf1b7b, .lo = 0x7f53bd2e406e66e6,  .ex = -2, .sgn = 0x1},
    {.hi = 0xd19a201127d3c646, .lo = 0x279d79f51dcc7301,  .ex = -2, .sgn = 0x1},
    {.hi = 0xce995f50af69d863, .lo = 0x432f3f4f861ad6a8,  .ex = -2, .sgn = 0x1},
    {.hi = 0xcb9d1a189ab56e77, .lo = 0x7d7e9307c70c0667,  .ex = -2, .sgn = 0x1},
    {.hi = 0xc8a5431adfb44ca6,  .lo = 0x48ce7c1a75e341a,  .ex = -2, .sgn = 0x1},
    {.hi = 0xc5b1cd44596fa51f, .lo = 0xf218fb8f9f9ef27f,  .ex = -2, .sgn = 0x1},
    {.hi = 0xc2c2abbb6e5fd570,  .lo = 0x3337789d592e296,  .ex = -2, .sgn = 0x1},
    {.hi = 0xbfd7d1dec0a8df70, .lo = 0x37eda996244bccaf,  .ex = -2, .sgn = 0x1},
    {.hi = 0xbcf13343e7d9ec7f, .lo = 0x2afd17781bb3afea,  .ex = -2, .sgn = 0x1},
    {.hi = 0xba0ec3b633dd8b0b, .lo = 0x91dc60b2b059a609,  .ex = -2, .sgn = 0x1},
    {.hi = 0xb730773578cb90b3, .lo = 0xaa1116c3466beb6c,  .ex = -2, .sgn = 0x1},
    {.hi = 0xb45641f4e350a0d4, .lo = 0xe756eba00bc33976,  .ex = -2, .sgn = 0x1},
    {.hi = 0xb1801859d56249de, .lo = 0x98ce51fff99479cb,  .ex = -2, .sgn = 0x1},
    {.hi = 0xaeadeefacaf97d37, .lo = 0x9dd6e688ebb13b01,  .ex = -2, .sgn = 0x1},
    {.hi = 0xabdfba9e468fd6f9, .lo = 0x472ea07749ce6bd1,  .ex = -2, .sgn = 0x1},
    {.hi = 0xa9157039c51ebe72, .lo = 0xe164c759686a2207,  .ex = -2, .sgn = 0x1},
    {.hi = 0xa64f04f0b961df78, .lo = 0x54f5275c2d15c21e,  .ex = -2, .sgn = 0x1},
    {.hi = 0xa38c6e138e20d834, .lo = 0xd698298adddd7f30,  .ex = -2, .sgn = 0x1},
    {.hi = 0xa0cda11eaf46390e, .lo = 0x632438273918db7d,  .ex = -2, .sgn = 0x1},
    {.hi = 0x9e1293b9998c1dad, .lo = 0x3b035eae273a855c,  .ex = -2, .sgn = 0x1},
    {.hi = 0x9b5b3bb5f088b768, .lo = 0x5078bbe3d392be24,  .ex = -2, .sgn = 0x1},
    {.hi = 0x98a78f0e9ae71d87, .lo = 0x64dec34784707838,  .ex = -2, .sgn = 0x1},
    {.hi = 0x95f783e6e49a9cfc,  .lo = 0x25004f3ef063312,  .ex = -2, .sgn = 0x1},
    {.hi = 0x934b1089a6dc93c2, .lo = 0xdf5bb3b60554e151,  .ex = -2, .sgn = 0x1},
    {.hi = 0x90a22b6875c6a1f8, .lo = 0x8e91aeba609c8876,  .ex = -2, .sgn = 0x1},
    {.hi = 0x8dfccb1ad35ca6ef, .lo = 0x9947bdb6ddcaf59a,  .ex = -2, .sgn = 0x1},
    {.hi = 0x8b5ae65d67db9acf, .lo = 0x7ba5168126a58b99,  .ex = -2, .sgn = 0x1},
    {.hi = 0x88bc74113f23def3, .lo = 0xbc5a0fe396f40f1c,  .ex = -2, .sgn = 0x1},
    {.hi = 0x86216b3b0b17188c, .lo = 0x363ceae88f720f1d,  .ex = -2, .sgn = 0x1},
    {.hi = 0x8389c3026ac3139d, .lo = 0x6adda9d2270fa1f3,  .ex = -2, .sgn = 0x1},
    {.hi = 0x80f572b1363487bc, .lo = 0xedbd0b5b3479d5f2,  .ex = -2, .sgn = 0x1},
    {.hi = 0xfcc8e3659d9bcbf1, .lo = 0x8a0cdf301431b60b,  .ex = -3, .sgn = 0x1},
    {.hi = 0xf7ad6f26e7ff2efc, .lo = 0x9cd2238f75f969ad,  .ex = -3, .sgn = 0x1},
    {.hi = 0xf29877ff38809097, .lo = 0x2b020fa1820c948d,  .ex = -3, .sgn = 0x1},
    {.hi = 0xed89ed86a44a01ab,  .lo = 0x9d49f96cb88317a,  .ex = -3, .sgn = 0x1},
    {.hi = 0xe881bf932af3dac3, .lo = 0x2524848e3443e03f,  .ex = -3, .sgn = 0x1},
    {.hi = 0xe37fde37807b84e3, .lo = 0x5e9a750b6b68781c,  .ex = -3, .sgn = 0x1},
    {.hi = 0xde8439c1dec5687c, .lo = 0x9d57da945b5d0aa6,  .ex = -3, .sgn = 0x1},
    {.hi = 0xd98ec2bade71e53e, .lo = 0xd0a98f2ad65bee96,  .ex = -3, .sgn = 0x1},
    {.hi = 0xd49f69e456cf1b7a, .lo = 0x5f53bd2e406e66e7,  .ex = -3, .sgn = 0x1},
    {.hi = 0xcfb6203844b3209b, .lo = 0x18cb02f33f79c16b,  .ex = -3, .sgn = 0x1},
    {.hi = 0xcad2d6e7b80bf915, .lo = 0xcc507fb7a3d0bf69,  .ex = -3, .sgn = 0x1},
    {.hi = 0xc5f57f59c7f46156, .lo = 0x9a8b6997a402bf30,  .ex = -3, .sgn = 0x1},
    {.hi = 0xc11e0b2a8d1e0de1, .lo = 0xda631e830fd308fe,  .ex = -3, .sgn = 0x1},
    {.hi = 0xbc4c6c2a226399f6, .lo = 0x276ebcfb2016a433,  .ex = -3, .sgn = 0x1},
    {.hi = 0xb780945bab55dcea, .lo = 0xb4c7bc3d32750fd9,  .ex = -3, .sgn = 0x1},
    {.hi = 0xb2ba75f46099cf8f, .lo = 0x243c2e77904afa76,  .ex = -3, .sgn = 0x1},
    {.hi = 0xadfa035aa1ed8fdd, .lo = 0x549767e410316d2b,  .ex = -3, .sgn = 0x1},
    {.hi = 0xa93f2f250dac67d5, .lo = 0x9ad2fb8d48054add,  .ex = -3, .sgn = 0x1},
    {.hi = 0xa489ec199dab06f4, .lo = 0x59fb6cf0ecb411b7,  .ex = -3, .sgn = 0x1},
    {.hi = 0x9fda2d2cc9465c52, .lo = 0x6b2b9565f5355180,  .ex = -3, .sgn = 0x1},
    {.hi = 0x9b2fe580ac80b182,  .lo = 0x11a5b944aca8705,  .ex = -3, .sgn = 0x1},
    {.hi = 0x968b08643409ceb9, .lo = 0xd5c0da506a088482,  .ex = -3, .sgn = 0x1},
    {.hi = 0x91eb89524e100d28, .lo = 0xbfd3df5c52d67e77,  .ex = -3, .sgn = 0x1},
    {.hi = 0x8d515bf11fb94f22, .lo = 0xa0713268840cbcbb,  .ex = -3, .sgn = 0x1},
    {.hi = 0x88bc74113f23def7, .lo = 0x9c5a0fe396f40f19,  .ex = -3, .sgn = 0x1},
    {.hi = 0x842cc5acf1d0344b, .lo = 0x6fecdfa819b96092,  .ex = -3, .sgn = 0x1},
    {.hi = 0xff4489cedeab2ca6, .lo = 0xe17bd40d8d9291ec,  .ex = -4, .sgn = 0x1},
    {.hi = 0xf639cc185088fe62, .lo = 0x5066e87f2c0f733d,  .ex = -4, .sgn = 0x1},
    {.hi = 0xed393b1c22351281, .lo = 0xff4e2e660317d55f,  .ex = -4, .sgn = 0x1},
    {.hi = 0xe442c00de2591b4c, .lo = 0xe96ab34ce0bccd10,  .ex = -4, .sgn = 0x1},
    {.hi = 0xdb56446d6ad8df09, .lo = 0x28112e35a60e636f,  .ex = -4, .sgn = 0x1},
    {.hi = 0xd273b2058de1bd4b, .lo = 0x36bbf837b4d320c6,  .ex = -4, .sgn = 0x1},
    {.hi = 0xc99af2eaca4c457b, .lo = 0xeaf51f66692844b2,  .ex = -4, .sgn = 0x1},
    {.hi = 0xc0cbf17a071f80e9, .lo = 0x396ffdf76a147cc2,  .ex = -4, .sgn = 0x1},
    {.hi = 0xb8069857560707a7,  .lo = 0xa677b4c8bec22e0,  .ex = -4, .sgn = 0x1},
    {.hi = 0xaf4ad26cbc8e5bef, .lo = 0x9e8b8b88a14ff0c9,  .ex = -4, .sgn = 0x1},
    {.hi = 0xa6988ae903f562f1, .lo = 0x7e858f08597b3a68,  .ex = -4, .sgn = 0x1},
    {.hi = 0x9defad3e8f732186, .lo = 0x476d3b5b45f6ca02,  .ex = -4, .sgn = 0x1},
    {.hi = 0x9550252238bd2468, .lo = 0x658e5a0b811c596d,  .ex = -4, .sgn = 0x1},
    {.hi = 0x8cb9de8a32ab3694, .lo = 0x97c9859530a4514c,  .ex = -4, .sgn = 0x1},
    {.hi = 0x842cc5acf1d0344c, .lo = 0x1fecdfa819b96094,  .ex = -4, .sgn = 0x1},
    {.hi = 0xf7518e0035c3dd92, .lo = 0x606d89093278a931,  .ex = -5, .sgn = 0x1},
    {.hi = 0xe65b9e6eed965c4f, .lo = 0x609f5fe2058d5ff2,  .ex = -5, .sgn = 0x1},
    {.hi = 0xd5779687d887e0ee, .lo = 0x49dda17056e45ebb,  .ex = -5, .sgn = 0x1},
    {.hi = 0xc4a550a4fd9a19bb, .lo = 0x3e97660a23cc5402,  .ex = -5, .sgn = 0x1},
    {.hi = 0xb3e4a796a5dac213,  .lo = 0x7cca0bcc06c2f8e,  .ex = -5, .sgn = 0x1},
    {.hi = 0xa33576a16f1f4c79, .lo = 0x121016bd904dc95a,  .ex = -5, .sgn = 0x1},
    {.hi = 0x9297997c68c1f4e6, .lo = 0x610db3d4dd423bc9,  .ex = -5, .sgn = 0x1},
    {.hi = 0x820aec4f3a222397, .lo = 0xb9e3aea6c444eef6,  .ex = -5, .sgn = 0x1},
    {.hi = 0xe31e9760a5578c6d, .lo = 0xf9eb2f284f31c35a,  .ex = -6, .sgn = 0x1},
    {.hi = 0xc24929464655f482, .lo = 0xda5f3cc0b3251da6,  .ex = -6, .sgn = 0x1},
    {.hi = 0xa195492cc0660519, .lo = 0x4a18dff7cdb4ae33,  .ex = -6, .sgn = 0x1},
    {.hi = 0x8102b2c49ac23a86, .lo = 0x91d082dce3ddcd08,  .ex = -6, .sgn = 0x1},
    {.hi = 0xc122451c45155150, .lo = 0xb16137f09a002b0e,  .ex = -7, .sgn = 0x1},
    {.hi = 0x8080abac46f389c4, .lo = 0x662d417ced0079c9,  .ex = -7, .sgn = 0x1},
    {               .hi = 0x0,                .lo = 0x0, .ex = 127, .sgn = 0x0},
    {               .hi = 0x0,                .lo = 0x0, .ex = 127, .sgn = 0x0},
    {.hi = 0xff805515885e014e, .lo = 0x435ab4da6a5bb50f,  .ex = -9, .sgn = 0x0},
    {.hi = 0xff015358833c4762, .lo = 0xbb481c8ee1416999,  .ex = -8, .sgn = 0x0},
    {.hi = 0xbee23afc0853b6a8, .lo = 0xa89782c20df350c2,  .ex = -7, .sgn = 0x0},
    {.hi = 0xfe054587e01f1e2b, .lo = 0xf6d3a69bd5eab72f,  .ex = -7, .sgn = 0x0},
    {.hi = 0x9e75221a352ba751, .lo = 0x452b7ea62f2198ea,  .ex = -6, .sgn = 0x0},
    {.hi = 0xbdc8d83ead88d518, .lo = 0x7faa638b5e00ee90,  .ex = -6, .sgn = 0x0},
    {.hi = 0xdcfe013d7c8cbfc5, .lo = 0x632dbac46f30d009,  .ex = -6, .sgn = 0x0},
    {.hi = 0xfc14d873c1980236, .lo = 0xc7e09e3de453f5fc,  .ex = -6, .sgn = 0x0},
    {.hi = 0x8d86cc491ecbfe03, .lo = 0xf1776453b7e82558,  .ex = -5, .sgn = 0x0},
    {.hi = 0x9cf43dcff5eafd2f, .lo = 0x2ad90155c8a7236a,  .ex = -5, .sgn = 0x0},
    {.hi = 0xac52dd7e4726a456, .lo = 0xa47a963a91bb3018,  .ex = -5, .sgn = 0x0},
    {.hi = 0xbba2c7b196e7e224, .lo = 0xe7950f7252c163cf,  .ex = -5, .sgn = 0x0},
    {.hi = 0xcae41876471f5bde, .lo = 0x91d00a417e330f8e,  .ex = -5, .sgn = 0x0},
    {.hi = 0xda16eb88cb8df5fb, .lo = 0x28a63ecfb66e94c0,  .ex = -5, .sgn = 0x0},
    {.hi = 0xe93b5c56d85a9083, .lo = 0xce2992bfea38e76b,  .ex = -5, .sgn = 0x0},
    {.hi = 0xf85186008b1532f9, .lo = 0xe64b8b7759978998,  .ex = -5, .sgn = 0x0},
    {.hi = 0x83acc1acc7238978, .lo = 0x5a5333c45b7f442e,  .ex = -4, .sgn = 0x0},
    {.hi = 0x8b29b7751bd7073b,  .lo = 0x2e0b9ee992f2372,  .ex = -4, .sgn = 0x0},
    {.hi = 0x929fb17850a0b7be, .lo = 0x5b4d3807660516a4,  .ex = -4, .sgn = 0x0},
    {.hi = 0x9a0ebcb0de8e848e, .lo = 0x2c1bb082689ba814,  .ex = -4, .sgn = 0x0},
    {.hi = 0xa176e5f5323781d2, .lo = 0xdcf935996c92e8d4,  .ex = -4, .sgn = 0x0},
    {.hi = 0xa8d839f830c1fb40, .lo = 0x4c7343517c8ac264,  .ex = -4, .sgn = 0x0},
    {.hi = 0xb032c549ba861d83, .lo = 0x774e27bc92ce3373,  .ex = -4, .sgn = 0x0},
    {.hi = 0xb78694572b5a5cd3, .lo = 0x24cdcf68cdb2067c,  .ex = -4, .sgn = 0x0},
    {.hi = 0xbed3b36bd8966419, .lo = 0x7c0644d7d9ed08b4,  .ex = -4, .sgn = 0x0},
    {.hi = 0xc61a2eb18cd907a1, .lo = 0xe5a1532f6d5a1ac1,  .ex = -4, .sgn = 0x0},
    {.hi = 0xcd5a1231019d66d7, .lo = 0x761e3e7b171e44b2,  .ex = -4, .sgn = 0x0},
    {.hi = 0xd49369d256ab1b1f, .lo = 0x9e9154e1d5263cda,  .ex = -4, .sgn = 0x0},
    {.hi = 0xdbc6415d876d0839, .lo = 0x3e33c0c9f8824f54,  .ex = -4, .sgn = 0x0},
    {.hi = 0xe2f2a47ade3a18a8, .lo = 0xa0bf7c0b0d8bb4ef,  .ex = -4, .sgn = 0x0},
    {.hi = 0xea189eb3659aeaeb, .lo = 0x93b2a3b21f448259,  .ex = -4, .sgn = 0x0},
    {.hi = 0xf1383b7157972f48, .lo = 0x543fff0ff4f0aaf1,  .ex = -4, .sgn = 0x0},
    {.hi = 0xf85186008b153302, .lo = 0x5e4b8b7759978993,  .ex = -4, .sgn = 0x0},
    {.hi = 0xff64898edf55d548, .lo = 0x428ccfc99271dffa,  .ex = -4, .sgn = 0x0},
    {.hi = 0x8338a89652cb714a, .lo = 0xb247eb86498c2ce7,  .ex = -3, .sgn = 0x0},
    {.hi = 0x86bbf3e68472cb2f,  .lo = 0xb8bd20615747126,  .ex = -3, .sgn = 0x0},
    {.hi = 0x8a3c2c233a156341, .lo = 0x9027c74fe0e6f64f,  .ex = -3, .sgn = 0x0},
    {.hi = 0x8db956a97b3d0143, .lo = 0xf023472cd739f9e1,  .ex = -3, .sgn = 0x0},
    {.hi = 0x913378c852d65be6, .lo = 0x977e3013d10f7525,  .ex = -3, .sgn = 0x0},
    {.hi = 0x94aa97c0ffa91a5d, .lo = 0x4ee3880fb7d34429,  .ex = -3, .sgn = 0x0},
    {.hi = 0x981eb8c723fe97f2, .lo = 0x1f1c134fb702d433,  .ex = -3, .sgn = 0x0},
    {.hi = 0x9b8fe100f47ba1d8,  .lo = 0x4b62af189fcba0d,  .ex = -3, .sgn = 0x0},
    {.hi = 0x9efe158766314e4f, .lo = 0x4d71827efe892fc8,  .ex = -3, .sgn = 0x0},
    {.hi = 0xa2695b665be8f338, .lo = 0x4eca87c3f0f06211,  .ex = -3, .sgn = 0x0},
    {.hi = 0xa5d1b79cd2af2aca, .lo = 0x8837986ceabfbed6,  .ex = -3, .sgn = 0x0},
    {.hi = 0xa9372f1d0da1bd10, .lo = 0x580eb71e58cd36e5,  .ex = -3, .sgn = 0x0},
    {.hi = 0xac99c6ccc1042e94, .lo = 0x3dd557528315838d,  .ex = -3, .sgn = 0x0},
    {.hi = 0xaff983853c9e9e40, .lo = 0x5f105039091dd7f5,  .ex = -3, .sgn = 0x0},
    {.hi = 0xb3566a13956a86f4, .lo = 0x471b1e1574d9fd55,  .ex = -3, .sgn = 0x0},
    {.hi = 0xb6b07f38ce90e463, .lo = 0x7bb2e265d0de37e1,  .ex = -3, .sgn = 0x0},
    {.hi = 0xba07c7aa01bd2648, .lo = 0x43f9d57b324bd05f,  .ex = -3, .sgn = 0x0},
    {.hi = 0xbd5c481086c848db, .lo = 0xbb596b5030403242,  .ex = -3, .sgn = 0x0},
    {.hi = 0xc0ae050a1abf56ad, .lo = 0x2f7f8c5fa9c50d76,  .ex = -3, .sgn = 0x0},
    {.hi = 0xc3fd03290648847d, .lo = 0x30480bee4cbbd698,  .ex = -3, .sgn = 0x0},
    {.hi = 0xc74946f4436a054e, .lo = 0xf4f5cb531201c0d3,  .ex = -3, .sgn = 0x0},
    {.hi = 0xca92d4e7a2b5a3ad, .lo = 0xc983a9c5c4b3b135,  .ex = -3, .sgn = 0x0},
    {.hi = 0xcdd9b173efdc1aaa, .lo = 0x8863e007c184a1e7,  .ex = -3, .sgn = 0x0},
    {.hi = 0xd11de0ff15ab18c6, .lo = 0xd88d83d4cc613f21,  .ex = -3, .sgn = 0x0},
    {.hi = 0xd45f67e44178c612, .lo = 0x5486e73c615158b4,  .ex = -3, .sgn = 0x0},
    {.hi = 0xd79e4a7405ff96c3, .lo = 0x1300c9be67ae5da0,  .ex = -3, .sgn = 0x0},
    {.hi = 0xdada8cf47dad236d, .lo = 0xdffb833c3409ee7e,  .ex = -3, .sgn = 0x0},
    {.hi = 0xde1433a16c66b14c, .lo = 0xde744870f54f0f18,  .ex = -3, .sgn = 0x0},
    {.hi = 0xe14b42ac60c60512, .lo = 0x4e38eb8092a01f06,  .ex = -3, .sgn = 0x0},
    {.hi = 0xe47fbe3cd4d10d5b, .lo = 0x2ec0f797fdcd125c,  .ex = -3, .sgn = 0x0},
    {.hi = 0xe7b1aa704e2ee240, .lo = 0xb40faab6d2ad0841,  .ex = -3, .sgn = 0x0},
    {.hi = 0xeae10b5a7ddc8ad8, .lo = 0x806b2fc9a8038790,  .ex = -3, .sgn = 0x0},
    {.hi = 0xee0de5055f63eb01, .lo = 0x90a33316df83ba5a,  .ex = -3, .sgn = 0x0},
    {.hi = 0xf1383b7157972f4a, .lo = 0xb43fff0ff4f0aaf1,  .ex = -3, .sgn = 0x0},
    {.hi = 0xf460129552d2ff41, .lo = 0xe62e3201bb2bbdce,  .ex = -3, .sgn = 0x0},
    {.hi = 0xf7856e5ee2c9b28a, .lo = 0x76f2a1b84190a7dc,  .ex = -3, .sgn = 0x0},
    {.hi = 0xfaa852b25bd9b833, .lo = 0xa6dbfa03186e0666,  .ex = -3, .sgn = 0x0},
    {.hi = 0xfdc8c36af1f15468,  .lo = 0xa3361bca696504a,  .ex = -3, .sgn = 0x0},
    {.hi = 0x8073622d6a80e631, .lo = 0xe897009015316073,  .ex = -2, .sgn = 0x0},
    {.hi = 0x82012ca5a68206d5, .lo = 0x8fde85afdd2bc88a,  .ex = -2, .sgn = 0x0},
    {.hi = 0x838dc2fe6ac868e7, .lo = 0x1a3fcbdef40100cb,  .ex = -2, .sgn = 0x0},
    {.hi = 0x851927139c871af8, .lo = 0x67bd00c38061c51f,  .ex = -2, .sgn = 0x0},
    {.hi = 0x86a35abcd5ba5901, .lo = 0x5481c3cbd925ccd2,  .ex = -2, .sgn = 0x0},
    {.hi = 0x882c5fcd7256a8c1, .lo = 0x39055a6598e7c29e,  .ex = -2, .sgn = 0x0},
    {.hi = 0x89b438149d4582f5, .lo = 0x34531dba493eb5a6,  .ex = -2, .sgn = 0x0},
    {.hi = 0x8b3ae55d5d30701a, .lo = 0xc63eab8837170480,  .ex = -2, .sgn = 0x0},
    {.hi = 0x8cc0696ea11b7b36, .lo = 0x94361c9a28d38a6a,  .ex = -2, .sgn = 0x0},
    {.hi = 0x8e44c60b4ccfd7dc, .lo = 0x1473aa01c7778679,  .ex = -2, .sgn = 0x0},
    {.hi = 0x8fc7fcf24517946a, .lo = 0x380cbe769f2c6793,  .ex = -2, .sgn = 0x0},
    {.hi = 0x914a0fde7bcb2d0e, .lo = 0xc429ed3aea197a60,  .ex = -2, .sgn = 0x0},
    {.hi = 0x92cb0086fbb1cf75, .lo = 0xa29d47c50b1182d0,  .ex = -2, .sgn = 0x0},
    {.hi = 0x944ad09ef4351af1, .lo = 0xa49827e081cb16ba,  .ex = -2, .sgn = 0x0},
    {.hi = 0x95c981d5c4e924ea, .lo = 0x45404f5aa577d6b4,  .ex = -2, .sgn = 0x0},
    {.hi = 0x974715d708e984dd, .lo = 0x6648d42840d9e6fb,  .ex = -2, .sgn = 0x0},
    {.hi = 0x98c38e4aa20c27d2, .lo = 0x846767ec990d7333,  .ex = -2, .sgn = 0x0},
    {.hi = 0x9a3eecd4c3eaa6ae, .lo = 0xdb3a7f6e6087b947,  .ex = -2, .sgn = 0x0},
    {.hi = 0x9bb93315fec2d790, .lo = 0x7f589fba0865790f,  .ex = -2, .sgn = 0x0},
    {.hi = 0x9d3262ab4a2f4e37, .lo = 0xa1ae6ba06846fae0,  .ex = -2, .sgn = 0x0},
    {.hi = 0x9eaa7d2e0fb87c35, .lo = 0xff472bc6ce648a7d,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa0218434353f1de4, .lo = 0xd493efa632530acc,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa197795027409daa, .lo = 0x1dd1d4a6df960357,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa30c5e10e2f613e4, .lo = 0x9bd9bd99e39a20b3,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa4803402004e865c, .lo = 0x31cbe0e8824116cd,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa5f2fcabbbc506d8, .lo = 0x68ca4fb7ec323d74,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa764b99300134d79,  .lo = 0xd04d10474301862,  .ex = -2, .sgn = 0x0},
    {.hi = 0xa8d56c396fc1684c,  .lo = 0x1eb067d578c4756,  .ex = -2, .sgn = 0x0},
    {.hi = 0xaa45161d6e93167b, .lo = 0x9b081cf72249f5b2,  .ex = -2, .sgn = 0x0},
    {.hi = 0xabb3b8ba2ad362a1, .lo = 0x1db6506cc17a01f5,  .ex = -2, .sgn = 0x0},
    {.hi = 0xad215587a67f0cdf, .lo = 0xe890422cb86b7cb1,  .ex = -2, .sgn = 0x0},
    {.hi = 0xae8dedfac04e5282, .lo = 0xac707b8ffc22b3e8,  .ex = -2, .sgn = 0x0},
    {.hi = 0xaff983853c9e9e3f, .lo = 0xc5105039091dd7f8,  .ex = -2, .sgn = 0x0},
    {.hi = 0xb1641795ce3ca978, .lo = 0xfaf915300e517393,  .ex = -2, .sgn = 0x0},
    {.hi = 0xb2cdab981f0f940b, .lo = 0xc857c77dc1df600f,  .ex = -2, .sgn = 0x0},
    {.hi = 0xb43640f4d8a5761f, .lo = 0xf5f080a71c34b25d,  .ex = -2, .sgn = 0x0},
    {.hi = 0xb59dd911aca1ec48, .lo = 0x1d2664cf09a0c1bf,  .ex = -2, .sgn = 0x0},
    {.hi = 0xb70475515d0f1c5e, .lo = 0x4c98c6b8be17818d,  .ex = -2, .sgn = 0x0},
    {.hi = 0xb86a1713c491aeaa, .lo = 0xd37ee2872a6f1cd6,  .ex = -2, .sgn = 0x0},
};

/* for 0 <= i < 64, T1_2[i] is a 128-bit nearest approximation of 2^(i/64),
   with error bounded by 2^-128 (both absolutely and relatively).
   Table generated by output_T1_2() from the accompanying dint.sage file. */
static const dint64_t T1_2[] = {
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0},
    {.hi = 0x8164d1f3bc030773, .lo = 0x7be56527bd14def5, .ex = 0, .sgn = 0x0},
    {.hi = 0x82cd8698ac2ba1d7, .lo = 0x3e2a475b46520bff, .ex = 0, .sgn = 0x0},
    {.hi = 0x843a28c3acde4046, .lo = 0x1af92eca13fd1582, .ex = 0, .sgn = 0x0},
    {.hi = 0x85aac367cc487b14, .lo = 0xc5c95b8c2154c1b2, .ex = 0, .sgn = 0x0},
    {.hi = 0x871f61969e8d1010, .lo = 0x3a1727c57b52a956, .ex = 0, .sgn = 0x0},
    {.hi = 0x88980e8092da8527, .lo = 0x5df8d76c98c67563, .ex = 0, .sgn = 0x0},
    {.hi = 0x8a14d575496efd9a, .lo = 0x80ca1d92c3680c2, .ex = 0, .sgn = 0x0},
    {.hi = 0x8b95c1e3ea8bd6e6, .lo = 0xfbe4628758a53c90, .ex = 0, .sgn = 0x0},
    {.hi = 0x8d1adf5b7e5ba9e5, .lo = 0xb4c7b4968e41ad36, .ex = 0, .sgn = 0x0},
    {.hi = 0x8ea4398b45cd53c0, .lo = 0x2dc0144c8783d4c6, .ex = 0, .sgn = 0x0},
    {.hi = 0x9031dc431466b1dc, .lo = 0x775814a8494e87e2, .ex = 0, .sgn = 0x0},
    {.hi = 0x91c3d373ab11c336, .lo = 0xfd6d8e0ae5ac9d8, .ex = 0, .sgn = 0x0},
    {.hi = 0x935a2b2f13e6e92b, .lo = 0xd339940e9d924ee7, .ex = 0, .sgn = 0x0},
    {.hi = 0x94f4efa8fef70961, .lo = 0x2e8afad12551de54, .ex = 0, .sgn = 0x0},
    {.hi = 0x96942d3720185a00, .lo = 0x48ea9b683a9c22c5, .ex = 0, .sgn = 0x0},
    {.hi = 0x9837f0518db8a96f, .lo = 0x46ad23182e42f6f6, .ex = 0, .sgn = 0x0},
    {.hi = 0x99e0459320b7fa64, .lo = 0xe43086cb34b5fcaf, .ex = 0, .sgn = 0x0},
    {.hi = 0x9b8d39b9d54e5538, .lo = 0xa2a817a2a3cc3f1f, .ex = 0, .sgn = 0x0},
    {.hi = 0x9d3ed9a72cffb750, .lo = 0xde494cf050e99b0b, .ex = 0, .sgn = 0x0},
    {.hi = 0x9ef5326091a111ad, .lo = 0xa0911f09ebb9fdd1, .ex = 0, .sgn = 0x0},
    {.hi = 0xa0b0510fb9714fc2, .lo = 0x192dc79edb0fd9a9, .ex = 0, .sgn = 0x0},
    {.hi = 0xa27043030c496818, .lo = 0x9b7a04ef80cfdea8, .ex = 0, .sgn = 0x0},
    {.hi = 0xa43515ae09e6809e, .lo = 0xd1db4831781e1ef, .ex = 0, .sgn = 0x0},
    {.hi = 0xa5fed6a9b15138ea, .lo = 0x1cbd7f621710701b, .ex = 0, .sgn = 0x0},
    {.hi = 0xa7cd93b4e9653569, .lo = 0x9ec5b4d5039f72af, .ex = 0, .sgn = 0x0},
    {.hi = 0xa9a15ab4ea7c0ef8, .lo = 0x541e24ec3531fa73, .ex = 0, .sgn = 0x0},
    {.hi = 0xab7a39b5a93ed337, .lo = 0x658023b2759e0079, .ex = 0, .sgn = 0x0},
    {.hi = 0xad583eea42a14ac6, .lo = 0x4980a8c8f59a2ec4, .ex = 0, .sgn = 0x0},
    {.hi = 0xaf3b78ad690a4374, .lo = 0xdf26101ccbb35033, .ex = 0, .sgn = 0x0},
    {.hi = 0xb123f581d2ac258f, .lo = 0x87d037e96d215d8e, .ex = 0, .sgn = 0x0},
    {.hi = 0xb311c412a9112489, .lo = 0x3ecf14dc798a519c, .ex = 0, .sgn = 0x0},
    {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn = 0x0},
    {.hi = 0xb6fd91e328d17791, .lo = 0x7165f0ddd541a5a, .ex = 0, .sgn = 0x0},
    {.hi = 0xb8fbaf4762fb9ee9, .lo = 0x1b879778566b65a2, .ex = 0, .sgn = 0x0},
    {.hi = 0xbaff5ab2133e45fb, .lo = 0x74d519d24593838c, .ex = 0, .sgn = 0x0},
    {.hi = 0xbd08a39f580c36be, .lo = 0xa8811fb66d0faf7a, .ex = 0, .sgn = 0x0},
    {.hi = 0xbf1799b67a731082, .lo = 0xe815d0abcbf0b851, .ex = 0, .sgn = 0x0},
    {.hi = 0xc12c4cca66709456, .lo = 0x7c457d59a50087b5, .ex = 0, .sgn = 0x0},
    {.hi = 0xc346ccda24976407, .lo = 0x20ec856128b83a42, .ex = 0, .sgn = 0x0},
    {.hi = 0xc5672a115506dadd, .lo = 0x3e2ad0c964dd9f37, .ex = 0, .sgn = 0x0},
    {.hi = 0xc78d74c8abb9b15c, .lo = 0xc13a2e3976c0277e, .ex = 0, .sgn = 0x0},
    {.hi = 0xc9b9bd866e2f27a2, .lo = 0x80e1f92a0511697e, .ex = 0, .sgn = 0x0},
    {.hi = 0xcbec14fef2727c5c, .lo = 0xf4907c8f45ebf6dd, .ex = 0, .sgn = 0x0},
    {.hi = 0xce248c151f8480e3, .lo = 0xe235838f95f2c6ed, .ex = 0, .sgn = 0x0},
    {.hi = 0xd06333daef2b2594, .lo = 0xd6d45c6559a4d502, .ex = 0, .sgn = 0x0},
    {.hi = 0xd2a81d91f12ae45a, .lo = 0x12248e57c3de4028, .ex = 0, .sgn = 0x0},
    {.hi = 0xd4f35aabcfedfa1f, .lo = 0x5921deffa6262c5b, .ex = 0, .sgn = 0x0},
    {.hi = 0xd744fccad69d6af4, .lo = 0x39a68bb9902d3fde, .ex = 0, .sgn = 0x0},
    {.hi = 0xd99d15c278afd7b5, .lo = 0xfe873deca3e12bac, .ex = 0, .sgn = 0x0},
    {.hi = 0xdbfbb797daf23755, .lo = 0x3d840d5a9e29aa64, .ex = 0, .sgn = 0x0},
    {.hi = 0xde60f4825e0e9123, .lo = 0xdd07a2d9e8466859, .ex = 0, .sgn = 0x0},
    {.hi = 0xe0ccdeec2a94e111, .lo = 0x65895048dd333ca, .ex = 0, .sgn = 0x0},
    {.hi = 0xe33f8972be8a5a51, .lo = 0x9bfe90795980eed, .ex = 0, .sgn = 0x0},
    {.hi = 0xe5b906e77c8348a8, .lo = 0x1e5e8f4a4edbb0ed, .ex = 0, .sgn = 0x0},
    {.hi = 0xe8396a503c4bdc68, .lo = 0x791790d0ac70c7de, .ex = 0, .sgn = 0x0},
    {.hi = 0xeac0c6e7dd24392e, .lo = 0xd02d75b3706e54fb, .ex = 0, .sgn = 0x0},
    {.hi = 0xed4f301ed9942b84, .lo = 0x600d2db6a64bfb12, .ex = 0, .sgn = 0x0},
    {.hi = 0xefe4b99bdcdaf5cb, .lo = 0x46561cf6948db913, .ex = 0, .sgn = 0x0},
    {.hi = 0xf281773c59ffb139, .lo = 0xe8980a9cc8f47a4b, .ex = 0, .sgn = 0x0},
    {.hi = 0xf5257d152486cc2c, .lo = 0x7b9d0c7aed980fc3, .ex = 0, .sgn = 0x0},
    {.hi = 0xf7d0df730ad13bb8, .lo = 0xfe90d496d60fb6eb, .ex = 0, .sgn = 0x0},
    {.hi = 0xfa83b2db722a033a, .lo = 0x7c25bb14315d7fcd, .ex = 0, .sgn = 0x0},
    {.hi = 0xfd3e0c0cf486c174, .lo = 0x853f3a5931e0ee03, .ex = 0, .sgn = 0x0},
};

/* for 0 <= i < 64, T2_2[i] is a 128-bit nearest approximation of 2^(i/2^12),
   with error bounded by 2^-128 (both absolutely and relatively).
   Table generated by output_T2_2() from the accompanying dint.sage file. */
static const dint64_t T2_2[] = {
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0},
    {.hi = 0x80058baf7fee3b5d, .lo = 0x1c718b38e549cb93, .ex = 0, .sgn = 0x0},
    {.hi = 0x800b179c82028fd0, .lo = 0x945e54e2ae18f2f0, .ex = 0, .sgn = 0x0},
    {.hi = 0x8010a3c708e73282, .lo = 0x2b96d62d51c15a07, .ex = 0, .sgn = 0x0},
    {.hi = 0x8016302f17467628, .lo = 0x3690dfe44d11d008, .ex = 0, .sgn = 0x0},
    {.hi = 0x801bbcd4afcacb08, .lo = 0xe23a986bd3e626f0, .ex = 0, .sgn = 0x0},
    {.hi = 0x802149b7d51ebefb, .lo = 0x7bdbadbc888aeb29, .ex = 0, .sgn = 0x0},
    {.hi = 0x8026d6d889ecfd69, .lo = 0xb904bbfb40d3a2b7, .ex = 0, .sgn = 0x0},
    {.hi = 0x802c6436d0e04f50, .lo = 0xff8ce94a6797b3ce, .ex = 0, .sgn = 0x0},
    {.hi = 0x8031f1d2aca39b43, .lo = 0xad9db772901d96b6, .ex = 0, .sgn = 0x0},
    {.hi = 0x80377fac1fe1e56a, .lo = 0x61cd0bffd7cfc683, .ex = 0, .sgn = 0x0},
    {.hi = 0x803d0dc32d464f85, .lo = 0x43456f71b96affd4, .ex = 0, .sgn = 0x0},
    {.hi = 0x80429c17d77c18ed, .lo = 0x49fc841afba9c3c6, .ex = 0, .sgn = 0x0},
    {.hi = 0x80482aaa212e9e95, .lo = 0x86f7b54f6c45c85e, .ex = 0, .sgn = 0x0},
    {.hi = 0x804db97a0d095b0c, .lo = 0x6c9f1f7d1efcfe68, .ex = 0, .sgn = 0x0},
    {.hi = 0x805348879db7e67d, .lo = 0x171eb1ceef1d1f28, .ex = 0, .sgn = 0x0},
    {.hi = 0x8058d7d2d5e5f6b0, .lo = 0x94d589f608ee4aa2, .ex = 0, .sgn = 0x0},
    {.hi = 0x805e675bb83f5f0f, .lo = 0x2ed38ab8472b2144, .ex = 0, .sgn = 0x0},
    {.hi = 0x8063f722477010a1, .lo = 0xb1652de1378af1a1, .ex = 0, .sgn = 0x0},
    {.hi = 0x8069872686241a12, .lo = 0xb4ad9233a0390cad, .ex = 0, .sgn = 0x0},
    {.hi = 0x806f17687707a7af, .lo = 0xe54ec5f966eb1872, .ex = 0, .sgn = 0x0},
    {.hi = 0x8074a7e81cc7036b, .lo = 0x4d204ecfc11f4aab, .ex = 0, .sgn = 0x0},
    {.hi = 0x807a38a57a0e94dc, .lo = 0x9bf3ef4d9be2d1e4, .ex = 0, .sgn = 0x0},
    {.hi = 0x807fc9a0918ae142, .lo = 0x7068ab2230585d13, .ex = 0, .sgn = 0x0},
    {.hi = 0x80855ad965e88b83, .lo = 0xa0cc0a49c10ea66b, .ex = 0, .sgn = 0x0},
    {.hi = 0x808aec4ff9d45430, .lo = 0x84099bf6830f2768, .ex = 0, .sgn = 0x0},
    {.hi = 0x80907e044ffb1984, .lo = 0x3aa8b9cbbc65a8ab, .ex = 0, .sgn = 0x0},
    {.hi = 0x80960ff66b09d765, .lo = 0xf7d88c0928ba3947, .ex = 0, .sgn = 0x0},
    {.hi = 0x809ba2264dada76a, .lo = 0x4a8a4f44bb703db6, .ex = 0, .sgn = 0x0},
    {.hi = 0x80a13493fa93c0d4, .lo = 0x6699dc50dd96b774, .ex = 0, .sgn = 0x0},
    {.hi = 0x80a6c73f74697897, .lo = 0x6e0472ed4ccfa2e0, .ex = 0, .sgn = 0x0},
    {.hi = 0x80ac5a28bddc4157, .lo = 0xba2dc7e0c72e51ba, .ex = 0, .sgn = 0x0},
    {.hi = 0x80b1ed4fd999ab6c, .lo = 0x25335719b6e6fd20, .ex = 0, .sgn = 0x0},
    {.hi = 0x80b780b4ca4f64df, .lo = 0x534dfa7417846aa4, .ex = 0, .sgn = 0x0},
    {.hi = 0x80bd145792ab3970, .lo = 0xfc41c5c2d5336ccc, .ex = 0, .sgn = 0x0},
    {.hi = 0x80c2a838355b1297, .lo = 0x34dc28baed8f3fde, .ex = 0, .sgn = 0x0},
    {.hi = 0x80c83c56b50cf77f, .lo = 0xb880575ea03548c1, .ex = 0, .sgn = 0x0},
    {.hi = 0x80cdd0b3146f0d11, .lo = 0x32c1f98704428c71, .ex = 0, .sgn = 0x0},
    {.hi = 0x80d3654d562f95ec, .lo = 0x890e222a5eb95372, .ex = 0, .sgn = 0x0},
    {.hi = 0x80d8fa257cfcf26e, .lo = 0x24628efd9ca9d59b, .ex = 0, .sgn = 0x0},
    {.hi = 0x80de8f3b8b85a0af, .lo = 0x3b13310f5ad57fb1, .ex = 0, .sgn = 0x0},
    {.hi = 0x80e4248f84783c87, .lo = 0x1a9dfefaeb616564, .ex = 0, .sgn = 0x0},
    {.hi = 0x80e9ba216a837f8c, .lo = 0x718d1151d109bf98, .ex = 0, .sgn = 0x0},
    {.hi = 0x80ef4ff140564116, .lo = 0x996709da2e25f04c, .ex = 0, .sgn = 0x0},
    {.hi = 0x80f4e5ff089f763e, .lo = 0xe0adc640acaa6b0b, .ex = 0, .sgn = 0x0},
    {.hi = 0x80fa7c4ac60e31e1, .lo = 0xd4eb5edc6b341283, .ex = 0, .sgn = 0x0},
    {.hi = 0x810012d47b51a4a0, .lo = 0x8ccd7223820719e3, .ex = 0, .sgn = 0x0},
    {.hi = 0x8105a99c2b191ce1, .lo = 0xf24ebd6eb9ca4292, .ex = 0, .sgn = 0x0},
    {.hi = 0x810b40a1d81406d4, .lo = 0xcef03ab14a66550, .ex = 0, .sgn = 0x0},
    {.hi = 0x8110d7e584f1ec6d, .lo = 0x4bf94297d1519822, .ex = 0, .sgn = 0x0},
    {.hi = 0x81166f673462756d, .lo = 0xd0d8372f966cf15e, .ex = 0, .sgn = 0x0},
    {.hi = 0x811c0726e9156760, .lo = 0xb97931db7b7be2ec, .ex = 0, .sgn = 0x0},
    {.hi = 0x81219f24a5baa59d, .lo = 0x6abd3b0eab9c7048, .ex = 0, .sgn = 0x0},
    {.hi = 0x812737606d023148, .lo = 0xdaf888e96508151a, .ex = 0, .sgn = 0x0},
    {.hi = 0x812ccfda419c2956, .lo = 0xdc8046821f46122e, .ex = 0, .sgn = 0x0},
    {.hi = 0x813268922638ca8b, .lo = 0x6846ad73a8d9027f, .ex = 0, .sgn = 0x0},
    {.hi = 0x813801881d886f7b, .lo = 0xe885724f14131287, .ex = 0, .sgn = 0x0},
    {.hi = 0x813d9abc2a3b9090, .lo = 0x83768490519df895, .ex = 0, .sgn = 0x0},
    {.hi = 0x8143342e4f02c405, .lo = 0x661b22b45e25de18, .ex = 0, .sgn = 0x0},
    {.hi = 0x8148cdde8e8ebdec, .lo = 0xf11430fef78c6ee, .ex = 0, .sgn = 0x0},
    {.hi = 0x814e67cceb90502c, .lo = 0x99775205944eadc4, .ex = 0, .sgn = 0x0},
    {.hi = 0x815401f968b86a87, .lo = 0x7de463a40d18261, .ex = 0, .sgn = 0x0},
    {.hi = 0x81599c6408b81a94, .lo = 0x8f4a0b6748df7960, .ex = 0, .sgn = 0x0},
    {.hi = 0x815f370cce408bc8, .lo = 0xe2404468cfe5ab9f, .ex = 0, .sgn = 0x0},
};

# ifdef CORE_MATH_POW
/* The following is a degree-9 polynomial generated by Sollya, with zero
   constant coefficient, which approximates log(1+z) for |z| < 0.0001221,
   see sollya/approximations_r2.sollya.
   The coefficients of largest degree are first.
   The relative error is bounded by 2^-128.316.
   Table generated by output_P2() in the accompanying dint.sage file.
*/
static const dint64_t P_2[] = {
    {.hi = 0xe38e3954a09e560e, .lo = 0x0, .ex = -4, .sgn = 0x0},
    {.hi = 0x800000399d09d767, .lo = 0x0, .ex = -3, .sgn = 0x1},
    {.hi = 0x9249249249248676, .lo = 0x0, .ex = -3, .sgn = 0x0},
    {.hi = 0xaaaaaaaaaaaa9fdd, .lo = 0x0, .ex = -3, .sgn = 0x1},
    {.hi = 0xcccccccccccccccc, .lo = 0xcccdc5fe0ef93b8d, .ex = -3, .sgn = 0x0},
    {.hi = 0x8000000000000000, .lo = 0x600135b960d8, .ex = -2, .sgn = 0x1},
    {.hi = 0xaaaaaaaaaaaaaaaa, .lo = 0xaaaaaaaaaaa77b5e, .ex = -2, .sgn = 0x0},
    {.hi = 0xffffffffffffffff, .lo = 0xfffffffffffe33ca, .ex = -2, .sgn = 0x1},
    {.hi = 0x8000000000000000, .lo = 0x0, .ex = 0, .sgn = 0x0},
};
# else
static const dint64_t P_2[] = {
    {.hi = 0x99df88a0430813ca, .lo = 0xa1cffb6e966a70f6, .ex = -4, .sgn = 0x0},
    {.hi = 0xaaa02d43f696c3e4, .lo = 0x4dbe754667b6bc48, .ex = -4, .sgn = 0x1},
    {.hi = 0xba2e7a1eaf856174, .lo = 0x70e5c5a5ebbe0226, .ex = -4, .sgn = 0x0},
    {.hi = 0xccccccb9ec017492, .lo = 0xf934e28d924e76d4, .ex = -4, .sgn = 0x1},
    {.hi = 0xe38e38e3807cfa4b, .lo = 0xc976e6cbd22e203f, .ex = -4, .sgn = 0x0},
    {.hi = 0xfffffffffff924cc,  .lo = 0x5b308e39fa7dfb5, .ex = -4, .sgn = 0x1},
    {.hi = 0x924924924924911d, .lo = 0x862bc3d33abb3649, .ex = -3, .sgn = 0x0},
    {.hi = 0xaaaaaaaaaaaaaaaa, .lo = 0x6637fd4b19743eec, .ex = -3, .sgn = 0x1},
    {.hi = 0xcccccccccccccccc, .lo = 0xccc2ca18b08fe343, .ex = -3, .sgn = 0x0},
    {.hi = 0xffffffffffffffff, .lo = 0xffffff2245823ae0, .ex = -3, .sgn = 0x1},
    {.hi = 0xaaaaaaaaaaaaaaaa, .lo = 0xaaaaaaaaa5c48b54, .ex = -2, .sgn = 0x0},
    {.hi = 0xffffffffffffffff, .lo = 0xffffffffffffebd8, .ex = -2, .sgn = 0x1},
    {.hi = 0x8000000000000000,                .lo = 0x0,  .ex = 0, .sgn = 0x0},
};
# endif /* CORE_MATH_POW */
#endif /* UINT128_T */

/* The following is a degree-7 polynomial generated by Sollya,
   which approximates exp(z) for |z| < 0.00016923,
   see sollya/approximations_r2.sollya.
   The coefficients of largest degree are first.
   The relative error is bounded by 2^-122.415.
   Table generated by output_Q2() in the accompanying dint.sage file.
*/
static const dint64_t Q_2[] = {
    {.hi = 0xd00d00cd98416862, .lo = 0x0, .ex = -13, .sgn = 0x0},
    {.hi = 0xb60b60b932146a54, .lo = 0x0, .ex = -10, .sgn = 0x0},
    {.hi = 0x8888888888888897, .lo = 0x0, .ex = -7, .sgn = 0x0},
    {.hi = 0xaaaaaaaaaaaaaaa3, .lo = 0x0, .ex = -5, .sgn = 0x0},
    {.hi = 0xaaaaaaaaaaaaaaaa, .lo = 0xaaaaaa6a1e0776ae, .ex = -3, .sgn = 0x0},
    {.hi = 0x8000000000000000, .lo = 0xc06f3cd29, .ex = -1, .sgn = 0x0},
    {.hi = 0x8000000000000000, .lo = 0x88, .ex = 0, .sgn = 0x0},
    {.hi = 0xffffffffffffffff, .lo = 0xffffffffffffffd0, .ex = -1, .sgn = 0x0},
};
