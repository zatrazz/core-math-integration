/* Correctly-rounded hyperbolic cosine function for binary32 value.

Copyright (c) 2022-2023 Alexei Sibidanov.

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

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdint.h>
#ifdef CORE_MATH_SUPPORT_ERRNO
#include <errno.h>
#endif

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

/* __builtin_roundeven was introduced in gcc 10:
   https://gcc.gnu.org/gcc-10/changes.html,
   and in clang 17 */
#if ((defined(__GNUC__) && __GNUC__ >= 10) || (defined(__clang__) && __clang_major__ >= 17)) && (defined(__aarch64__) || defined(__x86_64__) || defined(__i386__) || defined(__powerpc64__))
#define HAS_BUILTIN_ROUNDEVEN
#endif

#if !defined(HAS_BUILTIN_ROUNDEVEN) && (defined(__GNUC__) || defined(__clang__)) && (defined(__AVX__) || defined(__SSE4_1__) || (__ARM_ARCH >= 8))
static inline double __builtin_roundeven(double x){
   double ix;
#if defined __AVX__
   __asm__("vroundsd $0x8,%1,%1,%0":"=x"(ix):"x"(x));
#elif __ARM_ARCH >= 8
   __asm__ ("frintn %d0, %d1":"=w"(ix):"w"(x));
#else /* __SSE4_1__ */
   __asm__("roundsd $0x8,%1,%0":"=x"(ix):"x"(x));
#endif
   return ix;
}
#define HAS_BUILTIN_ROUNDEVEN
#endif

#ifndef HAS_BUILTIN_ROUNDEVEN
#include <math.h>
/* round x to nearest integer, breaking ties to even */
static double
__builtin_roundeven (double x)
{
  double y = round (x); /* nearest, away from 0 */
  if (fabs (y - x) == 0.5)
  {
    /* if y is odd, we should return y-1 if x>0, and y+1 if x<0 */
    union { double f; uint64_t n; } u, v;
    u.f = y;
    v.f = (x > 0) ? y - 1.0 : y + 1.0;
    if (__builtin_ctz (v.n) > __builtin_ctz (u.n))
      y = v.f;
  }
  return y;
}
#endif

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

float cr_coshf(float x){
  static const double c[] =
    {1, 0x1.62e42fef4c4e7p-6, 0x1.ebfd1b232f475p-13, 0x1.c6b19384ecd93p-20};
  static const double ch[] =
    {1, 0x1.62e42fefa39efp-6, 0x1.ebfbdff82c58fp-13, 0x1.c6b08d702e0edp-20, 0x1.3b2ab6fb92e5ep-27,
     0x1.5d886e6d54203p-35, 0x1.430976b8ce6efp-43};
  static const uint64_t tb[] =
    {0x3fe0000000000000, 0x3fe059b0d3158574, 0x3fe0b5586cf9890f, 0x3fe11301d0125b51,
     0x3fe172b83c7d517b, 0x3fe1d4873168b9aa, 0x3fe2387a6e756238, 0x3fe29e9df51fdee1,
     0x3fe306fe0a31b715, 0x3fe371a7373aa9cb, 0x3fe3dea64c123422, 0x3fe44e086061892d,
     0x3fe4bfdad5362a27, 0x3fe5342b569d4f82, 0x3fe5ab07dd485429, 0x3fe6247eb03a5585,
     0x3fe6a09e667f3bcd, 0x3fe71f75e8ec5f74, 0x3fe7a11473eb0187, 0x3fe82589994cce13,
     0x3fe8ace5422aa0db, 0x3fe93737b0cdc5e5, 0x3fe9c49182a3f090, 0x3fea5503b23e255d,
     0x3feae89f995ad3ad, 0x3feb7f76f2fb5e47, 0x3fec199bdd85529c, 0x3fecb720dcef9069,
     0x3fed5818dcfba487, 0x3fedfc97337b9b5f, 0x3feea4afa2a490da, 0x3fef50765b6e4540};
  const double iln2 = 0x1.71547652b82fep+5;
  b32u32_u t = {.f = x};
  double z = x;
  uint32_t ax = t.u<<1;
  if(__builtin_expect(ax>0x8565a9f8u, 0)){ // |x| >~ 89.4
    if(ax>=0xff000000u) {
      if(ax<<8) return x + x; // nan
      return __builtin_inff(); // +-inf
    }
    float r = 2.0f*0x1.fffffep127f;
#ifdef CORE_MATH_SUPPORT_ERRNO
    errno = ERANGE;
#endif
    return r;
  }
  if(__builtin_expect(ax<0x7c000000u, 0)){ // |x| < 0.125
    if(__builtin_expect(ax<0x74000000u, 0)){ // |x| < 0x1p-11
      if(__builtin_expect(ax<0x66000000u, 0)) // |x| < 0x1p-24
	return __builtin_fmaf(__builtin_fabsf(x), 0x1p-25, 1.0f);
      return (0.5f*x)*x + 1.0f;
    }
    static const double cp[] =
      {0x1.fffffffffffe3p-2, 0x1.55555555723cfp-5, 0x1.6c16bee4a5986p-10, 0x1.a0483fc0328f7p-16};
    double z2 = z*z, z4 = z2*z2;
    return 1 + z2*((cp[0] + z2*cp[1]) + z4*(cp[2] + z2*(cp[3])));
  }
  double a = iln2*z, ia = __builtin_roundeven(a), h = a - ia, h2 = h*h;
  b64u64_u ja = {.f = ia + 0x1.8p52};
  int64_t jp = ja.u, jm = -jp;
  b64u64_u sp = {.u = tb[jp&31] + ((jp>>5)<<52)}, sm = {.u = tb[jm&31] + ((jm>>5)<<52)};
  double te = c[0] + h2*c[2], to = (c[1] + h2*c[3]);
  double rp = sp.f*(te + h*to), rm = sm.f*(te - h*to), r = rp + rm;
  float ub = r, lb = r  - 1.45e-10*r;
  if(__builtin_expect(ub != lb, 0)){
    const double iln2h = 0x1.7154765p+5, iln2l = 0x1.5c17f0bbbe88p-26;
    h = (iln2h*z - ia) + iln2l*z;
    h2 = h*h;
    te = ch[0] + h2*ch[2] + (h2*h2)*(ch[4] + h2*ch[6]);
    to = ch[1] + h2*(ch[3] + h2*ch[5]);
    r = sp.f*(te + h*to) + sm.f*(te - h*to);
    ub = r;
  }
  return ub;
}
