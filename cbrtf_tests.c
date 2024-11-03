#include <fenv.h>
#include <math.h>
#include <mpfr.h>
#include <stdio.h>

static float
ref_cbrtf_round (float x, int rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cbrt (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

static float
ref_cbrtf_downward (float x)
{
  return ref_cbrtf_round (x, MPFR_RNDD);
}

static float
ref_cbrtf (float x)
{
  return ref_cbrtf_round (x, MPFR_RNDN);
}

void ref_init (void)
{
  mpfr_set_emin (-148);
  mpfr_set_emax (128);
}


extern float cr_cbrtf (float);

static float
cr_cbrtf_round (float x, int rnd)
{
  int old_rnd = fegetround ();
  fesetround (rnd);
  float f = cr_cbrtf (x);
  fesetround (old_rnd);
  return f;
}

static float
cr_cbrtf_downward (float x)
{
  return cr_cbrtf_round (x, FE_DOWNWARD);
}

static float
cbrtf_round (float x, int rnd)
{
  int old_rnd = fegetround ();
  fesetround (rnd);
  float f = cbrtf (x);
  fesetround (old_rnd);
  return f;
}

static float
cbrtf_downward (float x)
{
  return cbrtf_round (x, FE_DOWNWARD);
}


int main (int argc, char *argv[])
{
  ref_init ();

  {
    float i = -0x4.18937p-12;
    float r_mpfr  = ref_cbrtf (i);
    float r_core  = cr_cbrtf (i);
    float r_glibc = cbrtf (i);

    printf ("mpfr=%a core=%a flibc=%a\n",
	    r_mpfr,
	    r_core,
	    r_glibc);
  }
}
