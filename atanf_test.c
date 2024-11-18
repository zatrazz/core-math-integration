#include <math.h>
#include <fenv.h>
#include <stdio.h>
#include <signal.h>
#include <float.h>
#include <mpfr.h>
#include <setjmp.h>

static mpfr_rnd_t rnd2[] = { MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD };

static int rnd = 0; /* default is to nearest */

int ref_fesetround(int rounding_mode)
{
  rnd = rounding_mode;
  return 0;
}

void ref_init(void)
{
  mpfr_set_emin (-148);
  mpfr_set_emax (128);
}

float
ref_atanf (float x)
{
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, rnd2[rnd]);
  mpfr_subnormalize (y, inex, rnd2[rnd]);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

extern float cr_atanf (float);

static jmp_buf jb;

static void
sigfpe_handler (int sig, siginfo_t *info, void *closure)
{
  longjmp (jb, 0);
}

int main ()
{
#if 0
  {
    struct sigaction sa;
    sa.sa_sigaction = sigfpe_handler;
    sigemptyset (&sa.sa_mask);
    sa.sa_flags = SA_SIGINFO;
    sigaction (SIGFPE, &sa, NULL);
  }
  
  //ref_init ();

  feenableexcept (FE_UNDERFLOW);

  float x = FLT_MAX;
  int i = 0;

  while (1) {
    if (setjmp (jb) != 0)
      {
        printf ("[fail] %a (%f %g)\n", x, x, x);
	x = nexttowardf (x, 0.0);
	i++;
	continue;
      }
    
    float r = cr_atanf (x);
    x = nexttowardf (x, 0.0);
    if (r != 0x1.921fb6p+0)
      {
	printf ("[success (%d , %a)] %a (%f %g)\n", i, x, r, r, r);
	break;
      }
    i++;
  }
#else
  {
    float x = 0x1.e00a2cp+25;
    printf ("%a -> %a\n", x, cr_atanf (x));
  }
  {
    float x = nexttowardf (0x1.e00a2cp+25, INFINITY);
    for (int i=0; i<10; i++) {
      printf ("%a -> %a\n", x, cr_atanf (x));
      x = nexttowardf (x, INFINITY);
    }
  }
#endif
}
