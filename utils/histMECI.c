// compile:  gcc histMECI.c -o histMECI -Wall -lm

#include <stdio.h>
#include <math.h>
#include <strings.h>

//#define SCORE_THRESHOLD 100.0   // interpolated reduced-chi-squared for LC + color

#define HIST_SIZE 100 // number of histogram bins

#define HIST_PERIOD_SCALE 5.0   // how many histogram bins per unit (e.g. day, solarMass, etc.)
#define HIST_LOG_PERIOD_SCALE 40.0  
#define HIST_LOG_PERIOD_OFFSET 1.0
#define HIST_MASS_SCALE 20.0
#define HIST_Q_SCALE 100.0

#define Q_EPSILON 1.0e-6

int main (int argc, char **argv)
{
  FILE *finMECI, *fout = 0 ;
  double period, mass1MECI, mass2MECI, score, tmp ;
  int i ;
  int histPeriod[HIST_SIZE+1], histLogPeriod[HIST_SIZE+1], histMass1[HIST_SIZE+1], histMass2[HIST_SIZE+1], histq[HIST_SIZE+1] ;
  char str[256] ;
  
  if (argc != 2)
    {
      printf ("Usage: %s <intput MECI file>\n", argv[0]) ;
      return (1) ;
    }
  
  finMECI = fopen (argv[1], "rt") ;
  if (!finMECI)
    {
      printf ("ERROR: Couldn't open the intput MECI file ('%s')\n", argv[1]) ;
      return (2) ;
    }

  if (argc >= 3)
    {
      fout = fopen (argv[2], "wt") ;
      if (!fout)
	{
	  printf ("ERROR: Couldn't open the output MECI-fit list file ('%s')\n", argv[2]) ;
	  fclose (finMECI) ;
	  return (3) ;
	} 
    }

  for (i = 0 ; i <= HIST_SIZE ; i++)
    {
      histPeriod[i] = 0 ;
      histLogPeriod[i] = 0 ;
      histMass1[i] = 0 ;
      histMass2[i] = 0 ;
      histq[i] = 0 ;
    }

  while (4 == fscanf(finMECI, "%s %lf %*f %lf %*f %*f %*f %lf %*f %*f %*f %*f %*f %*s %*s %*f %*f %*f\n", str, &mass1MECI, &mass2MECI, &score))
    //   if (score < SCORE_THRESHOLD)
    {
      fscanf (finMECI, "J-filter: Mag=%*f mag=%*f   H-filter: Mag=%*f mag=%*f   K-filter: Mag=%*f mag=%*f   calculated V-filter: Mag=%*f mag=%*f\n") ;
      
	for (i = 1 ; ((str[i-1] != '_') || (str[i] != 'P')) && (i < strlen(str)) ; i++) ;
	period = atof (&(str[i+1])) ;
	
	if (mass1MECI < mass2MECI)
	  {
	    tmp = mass1MECI ;
	    mass1MECI = mass2MECI ;
	    mass2MECI = tmp ;
	  }

	i = (int)floor(HIST_PERIOD_SCALE * period) ;
	if ((i >= 0) && (i <= HIST_SIZE))
	  histPeriod[i]++ ;

	i = (int)floor(HIST_LOG_PERIOD_SCALE * (HIST_LOG_PERIOD_OFFSET + log10(period))) ;
	if ((i >= 0) && (i <= HIST_SIZE))
	  histLogPeriod[i]++ ;

	i = (int)floor(HIST_MASS_SCALE * mass1MECI) ;
	if ((i >= 0) && (i <= HIST_SIZE))
	  histMass1[i]++ ;

	i = (int)floor(HIST_MASS_SCALE * mass2MECI) ;
	if ((i >= 0) && (i <= HIST_SIZE))
	  histMass2[i]++ ;

	i = (int)floor((HIST_Q_SCALE * mass2MECI / mass1MECI) -  Q_EPSILON) ;
	if ((i >= 0) && (i <= HIST_SIZE))
	  histq[i]++ ;

	if (i > 100)
	  printf (" qqqqq   %d  %lf  %lf\n", i, mass2MECI, mass1MECI) ;
      }
 
  if (!feof(finMECI))
    printf ("\n\n ERROR: Didn't reach the end of the file\n") ;

  for (i = 0 ; i <= HIST_SIZE ; i++)
    printf ("%f %d   %f %d   %f  %d %d   %f %d\n", 
	    i / HIST_PERIOD_SCALE, histPeriod[i],
	    (i / HIST_LOG_PERIOD_SCALE) - HIST_LOG_PERIOD_OFFSET, histLogPeriod[i],
	    i / HIST_MASS_SCALE , histMass1[i], histMass2[i],
	    i / HIST_Q_SCALE, histq[i]) ;

  fclose (fout) ;
  fclose (finMECI) ;
  return (0) ; 
}
