// compile:  gcc compare_log.c -o compare_log -Wall -lm

#include <stdio.h>
#include <math.h>

#define HIST_SIZE 30

#define VERBOSE

double sqr (double x)
{
  return (x * x) ;
}


void swap (double *x1, double *x2)
{
  double temp = *x1 ;
  
  *x1 = *x2 ;
  *x2 = temp ;
} 



int main (int argc, char **argv)
{
  FILE *finMECI, *finLog ;
  double mass1Log, mass2Log, mass1MECI, mass2MECI, dMass, dMassSwap ;
  int count, i, hist[HIST_SIZE] ;
  char str[256] ;
  
  if (argc != 3)
    {
      printf ("Usage: %s  <log file> <meci solution>\n", argv[0]) ;
      return (1) ;
    }
  
  finLog = fopen (argv[1], "rt") ;
  if (!finLog)
    {
      printf ("ERROR: Couldn't open the log file ('%s')\n", argv[1]) ;
      return (2) ;
    }
  
  finMECI = fopen (argv[2], "rt") ;
  if (!finMECI)
    {
      printf ("ERROR: Couldn't open the MECI file ('%s')\n", argv[2]) ;
      fclose (finLog) ;
      return (3) ;
    }
  
  for (i = 0 ; i < HIST_SIZE ; i++)
    hist[i] = 0 ;


  while ((3 == fscanf(finLog, "%d %lf %lf %*f %*f %*d %*f %*f\n", &count, &mass1Log, &mass2Log)) &&
	 (3 == fscanf(finMECI, "%s %lf %*f %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*s %*s %*f %*f %*f\n", str, &mass1MECI, &mass2MECI)))
    {
      dMass = sqrt(sqr((mass1Log - mass1MECI) / mass1Log) +
		   sqr((mass2Log - mass2MECI) / mass2Log)) ;
      dMassSwap = sqrt(sqr((mass1Log - mass2MECI) / mass1Log) +
		       sqr((mass2Log - mass1MECI) / mass2Log)) ;
      
      if (dMass > dMassSwap)
	{
	  dMass = dMassSwap ;
	  swap (&mass1MECI, &mass2MECI) ;
	}

      if (mass1Log < mass2Log)
	{
	  swap (&mass1Log, &mass2Log) ;
	  swap (&mass1MECI, &mass2MECI) ;
	}

#ifdef VERBOSE
      printf ("%d  %s  %f   [%f %f] [%f %f]\n", count, str, dMass, mass1Log, mass1MECI, mass2Log, mass2MECI) ;
#endif

    
      if (dMass <= 1.0e-9)
	i = 0 ;
      else
	{
	  i = (int)ceil((4.0 + log10(dMass)) * 5.0) ;
	  if (i < 0) i = 0 ;
	}
      
      if (i < HIST_SIZE)
	hist[i]++ ;
    }
  
  for (i = 0 ; i < HIST_SIZE ; i++)
    printf ("%d  %f  %d\n", i, pow (10.0, (i / 5.0) - 4.0), hist[i]) ;

  fclose (finLog) ;
  fclose (finMECI) ;
  return (0) ; 
}
