/************************************************************
 *
 * This utility filters the DEBiL (31-param) database
 * above or below a given threshold for some column
 *
 * Compile:  gcc filterMECI.c -o filterMECI -Wall -lm
 * 
 * Run: filterMECI catalog.15x15w10.meci catalog.15x15w10.sortM1.info  catalog.15x15w10.sortM1.list
 *      ln -s /home/jdevor/work/MECI/TrES_MECI_catalog/run3_15x15w10/fitdat
 *      /home/jdevor/work/MECI/showMECI  catalog.15x15w10.sortM1.list  catalog.15x15w10.sortM1.list.sm
 *
 ************************************************************/

/* 
Line 1 - column description:
1.  string  File name
2.  float   period
3.  float   Mass1  [Msol]
4.  float   Error in mass1 [Msol]*
5.  float   Mass2 [Msol]
6.  float   Error in mass2 [Msol]*
7.  float   Age [Gyr]
8.  float   Error in age [Gyr]*
9.  float   Score   (low is better, should be about unity)
10. float   combined absolute magnitude
11. float   combined absolute magnitude error
12. float   score of absolute magnitude  (like #8, but for the magnitudes)
13. float   combined observes magnitude **
14. float   Distance [pc]   (from the distance modulus) **
15. string  Right-ascension (RA) ***
16. string  Declination (Dec) ***
17. float   Proper motion RA [mas/year] ***
18. float   Proper motion Dec [mas/year] ***
19. float   Coordinate error of proper motion catalogue [arcseconds]

Line 2 - column description:   (repeated)
1.  string  "x-filter"
2.  string/float  "Mag=xxx"   absolute magnitude in this filter
3.  string/float  "mag=xxx"   observed magnitude in this filter

Note: *If the best-fit was at the edge of the search range, values will be upper/lower limits,
       and the errors are ill-defined. There are three cases:
       -1.0 : unknown - one of the other parameters (mass/age) is upper/lower limit 
       -2.0 : upper limit
       -3.0 : lower limit
     **Assumes that the light curve magnitudes haven't been normalized, and that there is negligible reddening and blending.
    ***These values are taken from the 2mass and/or USNO-B catalogues, they may not be known ("-9999.999").

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>

#define MAX_SIZE 10000


struct EBinfo
{
  char name[256] ;
  double M1, M2, pm, Jmag, Vmag ;
} ;



double sqr (double x)
{
  return (x * x) ;
}


// A comparator for a pair-o-EBinfo
// Returns: -1 if x < y
//           0 if x = y
//           1 if x > y
int comparEB (const void *x, const void *y)
{
  if (((struct EBinfo*)x)->M1 < ((struct EBinfo*)y)->M1)
    return (-1) ;
  else
    return (((struct EBinfo*)x)->M1 > ((struct EBinfo*)y)->M1) ;
}


int main(int argc, char **argv)
{
  FILE *fin, *foutInfo, *foutList ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x17, x18, x19 ;
  double y1, y2, y3, y4, y5, y6, y7, y8, M1, M2 ;
  char x1[256], x15[256], x16[256] ;
  int i, size = 0 ;
  struct EBinfo *arr ;

  if (argc != 4)
    {
      printf (" Usage:  %s  <MECI DB>  <output info>  <output list>\n", argv[0]) ;
      return (1) ;
    }

  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (2) ;
    }  

  foutInfo = fopen (argv[2], "wt") ;

  if (!foutInfo)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[2]) ;
      return (3) ;
    }  

  foutList = fopen (argv[3], "wt") ;

  if (!foutList)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[3]) ;
      return (4) ;
    }  

  arr = (struct EBinfo*) malloc(MAX_SIZE * sizeof(struct EBinfo)) ;


  if (!arr)
    {
      printf ("Error: couldn't malloc()\n") ;
      return (5) ;
    }

  
  while (18 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s %lf %lf %lf\n",
		       x1, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, x15, x16, &x17, &x18, &x19))
    {
      fscanf (fin, "J-filter: Mag=%lf mag=%lf   H-filter: Mag=%lf mag=%lf   K-filter: Mag=%lf mag=%lf   calculated V-filter: Mag=%lf mag=%lf\n",
	      &y1, &y2, &y3, &y4, &y5, &y6, &y7, &y8) ;
      

      for (i = 1 ; ((x1[i-1] != '_') || (x1[i] != 'P')) && (i < strlen(x1)) ; i++) ;
      x2 = atof (&(x1[i+1])) ;
      
      if (x3 > x5)
	{
	  M1 = x3 ;
	  M2 = x5 ;
	}
      else
	{
	  M1 = x5 ;
	  M2 = x3 ;
	}

      if (y8 < 14.0)  ///xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	{
	  arr[size].M1 = M1 ;
	  arr[size].M2 = M2 ;
	  arr[size].pm = sqrt(sqr(x17) + sqr(x18)) ;
	  arr[size].Jmag = y2 ;
	  arr[size].Vmag = y8 ;
	  strcpy(arr[size].name, x1) ;
	  size++ ;
	}

      // printf ("%s   %lf    %lf %lf    %lf    %lf %lf\n",x1, x2, M1, M2, sqrt(sqr(x17) + sqr(x18)), y6, y8) ;      
    }

  qsort (arr, size, sizeof(struct EBinfo), comparEB) ;

  
  for (i = 0 ; i < size ; i++)
    {
      fprintf (foutInfo, "%30s  %3.2lf  %3.2lf   %5.1lf   %4.1lf  %4.1lf\n",arr[i].name, arr[i].M1, arr[i].M2, arr[i].pm, arr[i].Jmag, arr[i].Vmag) ;
      fprintf (foutList, "fitdat/%s.meci_fit\n", arr[i].name) ;
    }


  fclose (fin) ;
  fclose (foutInfo) ;
  fclose (foutList) ;
  free (arr) ;
  return (0) ;
}
