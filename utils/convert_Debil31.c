/************************************************************
 *
 * This utility converts debil database format (31-columns) to the MECI format
 * Warning: will create 3-column light curve files in the local directory
 *
 * Compile:  gcc convert_Debil31.c -o convert_Debil31 -Wall -lm
 * Run:       convert_Debil31 l2.dbl R > l2.dblc
 ************************************************************/

/* Column description:

1.  string  File name (with full path if not in local directory)
2.  float   Period (in days)
3.  float   Orbital eccentricity
4.  float   Absolute uncertainty in (3)
5.  float   Radius of large star (in units of semimajor axis)
6.  float   Absolute uncertainty in (5)
7.  float   Radius of small star (in units of semimajor axis)
8.  float   Absolute uncertainty in (7)
9.  float   Brightness of large star (magnitudes)
10. float   Absolute uncertainty in (9)
11. float   Brightness of small star (magnitudes)
12. float   Absolute uncertainty in (11)
13. float   Sine of inclination, sin(i)
14. float   Absolute uncertainty in (13)
15. float   Phased epoch of periastron
16. float   Absolute uncertainty in (15)
17. float   Argument of periastron (in degrees)
18. float   Absolute uncertainty in (17)
19. int     Number of used data points (not including outliers)
20. int     Number of outliers
21. float   Reduced chi squared of the best-fit model
22. float   Reduced chi squared of the average value
23. float   Reduced chi squared of a second order spline (parabolic fit within a sliding window)
24. float   Reduced chi squared of the sinusoidal best-fit
25. float   Significance of the secondary dip depth (in sigma)
26. float   Significance of the hump height (at midpoint between dips, in sigma)
27. float   Significance of the hump difference between the two humps (in sigma)
28. float   Waviness (see paper appendix for definition)
29. float   Scatter score (see paper appendix for definition)
30. float   Mean density (in grams per cm^3 ; see paper appendix for definition)
31. float   Max density (in grams per cm^3 ; see paper appendix for definition)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#define MIN_VALID_MAG_COLOR -100.0

#define UNKNOWN_PROPER_MOTION  -9999.9

int main(int argc, char **argv)
{
  FILE *fin ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  char x1[256] ;

  if (argc != 3)
    {
      printf (" Usage:  %s  <DEBiL database with color>  <LC filter>\n", argv[0]) ;
      return (1) ;
    }

  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (4) ;
    }  
  
  
  while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {
      printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
	      x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
	      x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
      
      printf ("%s  0  0 0  %.1f  %.1f  0\n", argv[2], UNKNOWN_PROPER_MOTION, UNKNOWN_PROPER_MOTION) ;
    }

fclose (fin) ;
return (0) ;
}
