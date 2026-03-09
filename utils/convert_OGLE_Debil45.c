/************************************************************
 *
 * This utility converts debil database format (45-columns) to the MECI format
 * Warning: will create 3-column light curve files in the local directory
 *
 * Compile:  gcc convert_OGLE_Debil45.c -o convert_OGLE_Debil45 -Wall -lm
 * Run:      convert_OGLE_Debil45 qq.dbl > qq.dblc
 ************************************************************/

/* Column description:

1.  OGLE II bulge field number
2.  Light curve number (unique to field)
3.  Period
4.  Orbital eccentricity
5.  Absolute error in (4)
6.  Radius of large star (in units of semi-major axis)
7.  Absolute error in (6)
8.  Radius of small star (in units of semi-major axis)
9.  Absolute error in (8)
10. Brightness of large star (I-band magnitudes)
11. Absolute error in (10)
12. Brightness of small star (I-band magnitudes)
13. Absolute error in (12)
14. Sine of inclination - sin(i)
15. Absolute error in (14)
16. Epoch of perihelion (modulo the period, in units of the period)
17. Absolute error in (16)
18. Argument of perihelion (in degrees)
19. Absolute error in (18)
20. Number of used data points (not including outliers)
21. Number of outliers
22. Chi2 of model
23. Chi2 relative to the average value (should be worse than the model)
24. Chi2 relative to a second order spline
25. Chi2 relative to a sinusoidal best fit
26. Fitness score (position of (22) between (23) and (24))
27. Significance of the secondary dip depth (in sigma)
28. Significance of the hump height (at midpoint between dips in sigma)
29. Significance of the hump difference between the two humps
30. Waviness (see manual)
31. Scatter score (see manual)
32. Mean density (grams per cm^3 ; see manual)
33. Max density (grams per cm^3 ; see manual)
34. Extinction corrected brightness of large star (I-band magnitudes)
35. Extinction corrected brightness of small star (I-band magnitudes))
36. Extinction corrected total binary brightness (I-band magnitudes)
37. Extinction corrected I-band magnitude, from Udalski's catalog - similar to (34)
38. Extinction corrected color (V-I), from Udalski's catalog
39. Uncorrected I-band magnitude
40. Uncorrected color (V-I)
41. Absolute error in V-band magnitude
42. Absolute error in I-band magnitude
43. Binary right-ascension (RA)
44. Binary declination (Dec)
45. Inverse of probablity
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define LOCAL_OGLE_DIR  "/data/4shooter/OGLE/Bulge"

#define MIN_VALID_MAG_COLOR -100.0

#define UNKNOWN_PROPER_MOTION  -9999.9

int main(int argc, char **argv)
{
  FILE *fin ;
  double x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x19, x22, x23, x24, x25, x26, x27, x28, x29, x30 ;
  double x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x45 ;
  int x1, x2, x20, x21 ;
  char x43[256], x44[256], str[1000] ;

  if (argc != 2)
    {
      printf (" Usage:  %s  <DEBiL database with color>\n", argv[0]) ;
      return (1) ;
    }
 
  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (4) ;
    }  
  
  
  while (45 == fscanf (fin, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s %lf\n",
		       &x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, 
		       &x31, &x32, &x33, &x34, &x35, &x36, &x37, &x38, &x39, &x40, &x41, &x42, x43, x44, &x45))
    {
      sprintf (str, "awk '{printf (\"%%f %%f %%f\\n\", $1, $4, $5)}' %s/BUL_SC%d/bul_sc%d_%d.dat > bul_sc%d_%d.dat", LOCAL_OGLE_DIR, x1, x1, x2, x1, x2) ;
      system (str) ;

      printf ("bul_sc%d_%d.dat %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
	      x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
	      x17, x18, x19, x20, x21, x22, x23, x24, x25, x27, x28, x29, x30, x31, x32, x33) ;
      
      if ((x40 > MIN_VALID_MAG_COLOR) && (x41 > MIN_VALID_MAG_COLOR) && (x42 > MIN_VALID_MAG_COLOR))
	{
	  printf ("I  1  0 0  %.1f  %.1f  0\n", UNKNOWN_PROPER_MOTION, UNKNOWN_PROPER_MOTION) ;
	  printf ("V-I  %f %f\n", x40, sqrt((x41 * x41) + (x42 * x42))) ;
	}
      else
	printf ("I  0  0 0  %.1f  %.1f  0\n", UNKNOWN_PROPER_MOTION, UNKNOWN_PROPER_MOTION) ;
    }

fclose (fin) ;
return (0) ;
}
