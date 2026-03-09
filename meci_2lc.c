/************************************************************************
 * Finds the best fit of binary system to a set of given observables 
 * (see astro-ph/0510067 as well as the accompanied ReadMe file for further details)
 *
 * This version of MECI is designed to use TWO light curves, in different color-filters.
 * MECI fits both light curves together.
 *
 * IMPORTANT: The resulting magnitudes are the combined magnitudes of both light curves,
 *            applied as if they were a seperate pair of stars. This approach assumes that 
 *            the magnitudes and errors of the two light curves are similar.
 *
 * Written by Jonathan Devor (jdevor@cfa.harvard.edu)
 *    Harvard-Smithsonian Center for Astrophysics 
 *    60 Garden St., Cambridge, MA 02138 USA
 *
 * compile:  gcc meci_2lc.c debil_lib.c -o meci_2lc -lm -Wall -O4
 *           gcc meci_2lc.c debil_lib.c -o meci_2lc_allParam  -lm -Wall -O4           // unremark:  #define SIMPLEX_SEARCH_ALL_PARAMS
 *           gcc meci_2lc.c debil_lib.c -o meci_2lc_contour -lm -Wall -O4             // unremark:  #define CONTOUR_OUTPUT_SUFFIX
 *           gcc meci_2lc.c debil_lib.c -o meci_2lc_allParam_contour  -lm -Wall -O4   // unremark:  #define CONTOUR_OUTPUT_SUFFIX and SIMPLEX_SEARCH_ALL_PARAMS
 *
 * run:    meci_2lc q1.dblc LimbDarkTable.txt SunIsochrone.txt q1.meci q1.meci_err
 *
 * plot:
 *
 *   idl
 *   IDL> .run listMECI.txt.idl
 *   IDL> exit
 *
 *   convertPS2jpg.s 
 *   convertPS2bmp.s
 *   gzip *.ps *.contour
 *************************************************************************/

#include <math.h>
#include <stdio.h>
#include <values.h>      // for MAXFLOAT
#include <stdlib.h>      // for malloc()
#include <string.h>      // for memset(), memmove(), strcmp(), strlen() and strcat()
#include "debil_lib.h"

#define VERSION 1.4

#define DEBUG

#define MAX_STRLEN 512

#define MIN_AGE_GYEARS 0.2
#define MAX_AGE_GYEARS 10.0

#define MAX_NUM_LIMB_DARK_Teff 100 // needed since limb darkening coefficients are read once
#define MAX_NUM_STAR_GROUPS 1000   // needed since the groups are read once, and their number is unknown

#define OUTLIER_VARIANCE_LIMIT 8.0  // "sigma clipping" 

// The larger the basin, the fewer minima qualify and are searched.
#define BASIN_RADIUS 1    // basin area= (2*BASIN_RADIUS + 1)^2

#define DEFAULT_NUM_RESOLUTION_ELEMENTS 10
#define DEFAULT_LC_WEIGHTING 20.0

#define MIN_LIGHT_CURVE_DATA_POINTS 10

#define DETACHED_MAX_r1r2 0.8

//#define MIN_OBSERVABLE_DIP_DEPTH 0.01  // in magnitudes

//#define SET_ZERO_ECCENTRICITY          // assumes circular orbits

//#define SIMPLEX_SEARCH_ALL_PARAMS 100  // number of iterations, ignored if SET_ZERO_ECCENTRICITY

//#define DO_SPLINE_SCORE_CORRECTION    // use only when you are not confident in the LC data errors

#define MIN_BRACKET_SEARCH_SINI 0.0001  // the error in sin(i) is about half that, should be about (<mag error> * <r1>)
#define GOLDEN_RATIO 0.381966011       // = (3 - sqrt(5)) / 2

#define CONTOUR_OUTPUT_SUFFIX ".contour"  // remark this line to prevent contour files (and an IDL script) from being created.
#define IDL_SCRIPT_SUFFIX ".idl"
#define NUM_CONTOUR_LEVELS 12
#define MAX_CONTOUR_LABEL_LEN 8
#define CONTOUR_MIN_MASS 0.4
#define CONTOUR_MAX_MASS 2.5

#define WRITE_MODEL_FIT_LC_fA  "_A.meci_fit"  // suffix for the best fit model
#define WRITE_MODEL_DAT_LC_fA  "_A.meci_dat"  // suffix for the light curve dat
#define WRITE_MODEL_FIT_LC_fB  "_B.meci_fit"  // suffix for the best fit model
#define WRITE_MODEL_DAT_LC_fB  "_B.meci_dat"  // suffix for the light curve dat

// enumeration:
#define FILTER_U 0
#define FILTER_B 1
#define FILTER_V 2
#define FILTER_R 3
#define FILTER_I 4
#define FILTER_J 5
#define FILTER_H 6
#define FILTER_K 7
#define NUM_FILTERS 8

#define NUM_LOGG 11    // log(g) is {0=0.0 ; 1=0.5 ; 2=1.0 ; 3=1.5 ; ... ; 10=5.0}, where g is in CGS units

// reasons for having a "bad" score
#define MAX_VALID_SCORE               1000000.0
#define RESULT_NONDETACHED            1000001.0
#define RESULT_UNOBSERVABLE_DIP       1000002.0
#define RESULT_SKIPPED                1000003.0
#define RESULT_OUT_OF_BOUND_LIMB_DARK 1000004.0

// for parabolaFit()
#define MASS_PARABOLA_FIT_NUM 7   // 1 + (2 * #dimentions)
#define MAGNITUDE_PARABOLA_FIT_NUM 3
#define MAX_MATRIX_SIZE 7   // the maximum of the above two

#define UNCERTAINTY_UNKNOWN     -1.0
#define UNCERTAINTY_UPPER_LIMIT -2.0
#define UNCERTAINTY_LOWER_LIMIT -3.0


#define DIM 4  // The number of dimensions (don't change!)
// list of dimensions:
#define D_ECC   0  // Eccentricity of orbit (e)
#define D_SIN_I 1  // sin(inclination)
#define D_TMO   2  // (time / 2) - (omega / 4pi)
#define D_TPO   3  // (time / 2) + (omega / 4pi)


//------------- data structures ------------------

struct StarGroupType
{
  float Y, Z, alpha, age ;   // age is in Gyr
  int arrLen ;
  float *mass, *Teff, *radius ;  // mass and radius are in solar units
  float *Mag[NUM_FILTERS] ;      // abs. mag in various filters
  double bestSliceScore ;        // set only at the end of each age-slice iteration
} ;


struct LimbDarkType
{
  int numTeff[NUM_LOGG] ;
  unsigned long Teff[NUM_LOGG][MAX_NUM_LIMB_DARK_Teff] ;
  float coef_a[NUM_LOGG][NUM_FILTERS][MAX_NUM_LIMB_DARK_Teff] ;
  float coef_b[NUM_LOGG][NUM_FILTERS][MAX_NUM_LIMB_DARK_Teff] ;
} ;


struct ColorInfoType
{
  int filterPos, filterNeg ;   // e.g.  (B-V)   ; filterPos=B ; filterNeg=V
  double colorVal, colorErr ;
} ;


//--------- Math utilities ----------------------


// Returns the square of a given x
double sqr (double x)
{
  return (x * x) ;
}


// Swaps two given values
void swap (double *x1, double *x2)
{
  double tmp = *x1 ;
  *x1 = *x2 ;
  *x2 = tmp ;
}


// Given two stars with given magnitudes, returns their combined magnitude
double combineMag (double mag1, double mag2)
{
  return (-2.5 * log10(pow(10.0, -0.4 * mag1) + pow(10.0, -0.4 * mag2))) ;
}


// Given two stars with given magnitudes, returns their combined magnitude in both filters
// (i.e. counted as two pairs of stars)
double combineMag_2lc (double mag1_fA, double mag2_fA, double mag1_fB, double mag2_fB)
{
  return (-2.5 * log10(pow(10.0, -0.4 * mag1_fA) + pow(10.0, -0.4 * mag2_fA) + 
		       pow(10.0, -0.4 * mag1_fB) + pow(10.0, -0.4 * mag2_fB))) ;
}


// Given their combined magnitude (combMag) and the difference between the magnitudes (Mag2-Mag1 = deltaMag)
// returns the magnitudes of the two stars.
void calcMags (double combMag, double deltaMag, double *mag1, double *mag2)
{
  *mag1 = combMag + (2.5 * log10(1.0 + pow(10.0, -0.4 * deltaMag))) ;
  *mag2 = *mag1 + deltaMag ;
} 



// returns -1 upon error
int colorConvert (char c)
{
  switch (c)
    {
    case 'U': return (FILTER_U) ;
    case 'B': return (FILTER_B) ;
    case 'V': return (FILTER_V) ;
    case 'R': return (FILTER_R) ;
    case 'I': return (FILTER_I) ;
    case 'J': return (FILTER_J) ;
    case 'H': return (FILTER_H) ;
    case 'K': return (FILTER_K) ;
    }

  return (-1) ;
}


// Makes the "correction" to combObsMag_fA and combObsMag_fB, so that the distance modulus
// in both filters will be the same. This is similar to the color information.
void combObsMagCorrection (int g, int i, int j, int indexLCcolor_fA, int indexLCcolor_fB,
			   double combObsMag_orig_fA, int size_fA, double combObsMag_orig_fB, int size_fB,
			   struct StarGroupType *starGroup,
			   double *pcombObsMag_corrected_fA, double *pcombObsMag_corrected_fB)
{
  double combAbsMag_fA, combAbsMag_fB, obsMag_minus_absMag ;

  combAbsMag_fA = combineMag(starGroup[g].Mag[indexLCcolor_fA][i],
			     starGroup[g].Mag[indexLCcolor_fA][j]) ;
  combAbsMag_fB = combineMag(starGroup[g].Mag[indexLCcolor_fB][i],
			     starGroup[g].Mag[indexLCcolor_fB][j]) ;
  
  obsMag_minus_absMag = ((size_fA * (combObsMag_orig_fA - combAbsMag_fA)) + 
			 (size_fB * (combObsMag_orig_fB - combAbsMag_fB))) / (size_fA + size_fB) ;
  
  *pcombObsMag_corrected_fA = combAbsMag_fA + obsMag_minus_absMag ;
  *pcombObsMag_corrected_fB = combAbsMag_fB + obsMag_minus_absMag ;
}

//--------------------- IDL contour --------------------------------

// Creates an array of NUM_CONTOUR_LEVELS label strings, to be used in the contour plot. 
// They range from the minimum to the maximum range, with approximately logorithmic spacing.
// The deviation from exact logorithmic spacing is so to shorten the length of the string 
// (look elegant), while garanteeing that each contour spacing is less than doubled in size..
// returns: 0 = success ; otherwise failure
int makeContourLabels (char labels[NUM_CONTOUR_LEVELS][MAX_CONTOUR_LABEL_LEN], double Xmin, double Xmax)
{
  int i, j, exponent, digitMin, digitMax ;
  double Smin, Smax, dX ;
  char strNum[2] ;

  if (Xmin <= 0.0) return (1) ;
  if (Xmin >= Xmax) return (2) ;

  dX = log(Xmax/Xmin) / (NUM_CONTOUR_LEVELS-1) ;

  for (i = 0 ; i < NUM_CONTOUR_LEVELS ; i++)
    {
      Smin = Xmin * exp(i * dX) ;
      Smax = Xmin * exp((i+1) * dX) ;

      exponent = (int)floor(log10(Smax)) ;  // need floor(), since it might be negative

      Smin /= pow(10.0, exponent) ;
      Smax /= pow(10.0, exponent) ;

      labels[i][0] = 0 ;

      if (Smin == Smax) return (3) ;  // Warning: this doesn't catch all cases 

      if (exponent < 0)
	{
	  if (exponent < -MAX_CONTOUR_LABEL_LEN+2) return (4) ;

	  strcpy (labels[i], "0.") ;
	  for (j = exponent+1 ; j < 0 ; j++)
	    strcat (labels[i], "0") ;
	}

      while ((digitMin = (int)Smin) == (digitMax = (int)Smax))
	{
	  if (strlen(labels[i]) >= MAX_CONTOUR_LABEL_LEN-3) return (5) ;

	  strNum[0] = '0' + digitMin ;
	  strNum[1] = 0 ;
	  strcat (labels[i], strNum) ;

	  if (exponent == 0)
	    strcat (labels[i], ".") ;

	  exponent-- ;
	  Smin -= digitMin ;
	  Smin *= 10.0 ;
	  Smax -= digitMax ;
	  Smax *= 10.0 ;
	}

      strNum[0] = '0' + (int)ceil(0.5 * (digitMin + digitMax)) ;
      strNum[1] = 0 ;
      strcat (labels[i], strNum) ;

      if (strlen(labels[i]) + exponent >= MAX_CONTOUR_LABEL_LEN-1) return (6) ;

      for (j = exponent ; j > 0 ; j--)
	strcat (labels[i], "0") ;

      //      printf ("%f %f   %s\n", Xmin * exp(i * dX), Smax = Xmin * exp((i+1) * dX), labels[i]) ;
    }

  return (0) ;
}


// return: 0 = success ; otherwise- failure
int contourHeader (char *nameStr, float age, FILE *foutIDL, char *filenameContour,
		   double *contourMatrix, int matrixLen, int *subsetArr, int subsetLen)
{
  char labels[NUM_CONTOUR_LEVELS][MAX_CONTOUR_LABEL_LEN] ;
  double score, minScore = MAX_VALID_SCORE, maxScore = 0.0 ;
  int i, j ;

  for (i = 0 ; i < subsetLen ; i++)
    for (j = 0 ; j < subsetLen ; j++)
      {
	score = contourMatrix[(matrixLen * subsetArr[i]) + subsetArr[j]] ;
	
	if (score < minScore)
	  minScore = score ;

	if ((score > maxScore) && (score <= MAX_VALID_SCORE))
	  maxScore = score ;
      }

  if (makeContourLabels (labels, minScore, maxScore))
    return (1) ;

  //-------------

  fprintf (foutIDL, "device, file=\"%s.ps\"\n", filenameContour) ;
  fprintf (foutIDL, "openr, 1, '%s'\n", filenameContour) ;
  fprintf (foutIDL, "nelement = 0\n") ;
  fprintf (foutIDL, "readf, 1, nelement\n") ;
  fprintf (foutIDL, "m = dblarr(nelement)\n") ;
  fprintf (foutIDL, "readf, 1, m\n") ;
  fprintf (foutIDL, "z = dblarr(nelement, nelement)\n") ;
  fprintf (foutIDL, "readf, 1, z\n") ;
  fprintf (foutIDL, "contour, z, m, m, /fill, xtitle='Star1 Mass  [M!D!Mn!N]', ytitle='Star2 Mass  [M!D!Mn!N]', title='%s   Age=%.3fGyr', xrange=[%.4f,%.4f], yrange=[%.4f,%.4f], ",
	   nameStr, age, CONTOUR_MIN_MASS, CONTOUR_MAX_MASS, CONTOUR_MIN_MASS, CONTOUR_MAX_MASS) ;
  fprintf (foutIDL, "levels=[0") ;
  for (i = 0 ; i < NUM_CONTOUR_LEVELS ; i++)
    fprintf (foutIDL, ",%s", labels[i]) ;
  fprintf (foutIDL, "]\n") ;
  
  fprintf (foutIDL, "contour, z, m, m, /follow, /overplot, xrange=[%.4f,%.4f], yrange=[%.4f,%.4f], ",
	   CONTOUR_MIN_MASS, CONTOUR_MAX_MASS, CONTOUR_MIN_MASS, CONTOUR_MAX_MASS) ;
  fprintf (foutIDL, "levels=[0") ;
  for (i = 0 ; i < NUM_CONTOUR_LEVELS ; i++)
    fprintf (foutIDL, ",%s", labels[i]) ;
  fprintf (foutIDL, "], c_annotation=['0'") ;
  for (i = 0 ; i < NUM_CONTOUR_LEVELS ; i++)
    fprintf (foutIDL, ",'%s'", labels[i]) ;
  fprintf (foutIDL, "]\n") ;

  fprintf (foutIDL, "device, /close\n") ;
  fprintf (foutIDL, "close, 1\n\n") ;
  fflush(foutIDL) ;
  return (0) ;
}

//------------------- stellar limb darkening -------------------------------------------

// Loads the stellar limb darkening table (quadratic coefficients).
// returns: 0 = success, otherwise failure
int loadLimbDark(char *filenameLimbDark, struct LimbDarkType *plimbDark)
{
  FILE *fLimbDark = fopen (filenameLimbDark, "rt") ;
  char c ;
  float logg, filterU, filterB, filterV, filterR, filterI, filterJ, filterH, filterK ;
  int index_logg, numTeff, num, isNoRead = 1 ;
  unsigned long Teff ;

  if (!fLimbDark)
    {
      printf ("ERROR: Couldn't open limb dark data file '%s'.\n", filenameLimbDark) ;
      return (1) ;
    }

  memset (plimbDark->numTeff, 0, NUM_LOGG * sizeof(int)) ;

  //                              a/b         u   v   b   y   U  B  V  R  I  J  H  K
  while (11 == fscanf (fLimbDark, "%c %f %lu %*f %*f %*f %*f %f %f %f %f %f %f %f %f\n",
		       &c, &logg, &Teff, &filterU, &filterB, &filterV, &filterR, &filterI, &filterJ, &filterH, &filterK))
    {
      isNoRead = 0 ;
      index_logg = (int)(logg + logg + 0.5) ;  // log(g) is {0=0.0 ; 1=0.5 ; 2=1.0 ; 3=1.5 ; ... ; 10=5.0}
      numTeff = plimbDark->numTeff[index_logg] ;
    
      if ((c != 'a') || (index_logg < 0) || (index_logg >= NUM_LOGG) || (numTeff >= MAX_NUM_LIMB_DARK_Teff))
	{
	  fclose (fLimbDark) ;
	  printf ("ERROR: Invalid limb dark file format (a).\n") ;
	  return (2) ;
	}

      if (fabs(index_logg - logg - logg) > 0.0001)
	{
	  printf ("Warning: Skipping limb dark value (log(g)=%.2f)\n", logg) ;
	  fscanf (fLimbDark, "%*c %*f %*u %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f\n") ;  // skip (b) too
	  continue ;
	}
      
      plimbDark->numTeff[index_logg]++ ;
      plimbDark->Teff[index_logg][numTeff] = Teff ;
      plimbDark->coef_a[index_logg][FILTER_U][numTeff] = filterU ;
      plimbDark->coef_a[index_logg][FILTER_B][numTeff] = filterB ;
      plimbDark->coef_a[index_logg][FILTER_V][numTeff] = filterV ;
      plimbDark->coef_a[index_logg][FILTER_R][numTeff] = filterR ;
      plimbDark->coef_a[index_logg][FILTER_I][numTeff] = filterI ;
      plimbDark->coef_a[index_logg][FILTER_J][numTeff] = filterJ ;
      plimbDark->coef_a[index_logg][FILTER_H][numTeff] = filterH ;
      plimbDark->coef_a[index_logg][FILTER_K][numTeff] = filterK ;
 
      if ((numTeff > 0) && (plimbDark->Teff[index_logg][numTeff-1] >= Teff))
	{
	  fclose (fLimbDark) ;
	  printf ("ERROR: Limb dark Teff are not sorted (log(g)=%.2f  Teff=%lu)\n", index_logg * 0.5, Teff) ;
	  return (3) ;
	}
    
      num = fscanf (fLimbDark, "%c %f %lu %*f %*f %*f %*f %f %f %f %f %f %f %f %f\n",
		    &c, &logg, &Teff, &filterU, &filterB, &filterV, &filterR, &filterI, &filterJ, &filterH, &filterK) ;
      
      if ((num != 11) || (c != 'b') || (Teff != plimbDark->Teff[index_logg][numTeff]))
	{
	  fclose (fLimbDark) ;
	  printf ("ERROR: Invalid limb dark file format (b).\n") ;
	  return (4) ;
	}
      
      plimbDark->coef_b[index_logg][FILTER_U][numTeff] = filterU ;
      plimbDark->coef_b[index_logg][FILTER_B][numTeff] = filterB ;
      plimbDark->coef_b[index_logg][FILTER_V][numTeff] = filterV ;
      plimbDark->coef_b[index_logg][FILTER_R][numTeff] = filterR ;
      plimbDark->coef_b[index_logg][FILTER_I][numTeff] = filterI ;
      plimbDark->coef_b[index_logg][FILTER_J][numTeff] = filterJ ;
      plimbDark->coef_b[index_logg][FILTER_H][numTeff] = filterH ;
      plimbDark->coef_b[index_logg][FILTER_K][numTeff] = filterK ;
    }

  if (!feof(fLimbDark))
    printf ("Warning: Limb dark file may not have ended properly.\n") ;

  if (isNoRead)
    {
      fclose (fLimbDark) ;
      printf ("ERROR: No limb dark database lines to read.\n") ;
      return (5) ;
    }

  fclose (fLimbDark) ;
  return (0) ;
}


// linearly interpolates the temperature scale of the limb darkening array
// returns the lower index of the interval, and the weighting within the interval (*w)
// In case of an error (out of bound), returns a negative value
// assumes that the given array is sorted (monotonically rising)
int getAvrWeight (double x, unsigned long dataLimbDarkSlice[MAX_NUM_LIMB_DARK_Teff],
		  int dataLimbDarkSliceSize, double *w)
{
  int i, imax = dataLimbDarkSliceSize - 1, imin = 0 ;

  if ((dataLimbDarkSliceSize < 2) || (dataLimbDarkSlice[0] > x) || (dataLimbDarkSlice[dataLimbDarkSliceSize - 1] < x))
    return (-1) ;

  while ((imax - imin) > 1)
    {
      i = (imax + imin) / 2 ;
      
      if (x > dataLimbDarkSlice[i])
	imin = i ;
      else
	imax = i ;
    }

  if (imax <= imin)  // this should never happen
    return (-2) ;
     
  *w = (dataLimbDarkSlice[imax] - x) / (dataLimbDarkSlice[imax] - dataLimbDarkSlice[imin]) ;

  return (imin) ;
}



// returns:  0 = success ; otherwise failure 
// Assumes that the given array is sorted (monotonically rising)
int interpolateLimbDark (int indexLCcolor, double xLogg, double xTeff, struct LimbDarkType *plimbDark,
			 double *pcoef_a, double *pcoef_b)
{
  int iTeff1, iTeff2, iLogg ;
  double wTeff1, wTeff2, wLogg, wTeff1comp, wTeff2comp, wLoggComp ;
    
  xLogg *= 2.0 ;
  iLogg = (int)xLogg ;
  wLoggComp = xLogg - iLogg ;
  wLogg = 1.0 - wLoggComp ;

  if ((iLogg < 0) || ((iLogg+1) >= NUM_LOGG))
    return (1) ;

  iTeff1 = getAvrWeight (xTeff, plimbDark->Teff[iLogg], plimbDark->numTeff[iLogg], &wTeff1) ;
  if (iTeff1 < 0)
    return (100 - iTeff1) ;
  
  iTeff2 = getAvrWeight (xTeff, plimbDark->Teff[iLogg+1], plimbDark->numTeff[iLogg+1], &wTeff2) ;
  if (iTeff2 < 0)
    return (200 - iTeff2) ;
  
  wTeff1comp = 1.0 - wTeff1 ;
  wTeff2comp = 1.0 - wTeff2 ;

  *pcoef_a =
    (wLogg * ((wTeff1 * plimbDark->coef_a[iLogg][indexLCcolor][iTeff1]) +
	      (wTeff1comp * plimbDark->coef_a[iLogg][indexLCcolor][iTeff1+1]))) +
    (wLoggComp * ((wTeff2 * plimbDark->coef_a[iLogg+1][indexLCcolor][iTeff2]) +
		  (wTeff2comp * plimbDark->coef_a[iLogg+1][indexLCcolor][iTeff2+1]))) ;
  
  *pcoef_b =
    (wLogg * ((wTeff1 * plimbDark->coef_b[iLogg][indexLCcolor][iTeff1]) +
	      (wTeff1comp * plimbDark->coef_b[iLogg][indexLCcolor][iTeff1+1]))) +
    (wLoggComp * ((wTeff2 * plimbDark->coef_b[iLogg+1][indexLCcolor][iTeff2]) + 
		  (wTeff2comp * plimbDark->coef_b[iLogg+1][indexLCcolor][iTeff2+1]))) ;

  return (0) ;
}


//--------------- Isochrones ---------------------------

// for parsing
void nextLine (FILE *fin)
{
  // char c ;
  // printf ("\n((( ") ;
  // while ((c = fgetc(fin)) != '\n') 
  //   printf ("%c", c) ;
  // printf (" )))\n") ;

  while ((fgetc(fin) != '\n') && !feof(fin)) ;
}


// reads the input file and constructs the star isochrone database
// returns: the number of star groups, a negative value refers to an error.
int loadIsochrones (char *isochroneFilename, struct StarGroupType *starGroup)
{
  int numGroups = 0, arrLen, numRead, i ;
  float Z, Y, alpha, age ;
  double mass, logT, logg ;
  double M_V, UmB, BmV, VmR, VmI, VmJ, VmH, VmK ;
  FILE *fin = fopen (isochroneFilename, "rt") ;
      
  if (!fin)
    {
      printf ("ERROR: Couldn't open isochrone data file '%s'.\n", isochroneFilename) ;
      return (-1) ;
    }

  if (3 != fscanf (fin, "Z=%f Y=%f OS=%*f l/Hp=%*f [Fe/H]= %*f [Alpha/Fe]= %f\n", &Z, &Y, &alpha))
    {
      printf ("ERROR: Invalid isochrone file header.\n") ;
      fclose (fin) ;
      return (-2) ;
    }

  nextLine (fin) ;

  while (2 == fscanf (fin, "age(Gyr)= %f  %d points\n", &age, &arrLen))
    {
      if ((age < MIN_AGE_GYEARS) || (age > MAX_AGE_GYEARS))
	{
#ifdef DEBUG
	  printf ("Skip isochrone  (age = %f Gyr)\n", age) ;
#endif
	  for (i = 0 ; i < arrLen ; i++)
	    fscanf (fin, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*E %*E %*E\n") ;
	      
	  continue ; 
	}

      starGroup[numGroups].Y = Y ;
      starGroup[numGroups].Z = Z ;
      starGroup[numGroups].alpha = alpha ;

      starGroup[numGroups].age = age ;
      starGroup[numGroups].mass = (float*)malloc(arrLen * sizeof(float)) ;
      starGroup[numGroups].Teff = (float*)malloc(arrLen * sizeof(float)) ;
      starGroup[numGroups].radius = (float*)malloc(arrLen * sizeof(float)) ;

      if ((!starGroup[numGroups].mass) || (!starGroup[numGroups].Teff) || (!starGroup[numGroups].radius))
	{
	  printf ("ERROR: Not enough memory for loadIsochrones().\n") ;
	  fclose (fin) ;
	  return (-3) ;
	}
	  
      for (i = 0 ; i < NUM_FILTERS ; i++)
	{
	  starGroup[numGroups].Mag[i] = (float*)malloc(arrLen * sizeof(float)) ;

	  if (!starGroup[numGroups].Mag[i])
	    {
	      fclose (fin) ;
	      return (-4) ;
	    }
	}

      for (i = 0 ; i < arrLen ; i++)
	{
	  //  M/Msun     logT  logL/Ls   logg    Mv     U-B    B-V    V-R    V-I    V-J    V-H    V-K    V-L    V-M     #(x=-1)    #(x=1.35)    #(x=3)
	  numRead = fscanf (fin, "%lf %lf %*f %lf %lf %lf %lf %lf %lf %lf %lf %lf %*f %*f %*E %*E %*E\n",
			    &mass, &logT, &logg, &M_V, &UmB, &BmV, &VmR, &VmI, &VmJ, &VmH, &VmK) ;
	    
	  if (numRead != 11)
	    {
	      printf ("ERROR: Invalid isochrone file format.\n") ;
	      fclose (fin) ;
	      return (-5) ;
	    }

	  if ((i > 0) && (starGroup[numGroups].mass[i-1] >= mass))
	    {
#ifdef DEBUG
	      printf ("Warning: Non-sorted masses, skipping: (age=%f  mass=%f)\n", age, mass) ;
#endif
	      arrLen-- ;
	      i-- ;
	      continue ;
	    }

	  starGroup[numGroups].mass[i] = mass ;
	  starGroup[numGroups].Teff[i] = pow (10.0, logT) ;
	  starGroup[numGroups].radius[i] = 165.48974 * sqrt(mass * pow(10.0, -logg)) ;  // 165.48974 =~ sqrt(G*Msol) / Rsol
	  starGroup[numGroups].Mag[FILTER_U][i] = UmB + BmV + M_V ;
	  starGroup[numGroups].Mag[FILTER_B][i] = BmV + M_V ;
	  starGroup[numGroups].Mag[FILTER_V][i] = M_V ;
	  starGroup[numGroups].Mag[FILTER_R][i] = M_V - VmR ;
	  starGroup[numGroups].Mag[FILTER_I][i] = M_V - VmI ;
	  starGroup[numGroups].Mag[FILTER_J][i] = M_V - VmJ ;
	  starGroup[numGroups].Mag[FILTER_H][i] = M_V - VmH ;
	  starGroup[numGroups].Mag[FILTER_K][i] = M_V - VmK ;
	}

      starGroup[numGroups].arrLen = arrLen ;

      numGroups++ ;

      if (numGroups >= MAX_NUM_STAR_GROUPS)
	{  
	  printf ("ERROR: Too many isochrone groups.\n") ;
	  fclose (fin) ;
	  return (-6) ;
	}
    }

  if (!feof (fin))
    printf ("Warning: Isochrone file may not have ended properly.\n") ;

  if (numGroups == 0)
    printf ("ERROR: Did not read any groups.\n") ;
   
  fclose (fin) ;
  return (numGroups) ;
}



// Given a sorted array (arr) of length arrLen, returns a
// subset of it, that is equally spaced (or something close to it).
// Returns a sorted index array (subset) of the given length subsetLen.
void calcEqualSpaced (int *subset, int subsetLen, float *arr, int arrLen)
{
  int foundTune ;
  int i, j, subindexSmallGap = -1, subindexBigGap = -1, indexBigGapSplit = -1 ;
  float sizeSmallGap, sizeBigGapSplit, sizeSplit1, sizeSplit2 ;

  // 1. initialize:
  for (i = 0 ; i < subsetLen ; i++)
    subset[i] = (int)(0.5 + ((float)i * (arrLen-1) / (subsetLen-1))) ;

  // 2. maximize the smallest gap
  do
    {
      // 2.1. find the smallest gap
      sizeSmallGap = arr[arrLen-1] - arr[0] ;
      for (i = 1 ; i < subsetLen ; i++)  // scan gaps
	if (arr[subset[i]] - arr[subset[i-1]] < sizeSmallGap)
	  {
	    sizeSmallGap = arr[subset[i]] - arr[subset[i-1]] ;
	    subindexSmallGap = i ;
	  }

      // 2.2. find the gap that can be split, and the smaller half is maximal
      //      note that we do this search, assuming the above smallerst gap
      //      has already been removed
      sizeBigGapSplit = 0.0 ;
      for (i = 1 ; i < subsetLen ; i++)   // scan gaps
	{
	  if (i - 1 == subindexSmallGap)
	    {
	      for (j = subset[i-1] + 1 ; j < subset[i] ; j++)  // split within gap
		{
		  sizeSplit1 = arr[j] - arr[subset[i-2]] ;
		  sizeSplit2 = arr[subset[i]] - arr[j] ;

		  if (sizeSplit2 < sizeSplit1) // find the smaller split half
		    sizeSplit1 = sizeSplit2 ;

		  if (sizeSplit1 > sizeBigGapSplit)
		    {
		      sizeBigGapSplit = sizeSplit1 ;
		      subindexBigGap = i ;
		      indexBigGapSplit = j ;
		    }
		}
	    }
	  else if (i + 1 == subindexSmallGap)
	    {
	      for (j = subset[i-1] + 1 ; j < subset[i] ; j++)  // split within gap
		{
		  sizeSplit1 = arr[j] - arr[subset[i-1]] ;
		  sizeSplit2 = arr[subset[i+1]] - arr[j] ;

		  if (sizeSplit2 < sizeSplit1) // find the smaller split half
		    sizeSplit1 = sizeSplit2 ;

		  if (sizeSplit1 > sizeBigGapSplit)
		    {
		      sizeBigGapSplit = sizeSplit1 ;
		      subindexBigGap = i ;
		      indexBigGapSplit = j ;
		    }
		}
	    }
	  else if (i != subindexSmallGap)
	    {
	      for (j = subset[i-1] + 1 ; j < subset[i] ; j++)  // split within gap
		{
		  sizeSplit1 = arr[j] - arr[subset[i-1]] ;
		  sizeSplit2 = arr[subset[i]] - arr[j] ;

		  if (sizeSplit2 < sizeSplit1) // find the smaller split half
		    sizeSplit1 = sizeSplit2 ;

		  if (sizeSplit1 > sizeBigGapSplit)
		    {
		      sizeBigGapSplit = sizeSplit1 ;
		      subindexBigGap = i ;
		      indexBigGapSplit = j ;
		    }
		}
	    }
	}

      // 2.3. shift, if it improves things
      if (sizeBigGapSplit > sizeSmallGap)
	{
	  if (subindexBigGap < subindexSmallGap)
	    {
	      memmove (&subset[subindexBigGap+1], &subset[subindexBigGap],
		       (subindexSmallGap - subindexBigGap - 1) * sizeof(int)) ;
	      subset[subindexBigGap] = indexBigGapSplit ;
	    }
	  else //  (subindexBigGap >= subindexSmallGap)
	    {
	      memmove (&subset[subindexSmallGap], &subset[subindexSmallGap+1],
		       (subindexBigGap - subindexSmallGap - 1) * sizeof(int)) ;
	      subset[subindexBigGap-1] = indexBigGapSplit ;
	    }
	}
    }
  while (sizeBigGapSplit > sizeSmallGap) ;

  // 3. fine tune positions  (first order only)
  do
    {
      foundTune = 0 ;

      for (i = 1 ; i < subsetLen-1 ; i++)
	{
	  sizeBigGapSplit = 0.0 ;
	  for (j = subset[i-1]+1 ; j < subset[i+1] ; j++)
	    {
	      sizeSplit1 = arr[j] - arr[subset[i-1]] ;
	      sizeSplit2 = arr[subset[i+1]] - arr[j] ;

	      if (sizeSplit2 < sizeSplit1) // find the smaller split half
		sizeSplit1 = sizeSplit2 ;

	      if (sizeSplit1 > sizeBigGapSplit)
		{
		  sizeBigGapSplit = sizeSplit1 ;
		  indexBigGapSplit = j ;
		}
	    }

	  if (subset[i] != indexBigGapSplit)
	    {
	      subset[i] = indexBigGapSplit ;
	      foundTune = 1 ;
	    }
	}
    }
  while (foundTune) ;
}


//------------------- Likelihood of observation --------------------------------


// given the stars' radii (in arbitrary units, and their absolute magnitudes, 
// with calculate their maximum primary and secondary dip depth, and sees if this is
// smaller than the minimum that can be observed (MIN_OBSERVABLE_DIP_DEPTH).
#ifdef MIN_OBSERVABLE_DIP_DEPTH
int isUnobservableDip(double R1, double Mag1, double R2, double Mag2)
{
  double maxDip1, maxDip2 ;

  if (R1 < R2)
    {
      swap (&R1, &R2) ;
      swap (&Mag1, &Mag2) ;
    }

  maxDip1 = 2.5 * log10 (1.0 + pow(10.0, 0.4 * (Mag1 - Mag2))) ;
  maxDip2 = 2.5 * log10 (1.0 + (sqr(R2/R1) / (1.0 - sqr(R2/R1) + pow(10.0, 0.4 * (Mag1 - Mag2))))) ;

  return ((maxDip1 < MIN_OBSERVABLE_DIP_DEPTH) || (maxDip2 < MIN_OBSERVABLE_DIP_DEPTH)) ;
}
#endif

//-------------------------- 2lc -----------------------------------------

// Returns the reduced chi square of a given model, with given two light curves' data
double getChi2_2lc (float *time_fA, float *amp_fA, float *err_fA, int size_fA,
		    float *time_fB, float *amp_fB, float *err_fB, int size_fB,
		    double e, double r1, double r2, 
		    double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
		    double sin_i, double time0, double omega,
		    double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
		    double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB)
{
  double score_fA, score_fB ;
  
  score_fA = getChi2 (time_fA, amp_fA, err_fA, size_fA, e, r1, r2, mag1_fA, mag2_fA, sin_i, time0, omega,
		      quadA1_fA, quadB1_fA, quadA2_fA, quadB2_fA) ;
  
  score_fB = getChi2 (time_fB, amp_fB, err_fB, size_fB, e, r1, r2, mag1_fB, mag2_fB, sin_i, time0, omega,
		      quadA1_fB, quadB1_fB, quadA2_fB, quadB2_fB) ;

  return (((size_fA * score_fA) + (size_fB * score_fB)) / (size_fA + size_fB)) ;
}


// Returns the reduced chi square of a given model, with given two light curves' data
// Assumes a circular orbit (e=0, omega=0)
double getChi2_noEcc_2lc (float *time_fA, float *amp_fA, float *err_fA, int size_fA,
			  float *time_fB, float *amp_fB, float *err_fB, int size_fB,
			  double r1, double r2, 
			  double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
			  double sin_i, double time0,
			  double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
			  double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB)
{
  double score_fA, score_fB ;
  
  score_fA = getChi2noEcc (time_fA, amp_fA, err_fA, size_fA, r1, r2, mag1_fA, mag2_fA, sin_i, time0,
			   quadA1_fA, quadB1_fA, quadA2_fA, quadB2_fA) ;
  
  score_fB = getChi2noEcc (time_fB, amp_fB, err_fB, size_fB, r1, r2, mag1_fB, mag2_fB, sin_i, time0,
			   quadA1_fB, quadB1_fB, quadA2_fB, quadB2_fB) ;

  return (((size_fA * score_fA) + (size_fB * score_fB)) / (size_fA + size_fB)) ;
}


// Writes the two observed and model light curves
void writeLC_2lc (FILE *foutFit_fA, FILE *foutDat_fA, FILE *foutFit_fB, FILE *foutDat_fB,
		  float *time_fA, float *amp_fA, float *err_fA, int size_fA,
		  float *time_fB, float *amp_fB, float *err_fB, int size_fB,
		  double e, double r1, double r2, 
		  double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
		  double sin_i, double time0, double omega,
		  double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
		  double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB)
{
  writeLC (foutFit_fA, foutDat_fA, time_fA, amp_fA, err_fA, size_fA, e, r1, r2, mag1_fA, mag2_fA, sin_i, time0, omega,
	   quadA1_fA, quadB1_fA, quadA2_fA, quadB2_fA) ;

  writeLC (foutFit_fB, foutDat_fB, time_fB, amp_fB, err_fB, size_fB, e, r1, r2, mag1_fB, mag2_fB, sin_i, time0, omega,
	   quadA1_fB, quadB1_fB, quadA2_fB, quadB2_fB) ;
}


// Writes the two observed and model light curves
// Assumes a circular orbit (e=0, omega=0)
void  writeLC_noEcc_2lc (FILE *foutFit_fA, FILE *foutDat_fA, FILE *foutFit_fB, FILE *foutDat_fB,
			 float *time_fA, float *amp_fA, float *err_fA, int size_fA,
			 float *time_fB, float *amp_fB, float *err_fB, int size_fB,
			 double r1, double r2, 
			 double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
			 double sin_i, double time0,
			 double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
			 double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB)
{
  writeLC_noEcc (foutFit_fA, foutDat_fA, time_fA, amp_fA, err_fA, size_fA, r1, r2, mag1_fA, mag2_fA, sin_i, time0,
		 quadA1_fA, quadB1_fA, quadA2_fA, quadB2_fA) ;

  writeLC_noEcc (foutFit_fB, foutDat_fB, time_fB, amp_fB, err_fB, size_fB, r1, r2, mag1_fB, mag2_fB, sin_i, time0,
		 quadA1_fB, quadB1_fB, quadA2_fB, quadB2_fB) ;
}


//-------------------------------------------------------------------

// Prints the best-fit sin(i) result.
void printSini (double sini, char *errFilename, char *dataFilename)
{
  char strErr[MAX_STRLEN] ;

  sprintf (strErr, "sin(i) = %f", sini) ;

#ifdef SET_ZERO_ECCENTRICITY
  printError (errFilename, -1, dataFilename, "findBestSini_noEcc()", strErr) ;
#else
  printError (errFilename, -1, dataFilename, "findBestSini()", strErr) ;
#endif
}


// Runs the "golden section" search for the optimal sin(i)
// Returns the minimal reduced-chi-square.
// Note that Brent's method was tried, with poor results.
double findBestSini (float *time_fA, float *amp_fA, float *err_fA, int size_fA, 
		     float *time_fB, float *amp_fB, float *err_fB, int size_fB, double e,
		     double r1, double r2,
		     double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
		     double time0, double omega,
		     double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
		     double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB, 
		     FILE *foutFit_fA, FILE *foutDat_fA, FILE *foutFit_fB, FILE *foutDat_fB, char *errFilename, char *dataFilename)
{
  double sini_min, sini_max, sini_mid, sini_mid_chi2, sini_check_chi2, sini_check = MAXFLOAT ;
  double sini_1_chi2 = getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, e,
				    r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, 1.0, time0, omega,
				    quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;

  // check first if sin(i) = 1
  sini_max = 1.0 - MIN_BRACKET_SEARCH_SINI ;
  if (sini_1_chi2 < getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, e,
				 r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_max, time0, omega,
				 quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB))
    {
      if (foutFit_fA)
	writeLC_2lc (foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, time_fA, amp_fA, err_fA, size_fA,
		     time_fB, amp_fB, err_fB, size_fB, e, r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, 1.0, time0, omega,
		     quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;

      if (errFilename)
	printSini (1.0, errFilename, dataFilename) ;

      return (sini_1_chi2) ;
    }
  
  // bisection search
  sini_min = sqrt (1.0 - sqr(r1 + r2)) ;
  sini_mid = ((1.0 - GOLDEN_RATIO) * sini_max) + (GOLDEN_RATIO * sini_min) ; 
  sini_mid_chi2 = getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, e,
			       r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_mid, time0, omega, 
			       quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
  
  while ((sini_max - sini_min) > MIN_BRACKET_SEARCH_SINI)
    {
      if ((sini_max - sini_mid) < (sini_mid - sini_min))
	{
	  sini_check = ((1.0 - GOLDEN_RATIO) * sini_mid) + (GOLDEN_RATIO * sini_min) ;  
	  
	  sini_check_chi2 = getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, e,
					 r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,
					 sini_check, time0, omega, quadA1_fA,
					 quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
	  
	  if (sini_check_chi2 < sini_mid_chi2)  // in case equal, go for higher sin(i)
	    {
	      sini_max = sini_mid ;
	      sini_mid = sini_check ;
	      sini_mid_chi2 = sini_check_chi2 ;
	    }
	  else  // sini_check_chi2 >= sini_mid_chi2
	    sini_min = sini_check ;
	}
      else
	{
	  sini_check = ((1.0 - GOLDEN_RATIO) * sini_mid) + (GOLDEN_RATIO * sini_max) ;  
	  
	  sini_check_chi2 = getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, e, 
					 r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,
					 sini_check, time0, omega, quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
					 quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
	  
	  if (sini_check_chi2 <= sini_mid_chi2)  // in case equal, go for higher sin(i)  
	    {
	      sini_min = sini_mid ;
	      sini_mid = sini_check ;
	      sini_mid_chi2 = sini_check_chi2 ;
	    }
	  else  // sini_check_chi2 > sini_mid_chi2
	    sini_max = sini_check ;
	}
    }
  
  if (foutFit_fA)
    writeLC_2lc (foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
		 e, r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_mid, time0, omega,
		 quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
  
  if (errFilename)
    printSini (sini_mid, errFilename, dataFilename) ;

  return (sini_mid_chi2) ;
}


// same as findBestSini() only for circular orbits  (e=0, omega=0)
double findBestSini_noEcc (float *time_fA, float *amp_fA, float *err_fA, int size_fA,
			   float *time_fB, float *amp_fB, float *err_fB, int size_fB,
			   double r1, double r2, double mag1_fA, double mag1_fB,
			   double mag2_fA, double mag2_fB, double time0,
			   double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
			   double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB, 
			   FILE *foutFit_fA, FILE *foutDat_fA, FILE *foutFit_fB, FILE *foutDat_fB,
			   char *errFilename, char *dataFilename)
{
  double sini_min, sini_max, sini_mid, sini_mid_chi2, sini_check_chi2, sini_check = MAXFLOAT ;
  double sini_1_chi2 = getChi2_noEcc_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
					  r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, 1.0, time0, 
					  quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
					  quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
  
  // check first if sin(i) = 1
  sini_max = 1.0 - MIN_BRACKET_SEARCH_SINI ;
  if (sini_1_chi2 < getChi2_noEcc_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
				       r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_max, time0,
				       quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, 
				       quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB))
    {
      if (foutFit_fA)
	writeLC_noEcc_2lc (foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, time_fA, amp_fA, err_fA, size_fA, 
			   time_fB, amp_fB, err_fB, size_fB,
			   r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, 1.0, time0, 
			   quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
			   quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;

      if (errFilename)
	printSini (1.0, errFilename, dataFilename) ;
      
      return (sini_1_chi2) ;
    }
  
  // bisection search
  sini_min = sqrt (1.0 - sqr(r1 + r2)) ;
  sini_mid = ((1.0 - GOLDEN_RATIO) * sini_max) + (GOLDEN_RATIO * sini_min) ; 
  sini_mid_chi2 = getChi2_noEcc_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
				     r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_mid, time0,
				     quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
				     quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;

  while ((sini_max - sini_min) > MIN_BRACKET_SEARCH_SINI)
    {
      if ((sini_max - sini_mid) < (sini_mid - sini_min))
	{
	  sini_check = ((1.0 - GOLDEN_RATIO) * sini_mid) + (GOLDEN_RATIO * sini_min) ;  
	  
	  sini_check_chi2 = getChi2_noEcc_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
					       r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_check, time0,
					       quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, 
					       quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
	  
	  if (sini_check_chi2 < sini_mid_chi2)  // in case equal, go for higher sin(i)
	    {
	      sini_max = sini_mid ;
	      sini_mid = sini_check ;
	      sini_mid_chi2 = sini_check_chi2 ;
	    }
	  else  // sini_check_chi2 >= sini_mid_chi2
	    sini_min = sini_check ;
	}
      else
	{
	  sini_check = ((1.0 - GOLDEN_RATIO) * sini_mid) + (GOLDEN_RATIO * sini_max) ;  
	  
	  sini_check_chi2 = getChi2_noEcc_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
					       r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_check, time0,
					       quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, 
					       quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
	  
	  if (sini_check_chi2 <= sini_mid_chi2)  // in case equal, go for higher sin(i)
	    {
	      sini_min = sini_mid ;
	      sini_mid = sini_check ;
	      sini_mid_chi2 = sini_check_chi2 ;
	    }
	  else  // sini_check_chi2 > sini_mid_chi2
	    sini_max = sini_check ;
	}
    }

  if (foutFit_fA)
    writeLC_noEcc_2lc (foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, time_fA, amp_fA, err_fA, size_fA, 
		       time_fB, amp_fB, err_fB, size_fB, r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, sini_mid, time0,
		       quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
  
  if (errFilename)
    printSini (sini_mid, errFilename, dataFilename) ;
  
  return (sini_mid_chi2) ;
}


//xxxxxxxxxxxxxxxxx Nelder-Mead Rolling-Simplex Algorithm  [see Press (2002)] xxxxxxxxxxxxxxxxxxxxxxxxxx

// Part of the downhill simplex method [Press 2002]
// Checks out a given extrapolation factor for the highest point. If it is an improvement
// then it replaces the highest point (makes all the necessary updates).
// p[DIM+1][DIM] - a matrix of the (DIM+1) simplex points in the DIM-dimensional space
// y[DIM+1] - the score (function value) of the (DIM+1) simplex points
// psum[DIM] - the sum of the points' values in each of the DIM coordinates
// ihi - the index for the highest point
// fac - the factor for the highest point, where:
//  fac= 1   --> highest point
//  fac= 0   --> c.m. of the remaining points
//  fac= 0.5 --> mid point between the highest point and the c.m. of the remaining points
//  fac= -1  --> reflection of the highest point about the c.m. of the remaining points
// Returns the score (function value) for the extrapolated point
// Note that the arrays: p, y, psum - may get updated
double amotry (double p[DIM+1][DIM], double y[DIM+1], double psum[DIM], int ihi,
	       double fac, float *time_fA, float *amp_fA, float *err_fA, int size_fA,
	       float *time_fB, float *amp_fB, float *err_fB, int size_fB, 
	       double r1, double r2,
	       double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
	       double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB,
	       double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB)
{
  int j ;
  double ytry, ptry[DIM] ;
  double fac1 = (1.0 - fac) / DIM ;
  double fac2 = fac1 - fac ;

  for (j = 0 ; j < DIM ; j++)
    ptry[j] = (psum[j] * fac1) - (p[ihi][j] * fac2) ;

  ytry = getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, ptry[D_ECC], r1, r2,
		      mag1_fA, mag1_fB, mag2_fA, mag2_fB,  ptry[D_SIN_I], ptry[D_TPO] + ptry[D_TMO],
		      2.0 * M_PI * (ptry[D_TPO] - ptry[D_TMO]), quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
		      quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
  
  if (ytry < y[ihi])  // if the result is better than the worst, incorporate it
    {
      y[ihi] = ytry ;
      
      for (j = 0 ; j < DIM ; j++)
	psum[j] += (ptry[j] - p[ihi][j]) ;
      
      memcpy (p[ihi], ptry, DIM * sizeof(double)) ;
    }
  
  return (ytry) ;
}



// Based on the downhill simplex method [Press 2002]
// p[DIM+1][DIM] - a matrix of the (DIM+1) simplex points in the DIM-dimensional space
// y[DIM+1] - the score (function value) of the (DIM+1) simplex points
// numIterations = number of iteration to run the convergence algorithm
// Returns the index for the lowest (best) point in the simplex
int amoeba (double p[DIM+1][DIM], double y[DIM+1], float *time_fA, float *amp_fA, float *err_fA, int size_fA,
	    float *time_fB, float *amp_fB, float *err_fB, int size_fB, long numIterations,
	    double r1, double r2, 
	    double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB,
	    double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
	    double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB)
{
  int i, j, ihi, inhi, ilo ;
  double ytry, ysave, sum, psum[DIM] ;

  // Step 1: initializes the sum of values in each coordinate
  //      (to be used later for calculating the average, or center of mass)
  for (j = 0 ; j < DIM ; j++)
    {
      sum = p[0][j] ;
      for (i = 1 ; i <= DIM ; i++)
	sum += p[i][j] ;

      psum[j] = sum ;
    }

  
  do  // The main loop
    {
      // Step 2: find-
      //  ihi  = index of highest (worst) value
      //  ilo  = index of lowest (best) value

      // 2.1: sort the first two elements
      if (y[0] > y[1])
	{
	  ihi = 0 ;
	  inhi = ilo = 1 ;
	}
      else  // y[0] <= y[1]
	{
	  ihi = 1 ;
	  inhi = ilo = 0 ;
	}

      // 2.2: sort the rest
      for (i = 2 ; i <= DIM ; i++) 
	{
	  if (y[i] < y[ilo]) ilo = i ;
	  else 
	    if (y[i] > y[ihi]) 
	      {
		inhi = ihi ;
		ihi = i ;
	      }
	    else if (y[i] > y[inhi]) inhi = i ;
	}

      // Step 3: try various extrapolations (reflection, expansion and contraction)

      // 3.1: reflect the high point to the other side of the c.m. of the rest
      ytry = amotry (p, y, psum, ihi, -1.0, time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
		     r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, 
		     quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
      numIterations-- ;

      if (ytry <= y[ilo])        // 3.2: if the result was the best (lowest) so far, expand it more
	{
	  ytry = amotry (p, y, psum, ihi, 2.0, time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
			 r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,
			 quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
	  numIterations-- ;
	}
      else if (ytry >= y[inhi])  // 3.3: if the result is worse than the secnod-highest, try interpolating
	{
	  ysave = y[ihi] ;
	  ytry = amotry (p, y, psum, ihi, 0.5, time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			 r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, 
			 quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
	  numIterations-- ;
	  
	  if (ytry >= ysave)  //  3.3.1: if the result is still the worst, contract around the lowest (best) point
	    {
	      for (i = 0 ; i <= DIM ; i++)
		if (i != ilo) 
		  {
		    for (j = 0 ; j < DIM ; j++)
		      {
			psum[j] = 0.5 * (p[i][j] + p[ilo][j]) ;
			p[i][j] = psum[j] ;
		      }
		    
		    y[i] = getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, p[i][D_ECC],
					r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,  p[i][D_SIN_I],
					p[i][D_TPO] + p[i][D_TMO], 2.0 * M_PI * (p[i][D_TPO] - p[i][D_TMO]), 
					quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
		  }

	      numIterations -= DIM ;
	      
	      for (j = 0 ; j < DIM ; j++)   //  3.3.2: reinitialize the sum of values in each coordinate
		{
		  sum = p[0][j] ;
		  for (i = 1 ; i <= DIM ; i++)
		    sum += p[i][j] ;
		  
		  psum[j] = sum ;
		}
	    }
	}
    }
  while (numIterations > 0) ;

  return (ilo) ;
}



// Runs the "golden section" search for the optimal sin(i)
// Returns the minimal reduced-chi-square.
// Note that Brent's method was tried, with poor results.
#ifdef SIMPLEX_SEARCH_ALL_PARAMS
double findBestAllParam (float *time_fA, float *amp_fA, float *err_fA, int size_fA,
			 float *time_fB, float *amp_fB, float *err_fB, int size_fB, 
			 double e, double r1, double r2, 
			 double mag1_fA, double mag1_fB, double mag2_fA, double mag2_fB, double time0, double omega,
			 double quadA1_fA, double quadA1_fB, double quadB1_fA, double quadB1_fB, 
			 double quadA2_fA, double quadA2_fB, double quadB2_fA, double quadB2_fB, 
			 FILE *foutFit_fA, FILE *foutDat_fA, FILE *foutFit_fB, FILE *foutDat_fB, char *errFilename, char *dataFilename)
{
  char strErr[MAX_STRLEN] ;
  int i, ilo ;
  double p[DIM+1][DIM], y[DIM+1] ;
  double varySin_i = 0.5 * (1.0 - sqrt(1.0 - sqr(r1 + r2))) ;   // for an eclipse: sin(i) > sqrt(1.0 - sqr(r1 + r2))
  double vary_e = (e > 0.01) ? -0.5 * e : 0.01 ;   // makes sure the simplex isn't degenerate

  for (i = 0 ; i < (DIM+1) ; i++)
    {
      p[i][D_ECC] = e + ((i == 1) * vary_e) ;
      p[i][D_SIN_I] = 1.0 - ((i == 2) * varySin_i) ;
      p[i][D_TMO] = (0.5 * time0) - (omega / (4.0 * M_PI)) + ((i == 3) * 0.01) ;
      p[i][D_TPO] = (0.5 * time0) + (omega / (4.0 * M_PI)) + ((i == 4) * 0.1) ;
      
      y[i] =  getChi2_2lc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, p[i][D_ECC],
			   r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,  p[i][D_SIN_I],
			   p[i][D_TPO] + p[i][D_TMO], 2.0 * M_PI * (p[i][D_TPO] - p[i][D_TMO]), 
			   quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
    }
  
  ilo = amoeba(p, y, time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, SIMPLEX_SEARCH_ALL_PARAMS,
	       r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,
	       quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;

  if (foutFit_fA)
    writeLC_2lc (foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
		 p[ilo][D_ECC], r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB,
		 p[ilo][D_SIN_I], p[ilo][D_TPO] + p[ilo][D_TMO], 2.0 * M_PI * (p[ilo][D_TPO] - p[ilo][D_TMO]),
		 quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB) ;
  
  if (errFilename)
    {
      sprintf (strErr, "e=%f  sin(i)=%f  time0=%f  omega=%f",
	       p[ilo][D_ECC], p[ilo][D_SIN_I], p[ilo][D_TPO] + p[ilo][D_TMO], 2.0 * M_PI * (p[ilo][D_TPO] - p[ilo][D_TMO])) ;
      printError (errFilename, -1, dataFilename, "findBestAllParam()", strErr) ;
    }

  return (y[ilo]) ;
}
#endif

//-------------------------------------------------


double getScore (int g, int i, int j, int indexLCcolor_fA, int indexLCcolor_fB, double period, double splineScore,
		 float *time_fA, float *amp_fA, float *err_fA, int size_fA,
		 float *time_fB, float *amp_fB, float *err_fB, int size_fB,
		 double combObsMag_orig_fA, double combObsMag_orig_fB,
		 double e, double t0, double omega, float LCweighting,
		 struct LimbDarkType *plimbDark, struct ColorInfoType *colorInfoArr, int numColors,
		 struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS], 
		 FILE *foutFit_fA, FILE *foutDat_fA ,FILE *foutFit_fB, FILE *foutDat_fB,
		 char *errFilename, char *dataFilename)
{
  double mag1_fA, mag1_fB, mag2_fA, mag2_fB, combColor, score, inv_a, r1, r2 ;
  double xLogg, quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB ;
  double combObsMag_corrected_fA, combObsMag_corrected_fB ;
  int n ; 

  combObsMagCorrection (g, i,  j, indexLCcolor_fA, indexLCcolor_fB, combObsMag_orig_fA,	size_fA, combObsMag_orig_fB, size_fB, 
			starGroup, &combObsMag_corrected_fA, &combObsMag_corrected_fB) ;
  
  calcMags (combObsMag_corrected_fA, starGroup[g].Mag[indexLCcolor_fA][j] - starGroup[g].Mag[indexLCcolor_fA][i], 
	    &mag1_fA, &mag2_fA) ;
  calcMags (combObsMag_corrected_fB, starGroup[g].Mag[indexLCcolor_fB][j] - starGroup[g].Mag[indexLCcolor_fB][i],
	    &mag1_fB, &mag2_fB) ;

  // 0.237757 = R_sol * ((2pi)^2 * (day/sec)^-2 / (G*mass_sol))^(1/3)
  inv_a = 0.237757 * pow(period * period * (starGroup[g].mass[i] + starGroup[g].mass[j]), -0.3333333333333) ;
  r1 = inv_a * starGroup[g].radius[i] ;
  r2 = inv_a * starGroup[g].radius[j] ;	
  
#ifdef SET_ZERO_ECCENTRICITY
  if ((r1 + r2) > DETACHED_MAX_r1r2)
    return (RESULT_NONDETACHED) ;
#else
  if ((r1 + r2 + e) > DETACHED_MAX_r1r2)
    return (RESULT_NONDETACHED) ;
#endif

#ifdef MIN_OBSERVABLE_DIP_DEPTH
  if (isUnobservableDip(starGroup[g].radius[i], starGroup[g].Mag[indexLCcolor][i],
			starGroup[g].radius[j], starGroup[g].Mag[indexLCcolor][j]))
    return (RESULT_UNOBSERVABLE_DIP) ;
#endif

  //-----------  

  xLogg = 4.437542 + log10(starGroup[g].mass[i] / sqr(starGroup[g].radius[i])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor_fA, xLogg, starGroup[g].Teff[i], plimbDark, &quadA1_fA, &quadB1_fA))
    return (RESULT_OUT_OF_BOUND_LIMB_DARK) ;

  xLogg = 4.437542 + log10(starGroup[g].mass[j] / sqr(starGroup[g].radius[j])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor_fA, xLogg, starGroup[g].Teff[j], plimbDark, &quadA2_fA, &quadB2_fA))
    return (RESULT_OUT_OF_BOUND_LIMB_DARK) ;
  
  xLogg = 4.437542 + log10(starGroup[g].mass[i] / sqr(starGroup[g].radius[i])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor_fB, xLogg, starGroup[g].Teff[i], plimbDark, &quadA1_fB, &quadB1_fB))
    return (RESULT_OUT_OF_BOUND_LIMB_DARK) ;
  
  xLogg = 4.437542 + log10(starGroup[g].mass[j] / sqr(starGroup[g].radius[j])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor_fB, xLogg, starGroup[g].Teff[j], plimbDark, &quadA2_fB, &quadB2_fB))
    return (RESULT_OUT_OF_BOUND_LIMB_DARK) ;
  
  //-------------

  if (LCweighting != 0.0)
    {
#ifdef SET_ZERO_ECCENTRICITY
      score = findBestSini_noEcc (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
				  r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, t0 - (omega / 360.0), 
				  quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB, 
				  quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB,
				  foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, errFilename, dataFilename) ;
#else
#ifdef SIMPLEX_SEARCH_ALL_PARAMS
      score = findBestAllParam (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
				e, r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, t0, omega,
				quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
				quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB, 
				foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, errFilename, dataFilename) ;
#else
      score = findBestSini (time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    e, r1, r2, mag1_fA, mag1_fB, mag2_fA, mag2_fB, t0, omega,
			    quadA1_fA, quadA1_fB, quadB1_fA, quadB1_fB,
			    quadA2_fA, quadA2_fB, quadB2_fA, quadB2_fB,
			    foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, errFilename, dataFilename) ;
#endif
#endif
    }
  else score = 0.0 ;
  
  score /= splineScore ;  // normalize to take into account the model systematic error
  
  if (LCweighting < 0.0)
    LCweighting = -LCweighting * (size_fA + size_fB) ;
  
  score *= LCweighting ; // information weighting
  
  for (n = 0 ; n < numColors ; n++)
    {
      combColor = 
	combineMag(starGroup[g].Mag[colorInfoArr[n].filterPos][i],
		   starGroup[g].Mag[colorInfoArr[n].filterPos][j]) -
	combineMag(starGroup[g].Mag[colorInfoArr[n].filterNeg][i],
		   starGroup[g].Mag[colorInfoArr[n].filterNeg][j]) ;
     
      score += sqr((combColor - colorInfoArr[n].colorVal) / colorInfoArr[n].colorErr) ;
     
      // printf ("%f  ", sqr((combColor - colorInfoArr[n].colorVal) / colorInfoArr[n].colorErr)) ;
    }

  if ((score < 0.0) || (LCweighting + numColors <= 0.0))
    {
      printf ("ERROR: Invalid score.\n") ;
      return (MAX_VALID_SCORE) ;
    }
  

  score /= (LCweighting + numColors) ;

  if (score > MAX_VALID_SCORE)
    return (MAX_VALID_SCORE) ;

  return (score) ;
}


// Checks if (i,j) is the basin minimum of the contourMatrix subset.
// Basin  minimum is defined here as a point that is smaller (or equal) all the points around it, within a 
// given "radius" (BASIN_RADIUS). The maximal search is within an area of (2*BASIN_RADIUS + 1)^2 points.
// Returns:  1 = (i,j) is a basin minimum ; 0 = not a basin minimum
int isBasin (int i, int j, double *contourMatrix, int matrixLen, int *subsetArr, int subsetLen)
{
  double minimumVal = contourMatrix[(matrixLen * subsetArr[i]) + subsetArr[j]] ;
  int min_i = i - BASIN_RADIUS, max_i = i + BASIN_RADIUS ;
  int min_j = j - BASIN_RADIUS, max_j = j + BASIN_RADIUS ;

  if (minimumVal > MAX_VALID_SCORE)
    return (0) ;

  if (min_i < 0)
    min_i = 0 ;

  if (max_i >= subsetLen)
    max_i = subsetLen - 1 ;

  if (min_j < 0)
    min_j = 0 ;
  
  if (max_j >= subsetLen)
    max_j = subsetLen - 1 ;

  for (i = min_i ; i <= max_i ; i++)
    for (j = min_j ; j <= max_j ; j++)
      if (minimumVal > contourMatrix[(matrixLen * subsetArr[i]) + subsetArr[j]])
	return (0) ;
  
  return (1) ;
}



// Not the most efficient optimizing algorithm, but simple and good for when we are a small number of
// descrete steps away from the goal. Assumed to require fewer steps than the grid search.
// I'm not sure if the direction order needs to be optimized
double greedyFindMin (int g, int *pi_best, int *pj_best, int indexLCcolor_fA, int indexLCcolor_fB, double period, double splineScore,
		      float *time_fA, float *amp_fA, float *err_fA, int size_fA,
		      float *time_fB, float *amp_fB, float *err_fB, int size_fB, 
		      double combObsMag_fA, double combObsMag_fB,
		      double e, double t0, double omega, float LCweighting,
		      struct LimbDarkType *plimbDark, struct ColorInfoType *colorInfoArr, int numColors,
		      struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS], 
		      double *contourMatrix, int numLen)
{
  int isBetter, i = *pi_best, j = *pj_best ;
  double score, bestScore = contourMatrix[(numLen * i) + j] ;
  
  do
    {
      isBetter = 0 ;
      
    scan_left:
      if ((i > 0) && (contourMatrix[(numLen * (i-1)) + j] == RESULT_SKIPPED))
	{
	  score = getScore (g, i-1, j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			    time_fA, amp_fA, err_fA, size_fA,
			    time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * (i-1)) + j] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      i-- ;
	      isBetter = 1 ;
	      goto scan_left ;
	    }
       	}
      
    scan_down:
      if ((j > 0) && (contourMatrix[(numLen * i) + j - 1] == RESULT_SKIPPED))
	{
	  score = getScore (g, i, j-1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			    time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * i) + j - 1] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      j-- ;     
	      isBetter = 1 ;
	      goto scan_down ;
	    }
	}

    scan_right:
      if ((i < (numLen-1)) && (score = contourMatrix[(numLen * (i+1)) + j] == RESULT_SKIPPED)) 
	{
	  score = getScore (g, i+1, j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			    time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * (i+1)) + j] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      i++ ;     
	      isBetter = 1 ;
	      goto scan_right ;
	    }
	}
  
    scan_up:
      if ((j < (numLen-1)) && (score = contourMatrix[(numLen * i) + j + 1] == RESULT_SKIPPED))
	{
	  score = getScore (g, i, j+1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			    time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, 
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * i) + j + 1] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      j++ ;     
	      isBetter = 1 ;
	      goto scan_up ;
	    }
	}

      //----- diagonals  (not necessary) -----
     
    scan_left_up:
      if ((i > 0) && (j < (numLen-1)) && (contourMatrix[(numLen * (i-1)) + j + 1] == RESULT_SKIPPED))
	{
	  score = getScore (g, i-1, j+1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			    time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * (i-1)) + j + 1] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      i-- ;
	      j++ ;
	      isBetter = 1 ;
	      goto scan_left_up ;
	    }
       	}
      
    scan_left_down:
      if ((i > 0) && (j > 0) && (contourMatrix[(numLen * (i-1)) + j - 1] == RESULT_SKIPPED))
	{
	  score = getScore (g, i-1, j-1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			    time_fA, amp_fA, err_fA, size_fA, 
			    time_fB, amp_fB, err_fB, size_fB, 
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * (i-1)) + j - 1] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      i-- ;
	      j-- ;     
	      isBetter = 1 ;
	      goto scan_left_down ;
	    }
	}

    scan_right_down:
      if ((i < (numLen-1)) && (j > 0) && (score = contourMatrix[(numLen * (i+1)) + j - 1] == RESULT_SKIPPED)) 
	{
	  score = getScore (g, i+1, j-1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			    time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * (i+1)) + j - 1] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      i++ ;     
	      j-- ;
	      isBetter = 1 ;
	      goto scan_right_down ;
	    }
	}
      
    scan_right_up:
      if ((i < (numLen-1)) && (j < (numLen-1)) && (score = contourMatrix[(numLen * (i+1)) + j + 1] == RESULT_SKIPPED))
	{
	  score = getScore (g, i+1, j+1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			    time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			    combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting,
			    plimbDark, colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
	  
	  contourMatrix[(numLen * (i+1)) + j + 1] = score ;

	  if ((score < bestScore) && (score <= MAX_VALID_SCORE))
	    {
	      bestScore = score ;
	      i++ ;     
	      j++ ;   
	      isBetter = 1 ;
	      goto scan_right_up ;
	    }
	}
    }
  while (isBetter) ;

  *pi_best = i ;
  *pj_best = j ;

  return (bestScore) ;
}


//------------- paraboloid fitting ----------------


// Gauss-Jordan with full pivoting (adopted from "Numerical Recipes in C")
// returns:  0 = success, otherwise failure 
int gaussj (double a[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int n)
{
  int *indxc, *indxr, *ipiv ;
  int i, icol = 0, irow = 0, j, k, l, ll ;
  double big, dum, pivinv, tmp ;

  if (n > MAX_MATRIX_SIZE)
    return (1) ;

  indxc = (int*)malloc(n * sizeof(int)) ;
  indxr = (int*)malloc(n * sizeof(int)) ;
  ipiv = (int*)calloc(n, sizeof(int)) ;

  if (!indxc || !indxr || !ipiv)
    {
      printf ("ERROR: Memory allocation failure in gaussj().\n") ;
      if (ipiv) free (ipiv) ;
      if (indxr) free (indxr) ;
      if (indxc) free (indxc) ;
      return (2) ; 
    }

  for (i = 0 ; i < n ; i++) 
    {
      big = 0.0 ;

      for (j = 0 ; j < n ; j++)
	if (ipiv[j] != 1)
	  for (k = 0 ; k < n ; k++) 
	    {
	      if (ipiv[k] == 0)
		{
		  if (fabs(a[j][k]) >= big)
		    {
		      big = fabs(a[j][k]) ;
		      irow = j ;
		      icol = k ;
		    }
		}
	      else if (ipiv[k] > 1) 
		return (3) ; // Singular Matrix-1
	    }

      ipiv[icol]++ ;

      if (irow != icol)
	for (l = 0 ; l < n ; l++) 
	  {
	    tmp = a[irow][l] ;
	    a[irow][l] = a[icol][l] ;
	    a[icol][l] = tmp ;
	  }
      
      indxr[i] = irow ;
      indxc[i] = icol ;

      if (a[icol][icol] == 0.0) 
	return (4) ; // Singular Matrix-2
      
      pivinv = 1.0 / a[icol][icol] ;
      a[icol][icol] = 1.0 ; 

      for (l = 0 ; l < n ; l++)
	a[icol][l] *= pivinv ;

      for (ll = 0 ; ll < n ; ll++)
	if (ll != icol) 
	  {
	    dum = a[ll][icol] ;
	    a[ll][icol] = 0.0 ;
 
	    for (l = 0 ; l < n ; l++) 
	      a[ll][l] -= a[icol][l] * dum ;
	  }
    }
  
  for (l = n-1 ; l >= 0 ; l--)
    {
      if (indxr[l] != indxc[l])
	for (k = 0 ; k < n ; k++)
	  {
	    tmp = a[k][indxr[l]] ;
	    a[k][indxr[l]] = a[k][indxc[l]] ;
	    a[k][indxc[l]] = tmp ;
	  }
    }
  
  free (ipiv) ;
  free (indxr) ;
  free (indxc) ;
  return (0) ;
}



// returns:  0 = success, otherwise (partial) failure 
int parabolaFit (int best_g, int best_i, int best_j, int indexLCcolor_fA, int indexLCcolor_fB, 
		 double period, double splineScore,
		 float *time_fA, float *amp_fA, float *err_fA, int size_fA, 
		 float *time_fB, float *amp_fB, float *err_fB, int size_fB,
		 double combObsMag_fA, double combObsMag_fB,
		 double e, double t0, double omega, float LCweighting,
		 char *strRA, char *strDEC, float pmRA, float pmDEC, float ra_dec_err,
		 struct LimbDarkType *plimbDark, struct ColorInfoType *colorInfoArr, int numColors,
		 struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS],
		 int numGroups, char *errFilename, char *deb_filename, FILE *foutResults)
{
  int i, j ;
  double fitMatrix[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE] ;  // row num, col num
  double vec[MASS_PARABOLA_FIT_NUM] ;
  double coord[MASS_PARABOLA_FIT_NUM][4] ;   // 0=score 1=X  2=Y  3=Z  
  char strErr[MAX_STRLEN] ;
  double mass1, mass1Err, mass2, mass2Err, age, ageErr, combAbsMag, combAbsMagErr ;
  double minimumOfMassAge, minimumOfMag ;
  double combObsMag_fAB = combineMag (combObsMag_fA, combObsMag_fB) ;

  if ((best_i < 0) || (best_i >= starGroup[best_g].arrLen) ||
      (best_j < 0) || (best_j >= starGroup[best_g].arrLen) ||
      (best_g < 0) || (best_g >= numGroups))                      // shouldn't happen
    {
      printError (errFilename, 1, deb_filename, "parabolaFit()", "Out of bound parameter.") ;
      return (1) ;
    } 

  if ((best_i == 0) || (best_i == (starGroup[best_g].arrLen - 1)) ||
      (best_j == 0) || (best_j == (starGroup[best_g].arrLen - 1)) ||
      (best_g == 0) || (best_g == (numGroups - 1)))
    {
      if (best_i == 0)
	mass1Err = UNCERTAINTY_UPPER_LIMIT ;
      else
	if (best_i == (starGroup[best_g].arrLen - 1))
	  mass1Err = UNCERTAINTY_LOWER_LIMIT ;
	else
	  mass1Err = UNCERTAINTY_UNKNOWN ;

      if (best_j == 0)
	mass2Err = UNCERTAINTY_UPPER_LIMIT ;
      else
	if (best_j == (starGroup[best_g].arrLen - 1))
	  mass2Err = UNCERTAINTY_LOWER_LIMIT ;
	else
	  mass2Err = UNCERTAINTY_UNKNOWN ;
      
      if (best_g == 0) 
	ageErr = UNCERTAINTY_UPPER_LIMIT ;
      else 
	if (best_g == (numGroups - 1))
	  ageErr = UNCERTAINTY_LOWER_LIMIT ;
	else
	  ageErr = UNCERTAINTY_UNKNOWN ;

      printError (errFilename, 0, deb_filename, "parabolaFit()", "Best position is at the range edge.") ;

      goto mass_and_mag_error_est_failed ;
    } 


  // ------------ mass/age fit ---------------

  coord[0][1] = starGroup[best_g].mass[best_i] ;  
  coord[0][2] = starGroup[best_g].mass[best_j] ; 
  coord[0][3] = starGroup[best_g].age ;
  coord[0][0] = starGroup[best_g].bestSliceScore ;

  coord[1][1] = starGroup[best_g].mass[best_i+1] ;
  coord[1][2] = starGroup[best_g].mass[best_j] ;
  coord[1][3] = starGroup[best_g].age ;
  coord[1][0] = getScore (best_g, best_i+1, best_j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			  time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
			  combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, plimbDark,
			  colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;

  coord[2][1] = starGroup[best_g].mass[best_i-1] ;
  coord[2][2] = starGroup[best_g].mass[best_j] ;
  coord[2][3] = starGroup[best_g].age ;
  coord[2][0] = getScore (best_g, best_i-1, best_j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			  time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			  combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, plimbDark,
			  colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;

  coord[3][1] = starGroup[best_g].mass[best_i] ;
  coord[3][2] = starGroup[best_g].mass[best_j+1] ;
  coord[3][3] = starGroup[best_g].age ;
  coord[3][0] = getScore (best_g, best_i, best_j+1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			  time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			  combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, plimbDark,
			  colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;

  coord[4][1] = starGroup[best_g].mass[best_i] ;
  coord[4][2] = starGroup[best_g].mass[best_j-1] ;
  coord[4][3] = starGroup[best_g].age ;
  coord[4][0] = getScore (best_g, best_i, best_j-1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			  time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			  combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, plimbDark,
			  colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;

  // Note that this part is cheating a bit:
  coord[5][1] = starGroup[best_g].mass[best_i] ;
  coord[5][2] = starGroup[best_g].mass[best_j] ;
  coord[5][3] = starGroup[best_g+1].age ;
  coord[5][0] = starGroup[best_g+1].bestSliceScore ;
  
  coord[6][1] = starGroup[best_g].mass[best_i] ;
  coord[6][2] = starGroup[best_g].mass[best_j] ;
  coord[6][3] = starGroup[best_g-1].age ;
  coord[6][0] = starGroup[best_g-1].bestSliceScore ;


  for (i = 0 ; i < MASS_PARABOLA_FIT_NUM ; i++) 
    if (coord[i][0] >= MAX_VALID_SCORE)
      {
	sprintf (strErr, "Bad neighboring position (%d  %f)", i, coord[i][0]) ;
	printError (errFilename, 0, deb_filename, "parabolaFit() (mass/age)", strErr) ;

	goto mass_and_mag_with_coord_error_est_failed ;
      }
  
  
  for (i = 1 ; i < MASS_PARABOLA_FIT_NUM ; i++)    // shouldn't happen
    if (coord[i][0] < coord[0][0])
      {
	sprintf (strErr, "Best point is not the minimum (%d)", i) ;
	printError (errFilename, 0, deb_filename, "parabolaFit() (mass/age)", strErr) ;
      }


  // ---------- calc parabola fit (mass1-mass2-age) ----------

  for (i = 0 ; i < MASS_PARABOLA_FIT_NUM ; i++)
    {
      fitMatrix[i][0] = sqr(coord[i][1]) ;
      fitMatrix[i][1] = coord[i][1] ;
      fitMatrix[i][2] = sqr(coord[i][2]) ;
      fitMatrix[i][3] = coord[i][2] ;
      fitMatrix[i][4] = sqr(coord[i][3]) ;
      fitMatrix[i][5] = coord[i][3] ;
      fitMatrix[i][6] = 1.0 ;
    }

  if (gaussj(fitMatrix, MASS_PARABOLA_FIT_NUM)) 
    {
      printError (errFilename, 0, deb_filename, "gaussj() (mass/age)", "Degenerate (flat) paraboloid, or not enough memory.") ;
      
      goto mass_and_mag_with_coord_error_est_failed ;
    }
  
  for (i = 0 ; i < MASS_PARABOLA_FIT_NUM ; i++)
    {
      vec[i] = (fitMatrix[i][0] * coord[0][0]) ;

      for (j = 1 ; j < MASS_PARABOLA_FIT_NUM ; j++)
	vec[i] += (fitMatrix[i][j] * coord[j][0]) ;
    }

  if ((vec[0] <= 0.0) || (vec[2] <= 0.0) || (vec[4] <= 0.0))     // shouldn't happen
    {
      printError (errFilename, 0, deb_filename, "parabolaFit() (mass/age)", "One of the axes is convex (has a maximum instead of a minimum).") ;
      
      goto mass_and_mag_with_coord_error_est_failed ;
    }

  minimumOfMassAge = vec[6] - (0.25 * ((vec[1]*vec[1]/vec[0]) + (vec[3]*vec[3]/vec[2]) + (vec[5]*vec[5]/vec[4]))) ;

  if (minimumOfMassAge < 0.0)
    printError (errFilename, 0, deb_filename, "parabolaFit() (mass/age)", "Fitted parabola has a negative minimum score.") ;

  mass1 = -0.5 * vec[1] / vec[0] ;
  mass1Err = sqrt(1.0 / vec[0]) ;
  mass2 = -0.5 * vec[3] / vec[2] ;
  mass2Err = sqrt(1.0 / vec[2]) ;
  age = -0.5 * vec[5] / vec[4] ;
  ageErr = sqrt(1.0 / vec[4]) ;


  if ((mass1 < starGroup[best_g].mass[best_i-1]) ||
      (mass1 > starGroup[best_g].mass[best_i+1]) ||
      (mass2 < starGroup[best_g].mass[best_j-1]) ||
      (mass2 > starGroup[best_g].mass[best_j+1]) ||
      (age < starGroup[best_g-1].age) ||
      (age > starGroup[best_g+1].age))            // shouldn't happen
    printError (errFilename, 0, deb_filename, "parabolaFit() (mass/age)", "Parabola minimum is out of bound.") ;
    
  //------------- combined magnitude fit --------------
  
  coord[0][1] = combineMag_2lc(starGroup[best_g].Mag[indexLCcolor_fA][best_i],
			       starGroup[best_g].Mag[indexLCcolor_fA][best_j],
			       starGroup[best_g].Mag[indexLCcolor_fB][best_i],
			       starGroup[best_g].Mag[indexLCcolor_fB][best_j]) ;
  coord[0][0] = starGroup[best_g].bestSliceScore ;
  
  coord[1][1] = combineMag_2lc(starGroup[best_g].Mag[indexLCcolor_fA][best_i+1],
			       starGroup[best_g].Mag[indexLCcolor_fA][best_j+1],
			       starGroup[best_g].Mag[indexLCcolor_fB][best_i+1],
			       starGroup[best_g].Mag[indexLCcolor_fB][best_j+1]) ;
  coord[1][0] = getScore(best_g, best_i+1, best_j+1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			 time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
			 combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, plimbDark,
			 colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
  
  coord[2][1] = combineMag_2lc(starGroup[best_g].Mag[indexLCcolor_fA][best_i-1],
			       starGroup[best_g].Mag[indexLCcolor_fA][best_j-1],
			       starGroup[best_g].Mag[indexLCcolor_fB][best_i-1],
			       starGroup[best_g].Mag[indexLCcolor_fB][best_j-1]) ;
  coord[2][0] = getScore(best_g, best_i-1, best_j-1, indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
			 time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
			 combObsMag_fA, combObsMag_fB, e, t0, omega, LCweighting, plimbDark,
			 colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;

  // for (i = 0 ; i < MAGNITUDE_PARABOLA_FIT_NUM ; i++)
  //  printf ("*** %d  %f  %f\n", i, coord[i][0], coord[i][1]) ;
  

  for (i = 0 ; i < MAGNITUDE_PARABOLA_FIT_NUM ; i++) 
    if (coord[i][0] >= MAX_VALID_SCORE)
      {
	sprintf (strErr, "Bad neighboring position (%d  %f)", i, coord[i][0]) ;
	printError (errFilename, 0, deb_filename, "parabolaFit() (mag)", strErr) ;

	goto mag_error_est_failed ;
      }
  
  for (i = 1 ; i < MAGNITUDE_PARABOLA_FIT_NUM ; i++)
    if (coord[i][0] < coord[0][0])
      {
	sprintf (strErr, "The given best point is not the minimum (%f %f %f %d)",
		 starGroup[best_g].mass[best_i], starGroup[best_g].mass[best_j], starGroup[best_g].age, i) ;
	printError (errFilename, 0, deb_filename, "parabolaFit() (mag)", strErr) ;
      }

  // ---------- calc parabola fit (comb mag) ----------


  for (i = 0 ; i < MAGNITUDE_PARABOLA_FIT_NUM ; i++)
    {
      fitMatrix[i][0] = sqr(coord[i][1]) ;
      fitMatrix[i][1] = coord[i][1] ;
      fitMatrix[i][2] = 1.0 ;
    }

  if (gaussj(fitMatrix, MAGNITUDE_PARABOLA_FIT_NUM)) 
    {
      printError (errFilename, 0, deb_filename, "gaussj() (mag)", "Degenerate (flat) paraboloid, or not enough memory.") ;

      goto mag_error_est_failed ;
    }
  
  for (i = 0 ; i < MAGNITUDE_PARABOLA_FIT_NUM ; i++)
    {
      vec[i] = (fitMatrix[i][0] * coord[0][0]) ;
      
      for (j = 1 ; j < MAGNITUDE_PARABOLA_FIT_NUM ; j++)
	vec[i] += (fitMatrix[i][j] * coord[j][0]) ;
    }
  
  if (vec[0] <= 0.0)
    {
      printError (errFilename, 0, deb_filename, "parabolaFit() (mag)", "One of the axes is convex (has a maximum instead of a minimum).") ;

      goto mag_error_est_failed ;
    }

  minimumOfMag = vec[2] - (0.25*vec[1]*vec[1]/vec[0]) ;

  if (minimumOfMag < 0.0)
    printError (errFilename, 0, deb_filename, "parabolaFit() (mag)", "Fitted parabola has a negative minimum score.") ;

  combAbsMag = -0.5 * vec[1] / vec[0] ;
  combAbsMagErr = sqrt(1.0 / vec[0]) ;

  //------------------------------ all OK -----------------------------------

  fprintf (foutResults, "%s  %f %f  %f %f  %f %f  %f   %f %f  %f   %f %f   %s %s  %.3f %.3f %.3f\n",
	   deb_filename, mass1, mass1Err, mass2, mass2Err, age, ageErr, minimumOfMassAge,
	   combAbsMag, combAbsMagErr, minimumOfMag, 
	   combObsMag_fAB, 10.0 * pow(10.0, 0.2 * (combObsMag_fAB - combAbsMag)),
	   strRA, strDEC, pmRA, pmDEC, ra_dec_err) ;


  fflush (foutResults) ;

  return (0) ;

  //-------------------------- missing mag info -----------------------------
  
 mag_error_est_failed:    // need to set: mass1, mass1Err, mass2, mass2Err, age, ageErr, minimumOfMassAge
  
  combAbsMag = combineMag_2lc(starGroup[best_g].Mag[indexLCcolor_fA][best_i],
			      starGroup[best_g].Mag[indexLCcolor_fA][best_j],
			      starGroup[best_g].Mag[indexLCcolor_fB][best_i],
			      starGroup[best_g].Mag[indexLCcolor_fB][best_j]) ;
  
  fprintf (foutResults, "%s  %f %f  %f %f  %f %f  %f   %f %f  %f   %f %f   %s %s  %.3f %.3f %.3f\n",
	   deb_filename, mass1, mass1Err, mass2, mass2Err, age, ageErr, minimumOfMassAge,
	   combAbsMag, UNCERTAINTY_UNKNOWN, starGroup[best_g].bestSliceScore, combObsMag_fAB,
	   10.0 * pow(10.0, 0.2 * (combObsMag_fAB - combAbsMag)),
	   strRA, strDEC, pmRA, pmDEC, ra_dec_err) ;
  
  fflush (foutResults) ;
  
  return (2) ;

  //-------------------- missing mass and mag info --------------------------

 mass_and_mag_with_coord_error_est_failed:  // need to set:  coord[][]

  if ((coord[1][0] >= MAX_VALID_SCORE) && (coord[2][0] < MAX_VALID_SCORE))
    mass1Err = UNCERTAINTY_UPPER_LIMIT ;
  else
    if ((coord[2][0] >= MAX_VALID_SCORE) && (coord[1][0] < MAX_VALID_SCORE))
      mass1Err = UNCERTAINTY_LOWER_LIMIT ;
    else
      mass1Err = UNCERTAINTY_UNKNOWN ;
  
  if ((coord[3][0] >= MAX_VALID_SCORE) && (coord[4][0] < MAX_VALID_SCORE))
    mass2Err = UNCERTAINTY_UPPER_LIMIT ;
  else
    if ((coord[4][0] >= MAX_VALID_SCORE) && (coord[3][0] < MAX_VALID_SCORE))
      mass2Err = UNCERTAINTY_LOWER_LIMIT ;
    else
      mass2Err = UNCERTAINTY_UNKNOWN ;
  
  if ((coord[5][0] >= MAX_VALID_SCORE) && (coord[6][0] < MAX_VALID_SCORE))
    ageErr = UNCERTAINTY_UPPER_LIMIT ;
  else 
    if ((coord[6][0] >= MAX_VALID_SCORE) && (coord[5][0] < MAX_VALID_SCORE))
      ageErr = UNCERTAINTY_LOWER_LIMIT ;
    else
      ageErr = UNCERTAINTY_UNKNOWN ;


 mass_and_mag_error_est_failed:      // need to set:  mass1Err, mass2Err, ageErr
  
  combAbsMag = combineMag_2lc(starGroup[best_g].Mag[indexLCcolor_fA][best_i],
			      starGroup[best_g].Mag[indexLCcolor_fA][best_j],
			      starGroup[best_g].Mag[indexLCcolor_fB][best_i],
			      starGroup[best_g].Mag[indexLCcolor_fB][best_j]) ;
  
  fprintf (foutResults, "%s  %f %f  %f %f  %f %f  %f   %f %f  %f   %f %f   %s %s  %.3f %.3f %.3f\n",
	   deb_filename, starGroup[best_g].mass[best_i], mass1Err, 
	   starGroup[best_g].mass[best_j], mass2Err, starGroup[best_g].age, ageErr, 
	   starGroup[best_g].bestSliceScore, combAbsMag, UNCERTAINTY_UNKNOWN,
	   starGroup[best_g].bestSliceScore, combObsMag_fAB,
	   10.0 * pow(10.0, 0.2 * (combObsMag_fAB - combAbsMag)),
	   strRA, strDEC, pmRA, pmDEC, ra_dec_err) ;
  fflush (foutResults) ;
  
  return (3) ;
}



int main (int argc, char **argv)
{
  int g, i, j, n, numGroups, size_fA, size_fB, numColors, indexLCcolor_fA, indexLCcolor_fB ;
  int isInvalidColorFormat, numElements, numElements_orig ;
  char deb_filename_fA[MAX_STRLEN], deb_filename_fB[MAX_STRLEN] ;
  char strErr[MAX_STRLEN], deb_filenamePeriod_fA[MAX_STRLEN], deb_filenamePeriod_fB[MAX_STRLEN] ;
  double score, splineScore, splineScore_fA, splineScore_fB, combObsMag_fA, combObsMag_fB ;
  float LCweighting, pmRA, pmDEC, ra_dec_err ;
  float *time_fA, *amp_fA, *err_fA, *time_fB, *amp_fB, *err_fB ;
  double period_fA, ecc_fA, deb_mag1_fA, deb_mag2_fA, t0_fA, omega_fA, chi2_fA ;
  double period_fB, ecc_fB, deb_mag1_fB, deb_mag2_fB, t0_fB, omega_fB, chi2_fB ;
  double period, ecc, t0, omega ;
  double bestBasinScore, bestSliceScore, bestScore ;
  char colorLC_fA, colorLC_fB, colorFilterPos, colorFilterNeg ;
  char strRA[MAX_STRLEN], strDEC[MAX_STRLEN] ;
  double *contourMatrix ;
  int *subsetArr, search_i, search_j,  bestSlice_i, bestSlice_j, best_g, best_i, best_j ;
  unsigned long countScoreGood, countScoreMax, countScoreNonDetached, countScoreUnobservable ;
  unsigned long countScoreOutBoundLimbDark ;
  FILE *finDebil, *foutResults ; 
  struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS] ;
  struct ColorInfoType *colorInfoArr ;
  struct LimbDarkType limbDark ;
#ifdef CONTOUR_OUTPUT_SUFFIX
  char filenameContour[MAX_STRLEN], filenameIDL[MAX_STRLEN] ;
  FILE *foutContourData, *foutIDL ;
#endif
#ifdef WRITE_MODEL_FIT_LC_fA
  char deb_filenamePeriodFit[MAX_STRLEN] ;
  FILE *foutFit_fA, *foutFit_fB, *foutDat_fA, *foutDat_fB ; 
#endif

  if ((argc != 6) && (argc != 7) && (argc != 8))
    {
      printf ("\n Usage:  %s  <2LC DEBiL DB (31-param) /w color>  <limb-darkening DB>  <isochrone DB>  <output result file>  <output error file>  [num elements (default=%d)]   [LC weighting* (default=%.2f)]\n", argv[0], DEFAULT_NUM_RESOLUTION_ELEMENTS, DEFAULT_LC_WEIGHTING) ;
      printf ("  {version %.1f}\n", VERSION) ;
      printf ("*A negative sign in the LC weighting value, multiplies it by the number of data points.\n") ;
      return (1) ;
    }

  if (!strcmp(argv[1], argv[2]) ||
      !strcmp(argv[1], argv[3]) ||
      !strcmp(argv[1], argv[4]) ||
      !strcmp(argv[1], argv[5]) ||
      !strcmp(argv[2], argv[3]) ||
      !strcmp(argv[2], argv[4]) ||
      !strcmp(argv[2], argv[5]) ||
      !strcmp(argv[3], argv[4]) ||
      !strcmp(argv[3], argv[5]) ||
      !strcmp(argv[4], argv[5]))
    {
      printf ("ERROR: All the input/output filenames must be different.\n") ;
      return (2) ;
    }

  if (argc >= 7)
    {
      numElements_orig = atoi (argv[6]) ;

      if (numElements_orig <= 1)
	{
	  printf ("ERROR: Invalid number of elements (%d <= 1).\n", numElements_orig) ;
	  return (3) ;
	}
    }
  else numElements_orig = DEFAULT_NUM_RESOLUTION_ELEMENTS ;


  if (argc >= 8)
    {
      LCweighting = atof (argv[7]) ;

      if (LCweighting > 0.0)
	printf ("The LC weighting = %f.\n", LCweighting) ;
      else if (LCweighting < 0.0)
	printf ("The LC weighting = %fx the number of non-outlier data points.\n", -LCweighting) ;
      else // (LCweighting == 0.0)
	printf ("Warning: No LC information will be used, only color.\n") ;
    }  
  else LCweighting = DEFAULT_LC_WEIGHTING ;

  finDebil = fopen (argv[1], "rt") ;
  if (!finDebil)
    {
      printf ("ERROR: Couldn't open the DEBiL input file ('%s').\n", argv[1]) ;
      return (4) ;
    }  

  foutResults = fopen (argv[4], "at") ;
  if (!foutResults)
    {
      printf ("ERROR: Couldn't open the results output file ('%s').\n", argv[4]) ;
      fclose (finDebil) ;
      return (5) ;
    } 
 
  //-------------------

  n = loadLimbDark (argv[2], &limbDark) ;
  if (n)
    {
      fclose (finDebil) ;
      fclose (foutResults) ;
      return (6) ;
    }

  numGroups = loadIsochrones (argv[3], starGroup) ;
  if (numGroups <= 0)
    {
      fclose (finDebil) ;
      fclose (foutResults) ;
      return (7) ;
    }

#ifdef CONTOUR_OUTPUT_SUFFIX
  sprintf (filenameIDL, "%s%s", argv[4], IDL_SCRIPT_SUFFIX) ;
  foutIDL = fopen (filenameIDL, "wt") ;
  if (!foutIDL)
    {
      printf ("ERROR: Couldn't create IDL output file '%s'.", filenameIDL) ;
      fclose (finDebil) ;
      fclose (foutResults) ;
      return (8) ;
    }
  
  fprintf (foutIDL, "set_plot, 'ps'\n\n") ; 
#endif

  //-----------------------------------------------------
 
  while ((8 == fscanf (finDebil, "%s %lf %lf %*f %*f %*f %*f %*f %lf %*f %lf %*f %*f %*f %lf %*f %lf %*f %*d %*d %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f\n", deb_filename_fA, &period_fA, &ecc_fA, &deb_mag1_fA, &deb_mag2_fA, &t0_fA, &omega_fA, &chi2_fA)) &&
	 (8 == fscanf (finDebil, "%s %lf %lf %*f %*f %*f %*f %*f %lf %*f %lf %*f %*f %*f %lf %*f %lf %*f %*d %*d %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f\n", deb_filename_fB, &period_fB, &ecc_fB, &deb_mag1_fB, &deb_mag2_fB, &t0_fB, &omega_fB, &chi2_fB)))
    {
      if (chi2_fA <= chi2_fB)
	{
	  printf ("Using the orbital parameters of line 1\n") ;
	  period = period_fA ;      
	  ecc = ecc_fA ;
	  t0 = t0_fA ;
	  omega = omega_fA ;
	}
      else  //  (chi2_fA > chi2_fB)
	{
	  printf ("Using the orbital parameters of line 2\n") ;
	  period = period_fB ;      
	  ecc = ecc_fB ;
	  t0 = t0_fB ;
	  omega = omega_fB ;
	}

      combObsMag_fA = combineMag(deb_mag1_fA, deb_mag2_fA) ;
      combObsMag_fB = combineMag(deb_mag1_fB, deb_mag2_fB) ;

      printf ("LCs= '%s' + '%s'  ; P=%.3f days.\n", deb_filename_fA, deb_filename_fB, period) ;
      
      for (j = 0, i = 0 ; i < strlen(deb_filename_fA) ; i++)
	if (deb_filename_fA[i] == '/')
	  j = i+1 ;

      sprintf (deb_filenamePeriod_fA, "%s_P%.3f", &deb_filename_fA[j], period) ;
      sprintf (deb_filenamePeriod_fB, "%s_P%.3f", &deb_filename_fB[j], period) ;

      // ---------- read color info: ------------
      if (8 != fscanf (finDebil, "%c %c %d  %s  %s  %f %f %f\n",
		       &colorLC_fA, &colorLC_fB, &numColors, strRA, strDEC, &pmRA, &pmDEC, &ra_dec_err))
	{
	  printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Invalid LC color header.") ;
	  break ;
	}

      indexLCcolor_fA = colorConvert (colorLC_fA) ;    
      indexLCcolor_fB = colorConvert (colorLC_fB) ;   
     
      if ((indexLCcolor_fA < 0) || (indexLCcolor_fB < 0))
	{
	  printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Invalid LC color.") ;
	  for (i = 0 ; i < numColors ; i++)
	    nextLine (finDebil) ;
	  continue ;
	}

      if (numColors > 0)
	{
	  colorInfoArr = (struct ColorInfoType*)malloc (numColors * sizeof(struct ColorInfoType)) ;
	  if (!colorInfoArr)
	    {
	      printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Couldn't allocate colorInfoArr.") ;
	      for (i = 0 ; i < numColors ; i++)
		nextLine (finDebil) ;
	      continue ;
	    }

	  isInvalidColorFormat = 0 ;

	  for (i = 0 ; (i < numColors) && !feof(finDebil) ; i++)
	    {
	      n = fscanf (finDebil, "%c-%c %lf %lf\n", &colorFilterPos, &colorFilterNeg, 
			  &(colorInfoArr[i].colorVal), &(colorInfoArr[i].colorErr)) ;
	      
	      printf ("(%c-%c) = %f +/- %f\n",
		      colorFilterPos, colorFilterNeg, colorInfoArr[i].colorVal, colorInfoArr[i].colorErr) ;

	      if ((n != 4) || (colorInfoArr[i].colorErr <= 0.0))
		isInvalidColorFormat = 1 ;    // you don't want to break and lose the file sync.

	      colorInfoArr[i].filterPos = colorConvert (colorFilterPos) ;
	      colorInfoArr[i].filterNeg = colorConvert (colorFilterNeg) ;

	      if ((colorInfoArr[i].filterPos < 0) || (colorInfoArr[i].filterNeg < 0))
		{
		  sprintf (strErr, "Invalid filter color (%c-%c).", colorFilterPos, colorFilterNeg) ;
		  printError (argv[5], 0, deb_filenamePeriod_fA, "main()", strErr) ;
		  i-- ;
		  numColors-- ;
		}
	    }

	  if (isInvalidColorFormat)
	    {
	      printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Invalid filter color format.") ;
	      if (colorInfoArr) free (colorInfoArr) ;
	      continue ;
	    }
	}
      else if (numColors < 0)
	{
	  printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Invalid number of colors.") ;
	  continue ;
	}
      else  // (numColors == 0)
	{
	  if (LCweighting == 0.0) continue ;
	  colorInfoArr = 0 ;
	}

      //--------------------------------------------
      
      size_fA = loadLC (deb_filename_fA, period, argv[5], &time_fA, &amp_fA, &err_fA) ; 
      if (size_fA < MIN_LIGHT_CURVE_DATA_POINTS)   // if loadLC() fails then (size < 0)
	{
	  if (size_fA >= 0)
	    printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Too few light curve data points.") ;

	  if (colorInfoArr) free (colorInfoArr) ;
	  continue ;
	}

      size_fB = loadLC (deb_filename_fB, period, argv[5], &time_fB, &amp_fB, &err_fB) ; 
      if (size_fB < MIN_LIGHT_CURVE_DATA_POINTS)   // if loadLC() fails then (size < 0)
	{
	  if (size_fB >= 0)
	    printError (argv[5], 1, deb_filenamePeriod_fB, "main()", "Too few light curve data points.") ;

	  if (colorInfoArr) free (colorInfoArr) ;
	  free(time_fA) ;
	  free(amp_fA) ;
	  free(err_fA) ;
	  continue ;
	}

      if (sortSamples(time_fA, amp_fA, err_fA, size_fA))
	{
	  printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "sortSamples() malloc failure.") ;
	  free (colorInfoArr) ;
	  free(time_fA) ;
	  free(amp_fA) ;
	  free(err_fA) ;
	  free(time_fB) ;
	  free(amp_fB) ;
	  free(err_fB) ;
	  continue ;	
	}

      if (sortSamples(time_fB, amp_fB, err_fB, size_fB))
	{
	  printError (argv[5], 1, deb_filenamePeriod_fB, "main()", "sortSamples() malloc failure.") ;
	  free (colorInfoArr) ;
	  free(time_fA) ;
	  free(amp_fA) ;
	  free(err_fA) ;
	  free(time_fB) ;
	  free(amp_fB) ;
	  free(err_fB) ;
	  continue ;	
	}
     
      i = size_fA ;
      splineScore_fA = ridOutliers(time_fA, amp_fA, err_fA, &size_fA, OUTLIER_VARIANCE_LIMIT) ;
      
      if (i != size_fA)
	{
	  if ((i - size_fA) == 1)
	    sprintf (strErr, "removed 1 data point.") ;
	  else
	    sprintf (strErr, "removed %d data points.", i - size_fA) ;
	  
	  printError (argv[5], 0, deb_filenamePeriod_fA, "ridOutliers()", strErr) ;
	}

      i = size_fB ;
      splineScore_fB = ridOutliers(time_fB, amp_fB, err_fB, &size_fB, OUTLIER_VARIANCE_LIMIT) ;
      
      if (i != size_fB)
	{
	  if ((i - size_fB) == 1)
	    sprintf (strErr, "removed 1 data point.") ;
	  else
	    sprintf (strErr, "removed %d data points.", i - size_fB) ;
	  
	  printError (argv[5], 0, deb_filenamePeriod_fA, "ridOutliers()", strErr) ;
	}

      splineScore = ((splineScore_fA * size_fA) + (splineScore_fB * size_fB)) / (size_fA + size_fB) ;

      // Note: the Spline score is very approximate. Too low when the phased light curve is flat,
      // and too high when it has high frequency features.
      // The errors are assumed to be underestimated by the square-root of this value.
#ifdef DEBUG
      printf ("Spline score = %f\n", splineScore) ;
#endif

#ifdef DO_SPLINE_SCORE_CORRECTION
      if (splineScore <= 0.0)
	{
	  printError (argv[5], 1, deb_filenamePeriod_fA, "ridOutliers()", "invalid spline score.") ;
	  free (colorInfoArr) ;
	  free(time) ;
	  free(amp) ;
	  free(err) ;
	  continue ;	
	}
      else
	{
	  sprintf (strErr, "the errors were normalized by a factor %f", sqrt(splineScore)) ;
	  printError (argv[5], 0, deb_filenamePeriod_fA, "main()", strErr) ;
	}
#else
      splineScore = 1.0 ;
#endif

      //--------------------------------------------------------

      countScoreGood = 0 ;
      countScoreMax = 0 ; 
      countScoreNonDetached = 0 ;
      countScoreUnobservable = 0 ; 
      countScoreOutBoundLimbDark = 0 ;
      
      bestScore = MAXFLOAT ;
      best_g = 0 ;
      best_i = 0 ;
      best_j = 0 ;

      // g = isochrone group number (each group hase a given age/metalicity)
      for (g = 0 ; g < numGroups ; g++)
	{
	  if (starGroup[g].arrLen <= 1)  // needed to protect the numElement calculation
	    {
	      sprintf (strErr, "Slice %d (age=%f) only has %d mass levels- skipping.", g, starGroup[g].age, starGroup[g].arrLen) ;
	      printError (argv[5], 1, deb_filenamePeriod_fA, "main()", strErr) ;
	      continue ;
	    }

	  contourMatrix = (double*)malloc (starGroup[g].arrLen * starGroup[g].arrLen * sizeof(double)) ;
	  subsetArr = (int*)malloc (starGroup[g].arrLen * sizeof(int)) ;

	  if (!contourMatrix || !subsetArr)
	    {
	      if (contourMatrix) free (contourMatrix) ;
	      if (subsetArr) free (subsetArr) ;
	      printError (argv[5], 1, deb_filenamePeriod_fA, "main()", "Couldn't allocate enough memory.") ;
	      continue ;
	    }

	  for (i = 0 ; i < (starGroup[g].arrLen * starGroup[g].arrLen) ; i++)
	    contourMatrix[i] = RESULT_SKIPPED ;

	  //----------------------

	  bestSliceScore = MAXFLOAT ;
	  bestSlice_i = 0 ;
	  bestSlice_j = 0 ;
	  
	  if (numElements_orig > starGroup[g].arrLen)
	    {
	      numElements = starGroup[g].arrLen ;
	      printError (argv[5], 0, deb_filenamePeriod_fA, "main()", "Decreased the number of elements to the size of the mass array.") ;
	    }
	  else numElements = numElements_orig ;

	  calcEqualSpaced (subsetArr, numElements, starGroup[g].mass, starGroup[g].arrLen) ;

	  // ----------- grid scan (not very efficient, but very robust) ----------------------

	  // i = isochrone index- column  (mass of star #1)
	  // j = isochrone index- row     (mass of star #2)
	  for (i = 0 ; i < numElements ; i++)
	    for (j = 0 ; j < numElements ; j++)
	      {		
		score = getScore (g, subsetArr[i], subsetArr[j], indexLCcolor_fA, indexLCcolor_fB, period, splineScore,
				  time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
				  combObsMag_fA, combObsMag_fB, ecc, t0, omega, LCweighting, &limbDark,
				  colorInfoArr, numColors, starGroup, 0, 0, 0, 0, 0, 0) ;
		
		contourMatrix[(starGroup[g].arrLen * subsetArr[i]) + subsetArr[j]] = score ;
	       
		// printf ("%f %f %f\n",  starGroup[g].mass[subsetArr[i]], starGroup[g].mass[subsetArr[j]], score) ;

		if (score < bestSliceScore)
		  {
		    bestSliceScore = score ;
		    bestSlice_i = subsetArr[i] ;
		    bestSlice_j = subsetArr[j] ;
		  }
		
		switch ((unsigned long)score)
		  {
		  case (unsigned long)MAX_VALID_SCORE:               countScoreMax++ ;              break ;
		  case (unsigned long)RESULT_NONDETACHED:            countScoreNonDetached++ ;      break ;
		  case (unsigned long)RESULT_UNOBSERVABLE_DIP:       countScoreUnobservable++ ;     break ;
		  case (unsigned long)RESULT_OUT_OF_BOUND_LIMB_DARK: countScoreOutBoundLimbDark++ ; break ;
		  default:                                           countScoreGood++ ;
		  }
	      }
	  
	  // --------- do a greedy search in each of the basins (significant local minima) -------
	  for (i = 0 ; i < numElements ; i++)
	    for (j = 0 ; j < numElements ; j++)
	      if (isBasin (i, j, contourMatrix, starGroup[g].arrLen, subsetArr, numElements))
		{		  
		  search_i = subsetArr[i] ;
		  search_j = subsetArr[j] ;
		  bestBasinScore = greedyFindMin(g, &search_i, &search_j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
						 time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
						 combObsMag_fA, combObsMag_fB, ecc, t0, omega, LCweighting,
						 &limbDark, colorInfoArr, numColors, starGroup, contourMatrix,
						 starGroup[g].arrLen) ;
		  
		  if (bestBasinScore < bestSliceScore)
		    {
		      bestSlice_i = search_i ;
		      bestSlice_j = search_j ;
		      bestSliceScore = bestBasinScore ;
		    }
		}

	  //-------------------------------------------------------------------
	  
	  starGroup[g].bestSliceScore = bestSliceScore ;
	  
	  if (bestSliceScore < bestScore)
	    {
	      bestScore = bestSliceScore ;
	      best_g = g ;
	      best_i = bestSlice_i ;
	      best_j = bestSlice_j ;
	    }
	  
#ifdef DEBUG
	  printf ("Age= %f  (%dx%d)  slice_best_score= %f   masses= %f %f  (%d,%d;%d)\n",
		  starGroup[g].age, numElements, numElements, bestSliceScore, 
		  starGroup[g].mass[bestSlice_i], starGroup[g].mass[bestSlice_j],
		  bestSlice_i, bestSlice_j, starGroup[g].arrLen) ;
#endif
	  
#ifdef CONTOUR_OUTPUT_SUFFIX
	  if (bestSliceScore < MAX_VALID_SCORE)
	    {	      
	      sprintf (filenameContour, "%s_age%.3f%s", deb_filenamePeriod_fA, starGroup[g].age, CONTOUR_OUTPUT_SUFFIX) ;
	      foutContourData = fopen(filenameContour, "wt") ;
	      if (!foutContourData)
		{
		  printError (argv[5], 0, deb_filenamePeriod_fA, "main()", "Couldn't open the output contour data file.") ;
		  free (contourMatrix) ;
		  free (subsetArr) ;
		  continue ;
		}

	      if (contourHeader(deb_filenamePeriod_fA, starGroup[g].age, foutIDL, filenameContour,
				contourMatrix, starGroup[g].arrLen, subsetArr, numElements))
		{
		  printError (argv[5], 0, deb_filenamePeriod_fA, "contourHeader()", "Failed to create contour plot.") ;
		  free (contourMatrix) ;
		  free (subsetArr) ;
		  continue ;
		}

	      fprintf (foutContourData,"%d\n", numElements) ;
	  
	      for (i = 0 ; i < numElements ; i++)
		fprintf (foutContourData,"%f ", starGroup[g].mass[subsetArr[i]]) ;
	      
	      fprintf (foutContourData,"\n") ;
	    
	      for (j = 0 ; j < numElements ; j++)
		{
		  for (i = 0 ; i < numElements ; i++)
		    fprintf (foutContourData,"%f ", contourMatrix[(starGroup[g].arrLen * subsetArr[i]) + subsetArr[j]]) ;
	      
		  fprintf (foutContourData,"\n") ;
		}

	      fclose(foutContourData) ;
	    }
	  else
	    {
	      sprintf (strErr, "No good pairings at age=%.3fGyr, didn't create a contour.", starGroup[g].age) ;
	      printError (argv[5], 0, deb_filenamePeriod_fA, "main()", strErr) ;
	    }
#endif

	  free (contourMatrix) ;
	  free (subsetArr) ;
	}  // end g loop (age)


      if (bestScore < MAX_VALID_SCORE)
	{
	  parabolaFit (best_g, best_i, best_j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
		       time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB, 
		       combObsMag_fA, combObsMag_fB, ecc, t0, omega, LCweighting, strRA, strDEC, pmRA, pmDEC, ra_dec_err,
		       &limbDark, colorInfoArr, numColors,  starGroup, numGroups, argv[5],
		       deb_filenamePeriod_fA, foutResults) ;
	}
      else
	{
#ifdef SET_ZERO_ECCENTRICITY
	  sprintf (strErr, "No good pairings at ALL ages (might be dwarfs)  [good=%lu max=%lu nonDet=%lu unObs=%lu outBoundLD=%lu]", countScoreGood, countScoreMax, countScoreNonDetached, countScoreUnobservable, countScoreOutBoundLimbDark) ;
#else
	  sprintf (strErr, "No good pairings at ALL ages (might be dwarfs or have an over-eccentric fit: e=%f)  [good=%lu max=%lu nonDet=%lu unObs=%lu outBoundLD=%lu]", ecc, countScoreGood, countScoreMax, countScoreNonDetached, countScoreUnobservable, countScoreOutBoundLimbDark) ;
#endif
	  
	  printError (argv[5], 1, deb_filenamePeriod_fA, "main()", strErr) ;
	}

#ifdef WRITE_MODEL_FIT_LC_fA
      if (LCweighting != 0.0)
	{
	  sprintf (deb_filenamePeriodFit, "%s%s", deb_filenamePeriod_fA, WRITE_MODEL_FIT_LC_fA) ;
	  foutFit_fA = fopen (deb_filenamePeriodFit, "wt") ;
	  sprintf (deb_filenamePeriodFit, "%s%s", deb_filenamePeriod_fA, WRITE_MODEL_DAT_LC_fA) ;
	  foutDat_fA = fopen (deb_filenamePeriodFit, "wt") ;

	  sprintf (deb_filenamePeriodFit, "%s%s", deb_filenamePeriod_fB, WRITE_MODEL_FIT_LC_fB) ;
	  foutFit_fB = fopen (deb_filenamePeriodFit, "wt") ;
	  sprintf (deb_filenamePeriodFit, "%s%s", deb_filenamePeriod_fB, WRITE_MODEL_DAT_LC_fB) ;
	  foutDat_fB = fopen (deb_filenamePeriodFit, "wt") ;
	  
	  if (foutFit_fA && foutDat_fA && foutFit_fB && foutDat_fB)
	    {
	      getScore (best_g, best_i, best_j, indexLCcolor_fA, indexLCcolor_fB, period, splineScore, 
			time_fA, amp_fA, err_fA, size_fA, time_fB, amp_fB, err_fB, size_fB,
			combObsMag_fA, combObsMag_fB, ecc, t0, omega, LCweighting, &limbDark, colorInfoArr, numColors, starGroup, 
			foutFit_fA, foutDat_fA, foutFit_fB, foutDat_fB, argv[5], deb_filenamePeriod_fA) ;
	      fclose (foutFit_fA) ;
	      fclose (foutDat_fA) ;
	      fclose (foutFit_fB) ;
	      fclose (foutDat_fB) ;
	    }
	  else
	    printError (argv[5], 0, deb_filenamePeriod_fA, "main()", "Couldn't write the model fit file.") ;
	}
      else // LCweighting == 0.0
	printError (argv[5], 0, deb_filenamePeriod_fA, "main()", "Since the light curve was ignored, no fit file was written.") ;
#endif

      //------------------------------
            
      if (colorInfoArr)	free (colorInfoArr) ;  // if there was no color info, the array was set to null
      free(time_fA) ;
      free(amp_fA) ;
      free(err_fA) ;
      free(time_fB) ;
      free(amp_fB) ;
      free(err_fB) ;

#ifdef DEBUG
      printf ("counts:  good=%lu  max=%lu  nonDet=%lu  unObs=%lu outBoundLD=%lu\n", 
	      countScoreGood, countScoreMax, countScoreNonDetached, countScoreUnobservable, countScoreOutBoundLimbDark) ;
#endif
    }  // end DEBiL loop
  
  //--------------- end ---------------
  
  for (i = 0 ; i < numGroups ; i++)
    {
      free (starGroup[i].mass) ;
      free (starGroup[i].Teff) ;
      free (starGroup[i].radius) ;

      for (j = 0 ; j < NUM_FILTERS ; j++)
	free (starGroup[i].Mag[j]) ;
    }

#ifdef CONTOUR_OUTPUT_SUFFIX
  fprintf (foutIDL, "end\n") ;
  fclose (foutIDL) ;
#endif

  fclose (finDebil) ;
  fclose (foutResults) ;
  return (0) ;
}
