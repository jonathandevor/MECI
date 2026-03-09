/************************************************************************
 * Simulates eclipsing binary light curves, for random mass-mass-age values, and given noise
 *
 * Written by Jonathan Devor (jdevor@cfa.harvard.edu)
 *    Harvard-Smithsonian Center for Astrophysics 
 *    60 Garden St., Cambridge, MA 02138 USA
 *
 * compile:  gcc meci_sim.c -o meci_sim  -lm -Wall -O4
 *           gcc meci_sim.c -o meci_sim_ecc  -lm -Wall -O4
 * run:      meci_sim ../LimbDarkTable.txt ../SunIsochrone.txt  test0.01.dblc test0.01.log  test0.01_ 0.01 10
 *           meci_sim ../LimbDarkTable.txt ../SunIsochrone.txt  test0.01.dblc test0.01.log  test0.01_ 0.01 10 0.05 H-K 
 *************************************************************************/

#include <math.h>
#include <stdio.h>
#include <values.h>      // for MAXFLOAT
#include <stdlib.h>      // for malloc()
#include <string.h>      // for memset(), memmove(), strcmp(), strlen() and strcat()
#include "debil_lib.c"

#define VERSION 1.2

//#define DO_ECCENTRICITY   // randomizes the value of the orbital eccentricity

#define NOISE_THRESHOLD 1.0  // the number of sigma both eclipses must protrude below the noise level

#define DEBUG

#define MAX_AGE_GYEARS 10.0
#define MIN_AGE_GYEARS 0.2

#define DETACHED_MAX_r1r2 0.8
#define TRUNCATE_RANGE 0.95    // this is meant to protect against evolving stars

#define MAX_NUM_LIMB_DARK_Teff 100 // needed since limb darkening coefficients are read once
#define MAX_NUM_STAR_GROUPS 1000   // needed since the groups are read once, and their number is unknown

#define MAX_STRLEN 512
#define MAX_NUM_COLORS 100

#define LC_FILTER 'R'

// write LC:
#define RAND_SEED 14564327
#define WRITE_SIZE 1000


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
} ;


//--------- math utilities ----------------------


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


// Given their combined magnitude (combMag) and the difference between the magnitudes (Mag2-Mag1 = deltaMag)
// returns the magnitudes of the two stars.
void calcMags (double combMag, double deltaMag, double *mag1_V, double *mag2_V)
{
  *mag1_V = combMag + (2.5 * log10(1.0 + pow(10.0, -0.4 * deltaMag))) ;
  *mag2_V = *mag1_V + deltaMag ;
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



//-------------------------------------------------

// returns values in a random guassian distribuion: mean=0, stddiv=1
double randGauss()
{
  double rho = (double)rand() / ((double)RAND_MAX+1.0) ;
  double sina = sin(2.0 * M_PI * (double)rand() / ((double)RAND_MAX+1.0)) ;

  return (sina * sqrt(-2.0 * log(1.0 - rho))) ;
}


// Writes the observed and model light curve
// Returns the number of points in the eclipse, above the noise level.
// If failed - returns a negative value
int writeLC_sim (FILE *foutFit, double period, double e, double r1, double r2, double mag1, double mag2,
		 double sin_i, double time0, double omega, float magErr,
		 double quadA1, double quadB1, double quadA2, double quadB2, double platFlux,
		 double *peclipse1Flux, double *peclipse2Flux)
{
  int i, countNumInEclipse = 0 ;
  double time, F, p[DIM] ;
  double prevFlux, prevprevFlux ;

  *peclipse1Flux = 0.0 ;
  *peclipse2Flux = 0.0 ;

  setLimbDarkening (quadA1, quadB1, quadA2, quadB2) ;
  
  omega /= 360.0 ;

  p[D_ECC] = e  ;
  p[D_R1] = r1 ;
  p[D_R2] = r2 ;
  p[D_B1] = pow(10.0, -0.4 * mag1) / integrateWholeDisk1(r1) ;
  p[D_B2] = pow(10.0, -0.4 * mag2) / integrateWholeDisk2(r2) ;
  p[D_SIN_I] = sin_i ;
  p[D_TMO] = 0.5 * (time0 - omega) ;
  p[D_TPO] = 0.5 * (time0 + omega) ;

  if ((p[D_ECC] < 0.0) || (p[D_ECC] >= 1.0)) return (-1) ; 

  if ((p[D_R1] <= 0.0) || (p[D_R1] >= (1.0 - p[D_ECC]))) return (-2) ;

  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_ECC] - p[D_R1]))) return (-3) ;

  if  (p[D_SIN_I] > 1.0) return (-4) ;  // (p[D_SIN_I] <= sqrt(1.0 - sqr(p[D_R1] + p[D_R2])))   <- using this causes a selection effect for large (massive/evolving) stars

  if (p[D_B1] <= 0.0) return (-5) ;

  if (p[D_B2] <= 0.0) return (-6) ;
  

  time = (double)(WRITE_SIZE-1) / WRITE_SIZE ;
  prevFlux = flux(time, p, MIN_INTEGRATION_STEP, 1) ; ;
  time = (double)(WRITE_SIZE-2) / WRITE_SIZE ;
  prevprevFlux = flux(time, p, MIN_INTEGRATION_STEP, 1) ;

  for (i = 0 ; i < WRITE_SIZE ; i++)
    {
      time = (double)i / WRITE_SIZE ;
      F = flux(time, p, MIN_INTEGRATION_STEP, 1) ;

      if (F < (platFlux - (magErr * F / MAGNITUDE_ERROR_FACTOR)))
	countNumInEclipse++ ;

      //------

      if ((prevprevFlux >= prevFlux) && (prevFlux < F))
	{
	  if (*peclipse1Flux == 0.0) 
	    *peclipse1Flux = prevFlux ;
	  else
	    if (*peclipse2Flux == 0.0) 
	      *peclipse2Flux = prevFlux ;
	    else
	      return (-1) ;
	}

      // printf ("%f %f %f     %f %f\n", F, prevFlux, prevprevFlux, *peclipse1Flux, *peclipse2Flux) ;

      prevprevFlux = prevFlux ;
      prevFlux = F ;

      //------
 
      F += (magErr * randGauss() * F / MAGNITUDE_ERROR_FACTOR) ;
      fprintf (foutFit, "%f  %f  %f\n", period * time, -2.5 * log10(F), magErr) ;
    }

  return (countNumInEclipse) ;
}


// Returns the number of points in the eclipse, above the noise level.
// If failed - returns a negative value
int getScore (int g, int i, int j, int indexLCcolor, double period, double combObsMag, double e,
	      double time0, double omega, double sin_i, float magErr, struct LimbDarkType *plimbDark,
	      struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS], FILE *foutFit, double platFlux,
	      double *peclipse1Flux, double *peclipse2Flux)
{
  double mag1, mag2, inv_a, r1, r2 ;
  double xLogg, quadA1, quadB1, quadA2, quadB2 ;

  calcMags (combObsMag, starGroup[g].Mag[indexLCcolor][j] - starGroup[g].Mag[indexLCcolor][i], &mag1, &mag2) ;

  // 0.237757 = R_sol * ((2pi)^2 * (day/sec)^-2 / (G*mass_sol))^(1/3)
  inv_a = 0.237757 * pow(period * period * (starGroup[g].mass[i] + starGroup[g].mass[j]), -0.3333333333333) ;
  r1 = inv_a * starGroup[g].radius[i] ;
  r2 = inv_a * starGroup[g].radius[j] ;	

  xLogg = 4.437542 + log10(starGroup[g].mass[i] / sqr(starGroup[g].radius[i])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor, xLogg, starGroup[g].Teff[i], plimbDark, &quadA1, &quadB1))
    {
      printf ("ERROR LD1  (%d %d %d)\n", i, j, g) ;
      return (-10) ;
    }

  xLogg = 4.437542 + log10(starGroup[g].mass[j] / sqr(starGroup[g].radius[j])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor, xLogg, starGroup[g].Teff[j], plimbDark, &quadA2, &quadB2))
    {
      printf ("ERROR LD2  (%d %d %d)\n", i, j, g) ;
      return (-20) ;
    }
  
  return (writeLC_sim (foutFit, period, e, r1, r2, mag1, mag2, sin_i, time0, omega, magErr, quadA1, quadB1, quadA2, quadB2, platFlux, peclipse1Flux, peclipse2Flux)) ;
}


double getLightCurveFlux (float time, int g, int i, int j, int indexLCcolor, double period,
			 double combObsMag, double e, double time0, double omega, double sin_i,
			 struct LimbDarkType *plimbDark, struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS])
{
  double mag1, mag2, inv_a, r1, r2 ;
  double xLogg, quadA1, quadB1, quadA2, quadB2 ;
  double p[DIM] ;
  
  calcMags (combObsMag, starGroup[g].Mag[indexLCcolor][j] - starGroup[g].Mag[indexLCcolor][i], &mag1, &mag2) ;

  // 0.237757 = R_sol * ((2pi)^2 * (day/sec)^-2 / (G*mass_sol))^(1/3)
  inv_a = 0.237757 * pow(period * period * (starGroup[g].mass[i] + starGroup[g].mass[j]), -0.3333333333333) ;
  r1 = inv_a * starGroup[g].radius[i] ;
  r2 = inv_a * starGroup[g].radius[j] ;	

  xLogg = 4.437542 + log10(starGroup[g].mass[i] / sqr(starGroup[g].radius[i])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor, xLogg, starGroup[g].Teff[i], plimbDark, &quadA1, &quadB1))
    {
      printf ("ERROR LD1  (%d %d %d)\n", i, j, g) ;
      return (1) ;
    }

  xLogg = 4.437542 + log10(starGroup[g].mass[j] / sqr(starGroup[g].radius[j])) ; // 4.437542 = log(g_sol)
  if (interpolateLimbDark (indexLCcolor, xLogg, starGroup[g].Teff[j], plimbDark, &quadA2, &quadB2))
    {
      printf ("ERROR LD2  (%d %d %d)\n", i, j, g) ;
      return (2) ;
    }
  
  setLimbDarkening (quadA1, quadB1, quadA2, quadB2) ;
  
  omega /= 360.0 ;

  p[D_ECC] = e  ;
  p[D_R1] = r1 ;
  p[D_R2] = r2 ;
  p[D_B1] = pow(10.0, -0.4 * mag1) / integrateWholeDisk1(r1) ;
  p[D_B2] = pow(10.0, -0.4 * mag2) / integrateWholeDisk2(r2) ;
  p[D_SIN_I] = sin_i ;
  p[D_TMO] = 0.5 * (time0 - omega) ;
  p[D_TPO] = 0.5 * (time0 + omega) ;

  if ((p[D_ECC] < 0.0) || (p[D_ECC] >= 1.0)) return (-1) ; 
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= (1.0 - p[D_ECC]))) return (-2) ;
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_ECC] - p[D_R1]))) return (-3) ;
  if (p[D_SIN_I] > 1.0) return (-4) ; 
  if (p[D_B1] <= 0.0) return (-5) ;
  if (p[D_B2] <= 0.0) return (-6) ;
 
  return (flux(mod1(time/period), p, MIN_INTEGRATION_STEP, 1)) ;
}


int getRandArrVal (float *arr, int len, float truncateRange)
{
  int index = 0 ;
  float goal = arr[0] + (truncateRange * (arr[len-1] - arr[0]) * (float)rand() / RAND_MAX) ;

  while ((index < len) && (arr[index] < goal))
    index++ ;

  if (index == len)
    return (len-1) ;

  if (index == 0)  // shouldn't happen
    return (0) ;

  if ((arr[index] - goal) > (goal - arr[index-1]))
    return (index-1) ;
  
  return (index) ;
}


int getRandArrAgeVal (struct StarGroupType *starGroup, int len)
{
  int index = 0 ;
  float goal = starGroup[0].age + ((starGroup[len-1].age - starGroup[0].age) * (float)rand() / RAND_MAX) ;

  while ((index < len) && (starGroup[index].age < goal))
    index++ ;

  if (index == len)
    return (len-1) ;

  if (index == 0)  // shouldn't happen
    return (0) ;

  if ((starGroup[index].age - goal) > (goal - starGroup[index-1].age))
    return (index-1) ;
  
  return (index) ;
}



int main (int argc, char **argv)
{
  int fileNum, g, i, j, n, numGroups, indexLCcolor, res, numOutputLCs, numColors ;
  char deb_filename[MAX_STRLEN] ;
  float magErr, colorErr = 0.0 ;
  double period, ecc, time0, omega, sin_i, combObsMag, combColor ;
  double platTime, platFlux, dMag1, dMag2 ;
  FILE *foutLog, *foutDBLC, *foutFit ; 
  struct StarGroupType starGroup[MAX_NUM_STAR_GROUPS] ;
  struct LimbDarkType limbDark ;
  struct ColorInfoType colorArr[MAX_NUM_COLORS] ;
  double eclipse1Flux, eclipse2Flux, r1, r2, inv_a, inclination ;


  if ((argc < 8) || (argc == 9) || (argc > 9 + MAX_NUM_COLORS))
    {
      if (argc == 9)
	printf ("\nMake sure you included the <color err> parameter\n\n") ;

      printf ("\n Usage:  %s  <limb-darkening DB>  <isochrone DB>  <output dblc file>  <output log file>  <output LC file path/prefix>  <mag error>  <num output LCs>  [<color err> <color1> [color2]...] \n", argv[0]) ;
      printf ("Note: each color must have the format: '[filter]-[filter]', e.g. 'K-H'\n") ;
      return (1) ;
    }

  magErr = atof(argv[6]) ;
  if (magErr <= 0.0)
    {
      printf ("ERROR: Invalid mag error (%f)\n", magErr) ;
      return (2) ;
    }
    
  numOutputLCs = atoi(argv[7]) ;
  if (numOutputLCs <= 0)
    {
      printf ("ERROR: Invalid num output LCs (%d)\n", numOutputLCs) ;
      return (3) ;
    }
  
  if (argc > 9)
    {
      colorErr = atof(argv[8]) ;
      if (colorErr <= 0.0)
	{
	  printf ("ERROR: Invalid color error ('%s'), make sure you included the <color err> parameter.\n", argv[8]) ;
	  return (4) ;
	}

      numColors = argc - 9 ;
      
      for (n = 0 ; n < numColors ; n++)
	{
	  colorArr[n].filterPos = colorConvert (argv[n+9][0]) ;
	  colorArr[n].filterNeg = colorConvert (argv[n+9][2]) ;
	}
    }
  else numColors = 0 ;

  if (!strcmp(argv[1], argv[2]) ||
      !strcmp(argv[1], argv[3]) ||
      !strcmp(argv[1], argv[4]) ||
      !strcmp(argv[2], argv[3]) ||
      !strcmp(argv[2], argv[4]) ||
      !strcmp(argv[3], argv[4]))
    {
      printf ("ERROR: All the input/output filenames must be different.\n") ;
      return (5) ;
    }
 
  //-------------------

  n = loadLimbDark (argv[1], &limbDark) ;
  if (n)
    return (6) ;

  numGroups = loadIsochrones (argv[2], starGroup) ;
  if (numGroups <= 0)
    return (7) ;

  foutDBLC = fopen (argv[3], "wt") ;
  if (!foutDBLC)
    return (8) ;

  foutLog = fopen (argv[4], "wt") ;
  if (!foutLog)
    return (9) ;
    
  //-----------------------------------------------------
 

  srand (RAND_SEED) ;

  for (fileNum = 0 ; fileNum < numOutputLCs ; fileNum++)
    {  
      period = 10.0 * (double)rand() / ((double)RAND_MAX+1.0) ;
      time0 = (double)rand() / ((double)RAND_MAX+1.0) ;
      combObsMag = 0.0 ;

#ifdef DO_ECCENTRICITY 
      ecc = 0.1 * (double)rand() / ((double)RAND_MAX+1.0) ;
      omega = 2.0 * M_PI * (double)rand() / ((double)RAND_MAX+1.0) ;
#else
      ecc = 0.0 ;
      omega = 0.0 ;
#endif

      indexLCcolor = colorConvert (LC_FILTER) ; 

      g = getRandArrAgeVal(starGroup, numGroups) ;
      i = getRandArrVal(starGroup[g].mass, starGroup[g].arrLen, TRUNCATE_RANGE) ;
      j = getRandArrVal(starGroup[g].mass, starGroup[g].arrLen, TRUNCATE_RANGE) ;

      //-- calc sin(a) ---
      // 0.237757 = R_sol * ((2pi)^2 * (day/sec)^-2 / (G*mass_sol))^(1/3)
      inv_a = 0.237757 * pow(period * period * (starGroup[g].mass[i] + starGroup[g].mass[j]), -0.3333333333333) ;
      r1 = inv_a * starGroup[g].radius[i] ;
      r2 = inv_a * starGroup[g].radius[j] ;

      if ((r1 + r2) >= DETACHED_MAX_r1r2)
	{
	  fileNum-- ;
	  printf ("skip contact\n") ;
	  continue ;
	}

      inclination = (0.5 * M_PI) - (((0.5 * M_PI) - acos(r1 + r2)) * (double)rand() / ((double)RAND_MAX+1.0)) ;

      sin_i = sin (inclination) ;

      //------------------

      sprintf (deb_filename, "%s_%02d.lc", argv[5], fileNum) ;
      foutFit = fopen (deb_filename, "wt") ;
      if (foutFit)
	{
	  platTime = period * (time0 - (omega / (2.0 * M_PI))) ;
	  platFlux = getLightCurveFlux (platTime, g, i, j, indexLCcolor, period, combObsMag, ecc, time0, omega, sin_i, &limbDark, starGroup) ;

	  res = getScore (g, i, j, indexLCcolor, period, combObsMag, ecc, time0, omega, sin_i, magErr, &limbDark, starGroup, foutFit, platFlux, &eclipse1Flux, &eclipse2Flux) ;
	  fclose (foutFit) ;
	     
	  if ((res < 0) || (eclipse1Flux == 0.0) || (eclipse2Flux == 0.0))
	    {
	      fileNum-- ;
	      printf ("skip LC doesn't have 2 eclipses  %d  %f %f\n", res, eclipse1Flux, eclipse2Flux) ;
	      printf ("%f %f %f %f %f %f\n", time0, omega, r1, r2, ecc, sin_i) ;
	      continue ;
	    }

	  dMag1 = -2.5 * log10 (eclipse1Flux / platFlux) ;
	  dMag2 = -2.5 * log10 (eclipse2Flux / platFlux) ;
	  
	  if (dMag1 < dMag2)
	    swap (&dMag1, &dMag2);
	  
	  if ((dMag1 < (NOISE_THRESHOLD * magErr)) || (dMag2 < (NOISE_THRESHOLD * magErr)))
	    {
	      fileNum-- ;
	      printf ("skip noise\n") ;
	      continue ;
	    }
	  
	  printf ("writing LC %d [%f %f %f %f %f %f %f]\n", fileNum, period, ecc, time0, omega, sin_i, r1, r2) ;
	  
	  fprintf (foutLog, "%d %f %f %f   %f  %d %f %f\n", fileNum, starGroup[g].mass[i], starGroup[g].mass[j], starGroup[g].age, sin_i, res, dMag1 / magErr, dMag2 / magErr) ;
	  
	  //------------------------------
	  
	  fprintf (foutDBLC, "%s %f %f 0 0 0 0 0 %f 0 %f 0 0 0 %f 0 %f 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n", deb_filename, period, ecc, combObsMag, 1000.0, time0, omega) ;
	  fprintf (foutDBLC, "%c %d 0 0 0 0 0\n", LC_FILTER, numColors) ;
	  
	  for (n = 0 ; n < numColors ; n++)
	    {
	      combColor = 
		combineMag(starGroup[g].Mag[colorArr[n].filterPos][i],
			   starGroup[g].Mag[colorArr[n].filterPos][j]) -
		combineMag(starGroup[g].Mag[colorArr[n].filterNeg][i],
			   starGroup[g].Mag[colorArr[n].filterNeg][j]) ;
	      
	      combColor += (colorErr * randGauss() * combColor / MAGNITUDE_ERROR_FACTOR) ;
	      
	      fprintf (foutDBLC,"%c-%c %f %f\n", argv[n+9][0], argv[n+9][2], combColor, colorErr) ;
	    }
	}
      else
	printf ("ERROR: Couldn't write the model fit file\n") ;      
    }


  for (i = 0 ; i < numGroups ; i++)
    {
      free (starGroup[i].mass) ;
      free (starGroup[i].Teff) ;
      free (starGroup[i].radius) ;
      
      for (j = 0 ; j < NUM_FILTERS ; j++)
	free (starGroup[i].Mag[j]) ;
    }
  
  fclose (foutDBLC) ;
  fclose (foutLog) ;
  return (0) ;
}
