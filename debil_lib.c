/************************************************************************
 * DEBiL library of assorted light curve utilities
 *
 * Written by Jonathan Devor (jdevor@cfa.harvard.edu)
 *    Harvard-Smithsonian Center for Astrophysics 
 *    60 Garden St., Cambridge, MA 02138 USA
 *
 ***********************************************************************/

// Smoothing kernel parameter:  for ridOutliers()
#define NUM_POINTS_HALF_KERNEL 8  // Must be smaller than a half a dip width, but at least 2

#define MIN_INTEGRATION_STEP 0.05   // Recomended range:  0.05-0.01 ; This should be about 0.1 * sqrt(amp error)

// Eccentric anomaly calculation parameter:
#define ECC_ANOMALY_MAX_ERROR 1.0e-4  // Convergence threshold (average error is less than a third of this) 

// Eccentricity calculation parameter:
#define ECCENTRICITY_MAX_ERROR 1.0e-6   

// General bisection implementation parameter:
// only needed for pathological cases (e.g. overflows or huge roundoff errors)
#define BISECTION_SAFETY_MAX_ITERATIONS 1000 

#define EPSILON 1.0e-9  // To avoid division by zero

// A useful constant for magnitude error conversion (2.5 / ln(10))
#define MAGNITUDE_ERROR_FACTOR 1.085736205

#define WRITE_CURVE_FIT 500   // The number of fitted points to write

#define ECLIPSE_SEARCH_STEP 0.0002
#define FLUX_EPSILON 1.0e-6

// ======================= Downhill simplex ================================

#define DIM 8  // The number of dimensions (don't change!)

// list of dimensions:
// physical-
#define D_ECC   0  // Eccentricity of orbit (e)
#define D_R1    1  // Star's radius (in units of the distance between them)
#define D_R2    2
#define D_B1    3  // Central surface brightness (flux per unit disk area)
#define D_B2    4
// orientation-
#define D_SIN_I 5  // sin(inclination)
#define D_TMO 6    // (time / 2) - (omega / 4pi)  [ or  (time - (omega / 2pi))  when no ecc. ]
#define D_TPO 7    // (time / 2) + (omega / 4pi)

// Time is the time (scan parameter) at perihelion and omega is the argument of perihelion
// (i.e. azimuth of eccentricity alignment
// Note that time and omega become degenerate for small eccentricities. For this
// reason I combined them into orthogonal parameters ,D_TMO and D_TPO, where at
// least the former can be accurately determined and therefor will probably converge
// much faster in the fitting algorithm.

// NOTE: [valid range of physical parameters]
//  0 <= time0 < 1       (not vital)
//  0 <= omega < 2*pi    (not vital)
//  0 < sin_i <= 1
//  0 <= e < 1           (closed orbit)
//  0 < r1 < 1-e         (prevent collisions)
//  0 < r2 < 1-e-r1      (prevent collisions)
//  0 < B1
//  0 < B2

// ==================== Includes ===================================

#include <math.h>
#include <stdlib.h> // For malloc() + qsort()
#include <values.h> // For MAXFLOAT
#include <stdio.h>  // For fopen()
#include <string.h> // For memmove()

// ================ General utility functions ======================


// Returns the square of a given x
static double sqr (double x)
{
  return (x * x) ;
}

// Returns the cube of a given x
static double cube (double x)
{
  return (x * x * x) ;
}


// Returns the modulo 1 value of x
static double mod1 (double x)
{
  return (x - floor(x)) ;
}


// -----------------------

// Prints a given error message (errStr) reported by a given function (functionName) about a
// given light curve file (dataFilename). It is added to a given error repository file (errFilename)
// isFatal: 1 = Error - go to next  ;  0 = Warning - try anyways  ; -1 = Note - just informative
void printError (char *errFilename, int typeID, char *dataFilename, char *functionName, char *errStr)
{
  FILE *ferr = fopen (errFilename, "at") ;

  if (!ferr)
    {
      printf ("%s: [ERROR] Couldn't open the output error file ('%s').\n", dataFilename, errFilename) ;
     
      switch (typeID)
	{
	case 1:	 printf ("%s: [ERROR] %s - %s\n", dataFilename, functionName, errStr) ; break ;
	case 0:	 printf ("%s: [Warning] %s - %s\n", dataFilename, functionName, errStr) ; break ;
	case -1: printf ("%s: [Note] %s - %s\n", dataFilename, functionName, errStr) ; break ;
	default: printf ("%s: [ID #%d)] %s - %s\n", dataFilename, typeID, functionName, errStr) ;
	}

      return ;
    }

  switch (typeID)
    {
    case 1:  fprintf (ferr, "%s: [ERROR] %s - %s\n", dataFilename, functionName, errStr) ; break ;
    case 0:  fprintf (ferr, "%s: [Warning] %s - %s\n", dataFilename, functionName, errStr) ; break ;
    case -1: fprintf (ferr, "%s: [Note] %s - %s\n", dataFilename, functionName, errStr) ; break ;
    default: fprintf (ferr, "%s: [ID #%d] %s - %s\n", dataFilename, typeID, functionName, errStr) ;
    }    
  
  fclose (ferr) ;
}


//--------------------------------------------------------

// Allocates and sets: *ptime, *pamp, *perr
// Returns the number of light curve measurements
// Returns a negative number upon error, and does not perform allocation
int loadLC (char *filename, double period, char *errfilename, float **ptime, float **pamp, float **perr)
{
  FILE *finLC ;
  int size, sizeOrig ;
  double tmpTime, tmpMag, tmpErr ;

  if (period <= 0.0)
    {
      printError (errfilename, 1, filename, "loadLC()", "Invalid period") ;
      return (-1) ;
    }
  
  if (!(finLC = fopen (filename, "rt")))
    {
      printError (errfilename, 1, filename, "loadLC()", "Couldn't open the input LC file") ;
      return (-2) ;
    }

  //----------------------------------------------------------------

  sizeOrig = 0 ;
  while (1 == fscanf(finLC, "%*f %*f %lf\n", &tmpErr))
    if (tmpErr > 0.0)  // Sign of an invalid magnitude
      sizeOrig++ ;

  if (sizeOrig == 0)
    {
      printError (errfilename, 1, filename, "loadLC()", "Not enough data points") ;
      fclose (finLC) ;
      return (-3) ;
    }

  *ptime = (float*)malloc (sizeOrig * sizeof(float)) ;
  *pamp = (float*)malloc (sizeOrig * sizeof(float)) ;
  *perr = (float*)malloc (sizeOrig * sizeof(float)) ;

  if (!(*ptime) || !(*pamp) || !(*perr))
    {
      printError (errfilename, 1, filename, "loadLC()", "Not enough memory") ;
      fclose(finLC) ;
      if (*ptime) free (*ptime) ;
      if (*pamp) free (*pamp) ;
      if (*perr) free (*perr) ;
      return (-4) ;
    }

  rewind(finLC) ;
  size = 0 ;
  while ((size < sizeOrig) && (3 == fscanf(finLC, "%lf %lf %lf\n", &tmpTime, &tmpMag, &tmpErr)))
    if (tmpErr > 0.0)  // Sign of an invalid magnitude
      {
	(*ptime)[size] = (float)mod1(tmpTime / period) ;
	(*pamp)[size] = (float)pow(10.0, -0.4 * tmpMag) ;  // Convert magnitude (logarithmic) to amplitude (linear)
	(*perr)[size] = (float)((*pamp)[size] * tmpErr / MAGNITUDE_ERROR_FACTOR) ;  // Convert to absolute error
	size++ ;
      }

  fclose(finLC) ;

  if (size != sizeOrig) 
    {
      printError (errfilename, 1, filename, "loadLC()", "Read number mismatch") ;
      free (*ptime) ;
      free (*pamp) ;
      free (*perr) ;
      return (-5) ;	
    }

  return (size) ;
}


// Same as loadLC(), but it also returns the midpoint of the timeStamps, without modulo (i.e. an epoch)
// Allocates and sets: *ptime, *pamp, *perr
// Returns the number of light curve measurements
// Returns a negative number upon error, and does not perform allocation
int loadLC_midTime (char *filename, double period, char *errfilename, float **ptime, float **pamp, float **perr, double *pmidTime)
{
  FILE *finLC ;
  int size, sizeOrig ;
  double tmpTime, tmpMag, tmpErr, maxTime = 0.0, minTime = 0.0 ;

  if (period <= 0.0)
    {
      printError (errfilename, 1, filename, "loadLC()", "Invalid period") ;
      return (-1) ;
    }
  
  if (!(finLC = fopen (filename, "rt")))
    {
      printError (errfilename, 1, filename, "loadLC()", "Couldn't open the input LC file") ;
      return (-2) ;
    }

  //----------------------------------------------------------------

  sizeOrig = 0 ;
  while (1 == fscanf(finLC, "%*f %*f %lf\n", &tmpErr))
    if (tmpErr > 0.0)  // Sign of an invalid magnitude
      sizeOrig++ ;

  if (sizeOrig == 0)
    {
      printError (errfilename, 1, filename, "loadLC()", "Not enough data points") ;
      fclose (finLC) ;
      return (-3) ;
    }

  *ptime = (float*)malloc (sizeOrig * sizeof(float)) ;
  *pamp = (float*)malloc (sizeOrig * sizeof(float)) ;
  *perr = (float*)malloc (sizeOrig * sizeof(float)) ;

  if (!(*ptime) || !(*pamp) || !(*perr))
    {
      printError (errfilename, 1, filename, "loadLC()", "Not enough memory") ;
      fclose(finLC) ;
      if (*ptime) free (*ptime) ;
      if (*pamp) free (*pamp) ;
      if (*perr) free (*perr) ;
      return (-4) ;
    }

  rewind(finLC) ;
  size = 0 ;
  while ((size < sizeOrig) && (3 == fscanf(finLC, "%lf %lf %lf\n", &tmpTime, &tmpMag, &tmpErr)))
    if (tmpErr > 0.0)  // Sign of an invalid magnitude
      {
	(*ptime)[size] = (float)mod1(tmpTime / period) ;
	(*pamp)[size] = (float)pow(10.0, -0.4 * tmpMag) ;  // Convert magnitude (logarithmic) to amplitude (linear)
	(*perr)[size] = (float)((*pamp)[size] * tmpErr / MAGNITUDE_ERROR_FACTOR) ;  // Convert to absolute error

	if (size == 0)
	  {
	    maxTime = tmpTime ;
	    minTime = tmpTime ;
	  }
	else
	  {
	    if (tmpTime > maxTime) maxTime = tmpTime ;
	    else if (tmpTime < minTime) minTime = tmpTime ;
	  }

	size++ ;
      }

  fclose(finLC) ;

  if (size != sizeOrig) 
    {
      printError (errfilename, 1, filename, "loadLC()", "Read number mismatch") ;
      free (*ptime) ;
      free (*pamp) ;
      free (*perr) ;
      return (-5) ;	
    }

  *pmidTime = 0.5 * (minTime + maxTime) ;
  return (size) ;
}



//--------------------------------------------------------

// Data structure of a single light curve observation
struct LCdata
{
  float primary ;
  float secondary ;
  float tertiary ;
};


// A comparator for a pair-o-float
// Returns: -1 if x.primary < y.primary
//           0 if x.primary = y.primary
//           1 if x.primary > y.primary
static int comparLCdata (const void *x, const void *y)
{
  if (((struct LCdata*)x)->primary < ((struct LCdata*)y)->primary)
    return (-1) ;
  else
    return (((struct LCdata*)x)->primary > ((struct LCdata*)y)->primary) ;
}


// Returns 0 = success ; 1 = failure (malloc)
int sortSamples (float *time, float *amp, float *err, int size)
{
  int i ;
  struct LCdata *arr = (struct LCdata*)malloc (size * sizeof(struct LCdata)) ;

  if (!arr)
    return (1) ;

  for (i = 0 ; i < size ; i++)
    {
      arr[i].primary = time[i] ;
      arr[i].secondary = amp[i] ;
      arr[i].tertiary = err[i] ;
    }

  qsort (arr, size, sizeof(struct LCdata), comparLCdata) ;

  for (i = 0 ; i < size ; i++)
    {
      time[i] = arr[i].primary ;
      amp[i] = arr[i].secondary ;
      err[i] = arr[i].tertiary ;
    }

  free (arr) ;
  return (0) ;
}

//-------------

// Performs one iteration of statistical summations
static void updateVars (float x, float y,
			double *A1, double *A2, double *A3, double *A4,
			double *B0, double *B1, double *B2, double *C0)
{
  const double x2 = x * x ;

  (*A1) += x ;
  (*A2) += x2 ;
  (*A3) += x * x2 ;
  (*A4) += x2 * x2 ;
  (*B0) += y ;
  (*B1) += y * x ;
  (*B2) += y * x2 ;
  (*C0) += y * y ;
}


// Calculates a second order regression (parabola:  a*x^2 + b*x + c) fit to given data (time, amp)
// around a given instant (sampTime) going up and down by a certain number of data points
// to minimize round off errors, time[] is centered around sampTime and the results are corrected
// for the shift at the end.
// Assumes that the data are sorted in time
// Output: pointers to: a b c  and the variance around the fit (pVariance, if not null)
static void regressionOrder2 (float sampTime, float *time, float *amp, int size,
			      double *pa, double *pb, double *pc, double *pVariance)
{
  int i, indexDown = 0, indexUp = size - 1 ;
  int A0 = NUM_POINTS_HALF_KERNEL + NUM_POINTS_HALF_KERNEL ;
  double A1 = 0.0, A2 = 0.0, A3 = 0.0, A4 = 0.0 ;
  double B0 = 0.0, B1 = 0.0, B2 = 0.0, C0 = 0.0 ;
  double denom, a, b, c ;

  // Step 1: find the closes index above it
  // Step 1.1: run the bisection algorithm
  while ((indexDown + 1) < indexUp)
    {
      i = (indexDown + indexUp) / 2 ;

      if (time[i] > sampTime)
	indexUp = i ;
      else
	indexDown = i ;
    }

  // Step 1.2: take care of the ends
  if (time[indexDown] > sampTime)
    indexUp = 0 ;
  else
    if (time[indexUp] <= sampTime)
      indexDown = size - 1 ;


  // Step 2: do statistics
  if (indexUp < NUM_POINTS_HALF_KERNEL)
    {
      for (i = indexUp + size - NUM_POINTS_HALF_KERNEL ; i < size ; i++)
	updateVars (time[i] - 1.0 - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;

      for (i = 0 ; i < indexUp + NUM_POINTS_HALF_KERNEL ; i++)
	updateVars (time[i] - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;
    }
  else if (indexDown >= (size - NUM_POINTS_HALF_KERNEL))
    {
      for (i = indexDown + 1 - NUM_POINTS_HALF_KERNEL ; i < size ; i++)
	updateVars (time[i] - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;

      for (i = 0 ; i <= (indexDown + NUM_POINTS_HALF_KERNEL - size) ; i++)
	updateVars (time[i] + 1.0 - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;
    }
  else  // Normal case
    {
      for (i = indexUp - NUM_POINTS_HALF_KERNEL ; i < (indexUp + NUM_POINTS_HALF_KERNEL) ; i++)
	updateVars (time[i] - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;
    }

  denom = (A2 * ((A0*A4) + (2.0*A1*A3) - (A2*A2))) - (A0*A3*A3) - (A1*A1*A4) ;  // Denominator

  if (denom == 0.0)
    denom = EPSILON ;    // Prevents division by zero

  *pa = ((A0 * ((A2*B2) - (A3*B1))) + (A1 * ((A3*B0) - (A1*B2))) + (A2 * ((A1*B1) - (A2*B0)))) / denom ;
  *pb = ((A0 * ((A4*B1) - (A3*B2))) + (A1 * ((A2*B2) - (A4*B0))) + (A2 * ((A3*B0) - (A2*B1)))) / denom ;
  *pc = ((A1 * ((A3*B2) - (A4*B1))) + (A2 * ((A4*B0) - (A2*B2))) + (A3 * ((A2*B1) - (A3*B0)))) / denom ;

  if (pVariance)
    {
      a = *pa ;
      b = *pb ;
      c = *pc ;
      *pVariance = ((A4*a*a) + (2.0*a*b*A3) + (A2 * ((b*b) + (2.0*a*c))) + (2.0*b*c*A1) + (A0*c*c) - (2.0 * ((B2*a) + (B1*b) + (B0*c))) + C0) / (A0-3) ;
    }

  // Moves back the "zero point":  a(x-x0)^2 + b(x-x0) + c  -->  ax^2 + bx + c
  *pc += (((*pa) * sampTime * sampTime) - ((*pb) * sampTime)) ;
  *pb -= (2.0 * (*pa) * sampTime) ;
}


// Removes the outliers from the input data array
// by comparing the data against a 2nd order regression spline.
// the data is assumed to be sorted by time (modulo period = phase)
// even when points are removes, the sorted order is maintained
// note that it will keep "lonely points" with no other data points near it
// time - input: array of time stamps
// amp - input: matching array of amplitudes
// size - size of input arrays
// *psize = outputs the new size of the array (will be reduced if outliers are found)
// varianceLimit = "sigma clipping" threshold, after which data points are considered outliers (in units of variance)
// Returns the reduced chi2 around the spline (-1.0 in case of error)
double ridOutliers (float *time, float *amp, float *err, int *psize, float varianceLimit)
{
  int index, doRid, numVariance ;
  double x, sumVariance, variance, varianceFit, sumChi2 ;
  double a, b, c ;

  do
    {
      doRid = 0 ;
      index = 0 ;
      numVariance = 0 ;
      sumVariance = 0.0 ;
      sumChi2 = 0.0 ;

      while (index < (*psize))
	{
	  x = time[index] ;

	  regressionOrder2 (x, time, amp, *psize, &a, &b, &c, &varianceFit) ;

	  variance = sqr((a * x * x) + (b * x) + c - amp[index]) ;

	  if (variance > (varianceFit * varianceLimit))
	    {
	      doRid = 1 ;
	      (*psize)-- ;

	      if (index < (*psize))
		{
		  memmove (&time[index], &time[index+1], ((*psize) - index) * sizeof(float)) ;
		  memmove (&amp[index], &amp[index+1], ((*psize) - index) * sizeof(float)) ;
		  memmove (&err[index], &err[index+1], ((*psize) - index) * sizeof(float)) ;
		}
	    }
	  else
	    {
	      sumVariance += variance ;
	      sumChi2 += variance / sqr(err[index]) ;
	      numVariance++ ;
	      index++ ;
	    }
	}
    }
  while (doRid) ;

  if (numVariance == 0)
    return (-1.0) ;

  return (sumChi2 / numVariance) ;
}


//========================= Eccentricity ======================================

// Uses a combined Newton-Raphson + bisection method for finding the root (E) of:
//  M = E - (e * sin(E))     [Kepler's equation]
// Given: M = mean anomaly ; e = eccentricity of orbit (0 <= e < 1)
// Returns: E = Eccentric Anomaly
static double eccentricAnomaly (double M, double e)
{
  int counter = BISECTION_SAFETY_MAX_ITERATIONS ;
  double Emin = M - 1.0, Emax = M + 1.0, E = M ;
  double f, dfm, dE ;

  do
    {
      f = M + (e * sin(E)) - E ;  // May be optimized by setting cos=sqrt(1-sin^2)
      dfm = 1.0 - (e * cos(E)) ;  // Minus differential of f (always non-negative)
      dE = dfm * ECC_ANOMALY_MAX_ERROR ;

      if (f > dE)
	{
	  Emin = E ;
	  dE = Emax - Emin ;

	  if ((dfm * dE) > f)
	    E += (f / dfm) ;
	  else
	    E = 0.5 * (Emax + Emin) ;
	}
      else if (f < (-dE))
	{
	  Emax = E ;
	  dE = Emax - Emin ;

	  if ((dfm * dE) > -f)
	    E += (f / dfm) ;
	  else
	    E = 0.5 * (Emax + Emin) ;
	}
      else return (E) ;
    }
  while ((dE > ECC_ANOMALY_MAX_ERROR) && ((--counter) > 0)) ;

  return (E) ;  // should almost never exits here
}

//========================== Limb Darkening ====================================

// Quadratic limb darkening coefficients for a solar-like star in I filter (Claret 1998)
// Note that a "square root" law give a better fit (especially for I filter),
// but the error in assuming a solar star, is much larger and the total CPU requirements
// would almost double (doing a linear fit, on the other hand, won't save much CPU).

// Note: Doubling the limbConst's and fucntions is inelegant, but needed for optimization.

static double limbConstA1 = 0.0 ;
static double limbConstB1 = 0.0 ;
static double limbConstC1 = 0.0 ;
static double limbConstD1 = 0.0 ;
static double limbConstE1 = 0.0 ;
static double limbConstF1 = 0.0 ;
static double limbConstG1 = 0.0 ;
static double limbConstH1 = 0.0 ;
static double limbConstA2 = 0.0 ;
static double limbConstB2 = 0.0 ;
static double limbConstC2 = 0.0 ;
static double limbConstD2 = 0.0 ;
static double limbConstE2 = 0.0 ;
static double limbConstF2 = 0.0 ;
static double limbConstG2 = 0.0 ;
static double limbConstH2 = 0.0 ;


// Sets the quadratic limb darkening parameters to arbitrary given values.
static void setLimbDarkening (double quadA1, double quadB1, double quadA2, double quadB2)
{
  limbConstA1 = 1.0 - quadA1 - quadB1 ;
  limbConstB1 = quadA1 + quadB1 + quadB1 ;
  limbConstC1 = -quadB1 ;
  limbConstD1 = M_PI * (1.0 - quadA1 - quadB1 - quadB1) ;
  limbConstE1 = 0.5 * M_PI * quadB1 ;
  limbConstF1 = 2.0 * M_PI * (quadA1 + quadB1 + quadB1) / 3.0 ;
  limbConstG1 = M_PI * (1.0 - quadA1 - quadB1) ;
  limbConstH1 = M_PI * (6.0 - quadA1 - quadA1 - quadB1) / 6.0 ;

  limbConstA2 = 1.0 - quadA2 - quadB2 ;
  limbConstB2 = quadA2 + quadB2 + quadB2 ;
  limbConstC2 = -quadB2 ;
  limbConstD2 = M_PI * (1.0 - quadA2 - quadB2 - quadB2) ;
  limbConstE2 = 0.5 * M_PI * quadB2 ;
  limbConstF2 = 2.0 * M_PI * (quadA2 + quadB2 + quadB2) / 3.0 ;
  limbConstG2 = M_PI * (1.0 - quadA2 - quadB2) ;
  limbConstH2 = M_PI * (6.0 - quadA2 - quadA2 - quadB2) / 6.0 ;
}


// cosAngle2 = the square of the cosine of the angle between the zenith of
// the star to the given point on the surface (relative to the star's center)
// Returns the solar limb darkening coefficient
// Note that it is normalized so that at the star's zenith (cosAngle2 = 1), it returns 1
static double limbDarkening1 (double cosAngle2)
{
  return (limbConstA1 + (limbConstB1 * sqrt(cosAngle2)) + (limbConstC1 * cosAngle2)) ;
}

static double limbDarkening2 (double cosAngle2)
{
  return (limbConstA2 + (limbConstB2 * sqrt(cosAngle2)) + (limbConstC2 * cosAngle2)) ;
}


// Analytically integrates the limb darkened disk, from r=0 to r=rBack
// rBack = the full radius of the star's disk
static double integrateWholeDisk1 (double rBack)
{
  return (limbConstH1 * rBack * rBack) ;
}

static double integrateWholeDisk2 (double rBack)
{
  return (limbConstH2 * rBack * rBack) ;
}


// Analytically integrates the limb darkened disk, from 0 to the given finish
// finish = where to end the integration
// rBack = the radius of the back star's disk (the one seen as a crescent)
static double integrateDiskFromStart1 (double finish, double rBack)
{
  const double x = sqr(finish / rBack) ;

  return (rBack * rBack * (((limbConstD1 + (limbConstE1 * x)) * x) + (limbConstF1 * (1.0 - cube(sqrt(1.0 - x)))))) ;
}

static double integrateDiskFromStart2 (double finish, double rBack)
{
  const double x = sqr(finish / rBack) ;

  return (rBack * rBack * (((limbConstD2 + (limbConstE2 * x)) * x) + (limbConstF2 * (1.0 - cube(sqrt(1.0 - x)))))) ;
}


// Analytically integrates the limb darkened disk, from the given start to rBack
// start = from where to begin the integration
// rBack = the radius of the back star's disk (the one seen as a crescent)
static double integrateDiskToFinish1 (double start, double rBack)
{
  const double x = 1.0 - sqr(start / rBack) ;

  return (rBack * rBack * (((limbConstG1 - (limbConstE1 * x)) * x) + (limbConstF1 * cube(sqrt(x))))) ;
}

static double integrateDiskToFinish2 (double start, double rBack)
{
  const double x = 1.0 - sqr(start / rBack) ;

  return (rBack * rBack * (((limbConstG2 - (limbConstE2 * x)) * x) + (limbConstF2 * cube(sqrt(x))))) ;
}


// Makes a simple trapezoid numerical integration with varying step size
// (smaller, the closer an asymptote) for finding the brightness of
// a partially hidden star. This method seems to works better than Simpson's.
// start = where to begin the integration
// finish = where to end the integration
// rFront = the radius of the front star's disk (the one seen as a disk)
// rBack = the radius of the back star's disk (the one seen as a crescent)
// D = the distance between the centers of the two stars
// dr0 = the infinitesimal integration step scale (unitless)
// Assumes:  D + rFront >= finish >= start >= |D - rFront| >= 0
// [num iterations] ~= 4 / dr0   ;  [abs. error] ~= 0.1 * rFront^2 * dr0^2
static double integrateCrescent1 (double start, double finish, double rFront, double rBack,
			  double D, double dr0)
{
  const double rBack_m2 = 1.0 / (rBack * rBack) ;
  const double rFront2mD2 = (rFront * rFront) - (D * D) ;
  const double middle = 0.5 * (start + finish) ;
  const double D2 = 2.0 * D ;
  const double limitUp = (((D + rFront) < rBack) ? D + rFront : rBack) ;
  const double limitDown = fabs(D - rFront) ;
  double dr, r, r2, rStep, sum = 0.0 ;

  if (D < EPSILON)  // To avoid division by zero
    return (0.0) ;

  // Since dr0 is unitless and we want dr to have units of "distance" (1 = sum semi-major axes)
  dr0 *= sqrt(rFront) ;  

  // Step 1: integrate from middle up
  dr = dr0 * sqrt(limitUp - middle) ;
  for (rStep = middle + dr ; rStep < finish ; rStep += dr)
    {
      r = rStep - (0.5 * dr) ;
      r2 = r * r ;
      sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening1(1.0 - (rBack_m2 * r2)) * dr) ;

      dr = dr0 * sqrt(limitUp - r) ;
    }

  // Step 2: add sliver at upper edge
  r = 0.5 * (finish + rStep - dr) ;
  r2 = r * r ;
  dr += (finish - rStep) ;
  sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening1(1.0 - (rBack_m2 * r2)) * dr) ;

  // Step 3: integrate from middle down
  dr = dr0 * sqrt(middle - limitDown) ;
  for (rStep = middle - dr ; rStep > start ; rStep -= dr)
    {
      r = rStep + (0.5 * dr) ;
      r2 = r * r ;
      sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening1(1.0 - (rBack_m2 * r2)) * dr) ;

      dr = dr0 * sqrt(r - limitDown) ;
    }

  // Step 4: add sliver at bottom edge
  r = 0.5 * (start + rStep + dr) ;
  r2 = r * r ;
  dr += (rStep - start) ;
  sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening1(1.0 - (rBack_m2 * r2)) * dr) ;

  return (2.0 * sum) ;
}


static double integrateCrescent2 (double start, double finish, double rFront, double rBack,
			  double D, double dr0)
{
  const double rBack_m2 = 1.0 / (rBack * rBack) ;
  const double rFront2mD2 = (rFront * rFront) - (D * D) ;
  const double middle = 0.5 * (start + finish) ;
  const double D2 = 2.0 * D ;
  const double limitUp = (((D + rFront) < rBack) ? D + rFront : rBack) ;
  const double limitDown = fabs(D - rFront) ;
  double dr, r, r2, rStep, sum = 0.0 ;

  if (D < EPSILON)  // To avoid division by zero
    return (0.0) ;

  // Since dr0 is unitless and we want dr to have units of "distance" (1 = sum semi-major axes)
  dr0 *= sqrt(rFront) ;  

  // Step 1: integrate from middle up
  dr = dr0 * sqrt(limitUp - middle) ;
  for (rStep = middle + dr ; rStep < finish ; rStep += dr)
    {
      r = rStep - (0.5 * dr) ;
      r2 = r * r ;
      sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening2(1.0 - (rBack_m2 * r2)) * dr) ;

      dr = dr0 * sqrt(limitUp - r) ;
    }

  // Step 2: add sliver at upper edge
  r = 0.5 * (finish + rStep - dr) ;
  r2 = r * r ;
  dr += (finish - rStep) ;
  sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening2(1.0 - (rBack_m2 * r2)) * dr) ;

  // Step 3: integrate from middle down
  dr = dr0 * sqrt(middle - limitDown) ;
  for (rStep = middle - dr ; rStep > start ; rStep -= dr)
    {
      r = rStep + (0.5 * dr) ;
      r2 = r * r ;
      sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening2(1.0 - (rBack_m2 * r2)) * dr) ;

      dr = dr0 * sqrt(r - limitDown) ;
    }

  // Step 4: add sliver at bottom edge
  r = 0.5 * (start + rStep + dr) ;
  r2 = r * r ;
  dr += (rStep - start) ;
  sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening2(1.0 - (rBack_m2 * r2)) * dr) ;

  return (2.0 * sum) ;
}

// =============== Calculate the observed radiation flux from the binary system ===============

// Calculates the crescent of the back star, when it is partially eclipsed by the front star
// rFront = the radius of the star in front
// rBack = the radius of the star in back
// D = the distance between the centers of the two stars
// dr = the infinitesimal integration step
// Returns the relative brightness of the crescent of the back star
static double integrateBackOverlap1 (double rFront, double rBack, double D, double dr)
{
  if (rFront < D)
    {
      if (rBack > (D + rFront))
	return (integrateDiskFromStart1 (D - rFront, rBack) +
		integrateCrescent1 (D - rFront, D + rFront, rFront, rBack, D, dr) +
		integrateDiskToFinish1 (D + rFront, rBack)) ;                  // E

      return (integrateDiskFromStart1 (D - rFront, rBack) +
	      integrateCrescent1 (D - rFront, rBack, rFront, rBack, D, dr)) ;  // B
    }

  if (rFront > (D + rBack))
    return (0.0) ;                                                             // D

  if (rBack > (D + rFront))
    return (integrateCrescent1 (rFront - D, rFront + D, rFront, rBack, D, dr) +
	    integrateDiskToFinish1 (rFront + D, rBack)) ;                      // F

  return (integrateCrescent1 (rFront - D, rBack, rFront, rBack, D, dr)) ;      // C
}

static double integrateBackOverlap2 (double rFront, double rBack, double D, double dr)
{
  if (rFront < D)
    {
      if (rBack > (D + rFront))
	return (integrateDiskFromStart2 (D - rFront, rBack) +
		integrateCrescent2 (D - rFront, D + rFront, rFront, rBack, D, dr) +
		integrateDiskToFinish2 (D + rFront, rBack)) ;                  // E

      return (integrateDiskFromStart2 (D - rFront, rBack) +
	      integrateCrescent2 (D - rFront, rBack, rFront, rBack, D, dr)) ;  // B
    }

  if (rFront > (D + rBack))
    return (0.0) ;                                                             // D

  if (rBack > (D + rFront))
    return (integrateCrescent2 (rFront - D, rFront + D, rFront, rBack, D, dr) +
	    integrateDiskToFinish2 (rFront + D, rBack)) ;                      // F

  return (integrateCrescent2 (rFront - D, rBack, rFront, rBack, D, dr)) ;      // C
}


/******************************************************************
 * Returns the total flux from both stars of the binary.
 * Ignores limb darkening, gravity darkening, etc.
 * time = scan parameter 0 <= time < 1
 * p[DIM] = parameter vector
 * dr = the infinitesimal integration step
 * doB12max = 1:recalculate the light curve "plateau" brightness ; 0:don't
 *
 * NOTE: [definition of 1 unit distance]
 *  The physical distance between the stars at perihelion (minimal) is (1-e) units
 *  and the physical distance between the stars at aphelion (maximal) is (1+e) units
 *  so, 1 unit is the sum of the distances of the two semi-major axes
 *
 ******************************************************************/
static double flux (double time, double p[DIM], double dr, int doB12max)
{
  static double B12max ;
  const double omega = 2.0 * M_PI * (p[D_TPO] - p[D_TMO]) ;
  const double meanAnomaly = 2.0 * M_PI * (time - p[D_TPO] - p[D_TMO]) ;  // 2pi * (t - t0)
  double D2, D, E ;

  if (doB12max)
    B12max = (p[D_B1] * integrateWholeDisk1(p[D_R1])) + (p[D_B2] * integrateWholeDisk2(p[D_R2])) ;

  E = eccentricAnomaly(meanAnomaly, p[D_ECC]) ;

  D2 = sqr(1.0 - (p[D_ECC] * cos(E))) - sqr((((cos(E) - p[D_ECC]) * sin(omega)) + (sqrt(1.0 - sqr(p[D_ECC])) * sin(E) * cos(omega))) * p[D_SIN_I]) ;

  if (D2 > (EPSILON * EPSILON))  // Prevents sqrt of negative or division by zero, later on
    D = sqrt (D2) ;
  else
    D = EPSILON ;

  if (D >= (p[D_R1] + p[D_R2]))   // Side by side (A)
    return (B12max) ;

  if (sin(E + omega) > 0.0)  // Is star 2 in front ?
    return ((p[D_B2] * integrateWholeDisk2(p[D_R2])) + (p[D_B1] * integrateBackOverlap1(p[D_R2], p[D_R1], D, dr))) ;

  // star 1 is in front
  return ((p[D_B1] * integrateWholeDisk1(p[D_R1])) + (p[D_B2] * integrateBackOverlap2(p[D_R1], p[D_R2], D, dr))) ;
}


// same as flux(), only is assumes a circular orbit  (e=0, omega=0)
static double flux_noEcc (double time, double p[DIM], double dr, int doB12max)
{
 static double B12max ;
 double D2, D, sinEplusOmega ;
 
  if (doB12max)
    B12max = (p[D_B1] * integrateWholeDisk1(p[D_R1])) + (p[D_B2] * integrateWholeDisk2(p[D_R2])) ;
  
  sinEplusOmega = sin(2.0 * M_PI * (time - p[D_TMO])) ;
  D2 = 1.0 - sqr(sinEplusOmega * p[D_SIN_I]) ;

  if (D2 > (EPSILON * EPSILON))  // Prevents sqrt of negative or division by zero, later on
    D = sqrt (D2) ;
  else
    D = EPSILON ;

  if (D >= (p[D_R1] + p[D_R2]))   // Side by side (A)
    return (B12max) ;

  if (sinEplusOmega > 0.0)  // Is star 2 in front ?
    return ((p[D_B2] * integrateWholeDisk2(p[D_R2])) + (p[D_B1] * integrateBackOverlap1(p[D_R2], p[D_R1], D, dr))) ;

  // star 1 is in front
  return ((p[D_B1] * integrateWholeDisk1(p[D_R1])) + (p[D_B2] * integrateBackOverlap2(p[D_R1], p[D_R2], D, dr))) ;
}


//------------------------------------------------------


// Returns the reduced chi squared of the model (a score for how good the fit is)
// If the parameter vector is out of range, returns MAXFLOAT
double getChi2 (float *time, float *amp, float *err, int size,
		double e, double r1, double r2, double mag1, double mag2,
		double sin_i, double time0, double omega,
		double quadA1, double quadB1, double quadA2, double quadB2)
{
  int i ;
  double score ;
  double p[DIM] ;

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

  // Note that I use MAXFLOAT in a double since is may need to be increased by annealingFluctuation()
  if ((p[D_SIN_I] <= 0.0) || (p[D_SIN_I] > 1.0)) return (MAXFLOAT) ;
  if ((p[D_ECC] < 0.0) || (p[D_ECC] >= 1.0)) return (MAXFLOAT) ; 
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= (1.0 - p[D_ECC]))) return (MAXFLOAT) ;
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_ECC] - p[D_R1]))) return (MAXFLOAT) ;
  if (p[D_B1] <= 0.0) return (MAXFLOAT) ;
  if (p[D_B2] <= 0.0) return (MAXFLOAT) ;

  score = sqr((amp[0] - flux(time[0], p, MIN_INTEGRATION_STEP, 1)) / err[0]) ;

  for (i = 1 ; i < size ; i++)
    score += sqr((amp[i] - flux(time[i], p, MIN_INTEGRATION_STEP, 0)) / err[i]) ;
  
  return (score / size) ;
}


// Returns the reduced chi squared of the model (a score for how good the fit is)
// If the parameter vector is out of range, returns MAXFLOAT
// Assumes a circular orbit (e=0, omega=0)
double getChi2noEcc (float *time, float *amp, float *err, int size,
		     double r1, double r2, double mag1, double mag2, double sin_i, double time0,
		     double quadA1, double quadB1, double quadA2, double quadB2)
{
  int i ;
  double score ;
  double p[DIM] ;

  setLimbDarkening (quadA1, quadB1, quadA2, quadB2) ;
  
  p[D_R1] = r1 ;
  p[D_R2] = r2 ;
  p[D_B1] = pow(10.0, -0.4 * mag1) / integrateWholeDisk1(r1) ;
  p[D_B2] = pow(10.0, -0.4 * mag2) / integrateWholeDisk2(r2) ;
  p[D_SIN_I] = sin_i ;
  p[D_TMO] = time0 ;

  // Note that I use MAXFLOAT in a double since is may need to be increased by annealingFluctuation()
  if ((p[D_SIN_I] <= 0.0) || (p[D_SIN_I] > 1.0)) return (MAXFLOAT) ;
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= 1.0)) return (MAXFLOAT) ; 
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_R1]))) return (MAXFLOAT) ;   
  if (p[D_B1] <= 0.0) return (MAXFLOAT) ;  
  if (p[D_B2] <= 0.0) return (MAXFLOAT) ;    

  score = sqr((amp[0] - flux_noEcc(time[0], p, MIN_INTEGRATION_STEP, 1)) / err[0]) ;

  for (i = 1 ; i < size ; i++)
    score += sqr((amp[i] - flux_noEcc(time[i], p, MIN_INTEGRATION_STEP, 0)) / err[i]) ;
  
  return (score / size) ;
}


//-----------

// Writes the observed and model light curve
void writeLC (FILE *foutFit, FILE *foutDat, float *time, float *amp, float *err, int size,
	      double e, double r1, double r2, double mag1, double mag2,
	      double sin_i, double time0, double omega,
	      double quadA1, double quadB1, double quadA2, double quadB2)
{
  int i ;
  double p[DIM] ;

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

  // Should not happen:
  if ((p[D_SIN_I] <= 0.0) || (p[D_SIN_I] > 1.0)) return ;
  if ((p[D_ECC] < 0.0) || (p[D_ECC] >= 1.0)) return ; 
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= (1.0 - p[D_ECC]))) return ;
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_ECC] - p[D_R1]))) return ;
  if (p[D_B1] <= 0.0) return ;
  if (p[D_B2] <= 0.0) return ;

  fprintf (foutFit, "%f  %f\n", 0.0, -2.5 * log10(flux(0.0, p, MIN_INTEGRATION_STEP, 1))) ;
  
  for (i = 1 ; i < WRITE_CURVE_FIT ; i++)
    fprintf (foutFit, "%f  %f\n", (float)i / WRITE_CURVE_FIT, 
	     -2.5 * log10(flux((float)i / WRITE_CURVE_FIT, p, MIN_INTEGRATION_STEP, 0))) ;
  
  for (i = 0 ; i < size ; i++)
    fprintf (foutDat, "%f  %f  %f  %f\n", time[i],
	     -2.5 * log10(amp[i]), MAGNITUDE_ERROR_FACTOR * err[i] / amp[i],
	     2.5 * log10(flux(time[i], p, MIN_INTEGRATION_STEP, 0) / amp[i])) ;
}


// Writes the observed and model light curve
// Assumes a circular orbit (e=0, omega=0)
void writeLC_noEcc (FILE *foutFit, FILE *foutDat, float *time, float *amp, float *err, int size,
		    double r1, double r2, double mag1, double mag2, double sin_i, double time0,
		    double quadA1, double quadB1, double quadA2, double quadB2)
{
  int i ;
  double p[DIM] ;

  setLimbDarkening (quadA1, quadB1, quadA2, quadB2) ;
  
  p[D_R1] = r1 ;
  p[D_R2] = r2 ;
  p[D_B1] = pow(10.0, -0.4 * mag1) / integrateWholeDisk1(r1) ;
  p[D_B2] = pow(10.0, -0.4 * mag2) / integrateWholeDisk2(r2) ;
  p[D_SIN_I] = sin_i ;
  p[D_TMO] = time0 ;

  // Should not happen:
  if ((p[D_SIN_I] <= 0.0) || (p[D_SIN_I] > 1.0)) return ;
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= 1.0)) return ; 
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_R1]))) return ;   
  if (p[D_B1] <= 0.0) return ;  
  if (p[D_B2] <= 0.0) return ;    


  fprintf (foutFit, "%f  %f\n", 0.0, -2.5 * log10(flux_noEcc(time[0], p, MIN_INTEGRATION_STEP, 1))) ;
  
  for (i = 1 ; i < WRITE_CURVE_FIT ; i++)
    fprintf (foutFit, "%f  %f\n", (float)i / WRITE_CURVE_FIT, -2.5 * log10(flux_noEcc(time[i], p, MIN_INTEGRATION_STEP, 0))) ;
  
  for (i = 0 ; i < size ; i++)
    fprintf (foutDat, "%f  %f  %f  %f\n", time[i], -2.5 * log10(amp[i]), MAGNITUDE_ERROR_FACTOR * err[i] / amp[i],
	     2.5 * log10(flux_noEcc(time[i], p, MIN_INTEGRATION_STEP, 0) / amp[i])) ;
}


// ---------------------------------------


// Finds the phase, and depth of the eclipses
// Note that the primary eclipse (pphase1, pdmag1) occures when star1 (r1, mag1) is being eclipsed (i.e. is behind),
// and that the secondary eclipse (pphase2, pdmag2) occures when star2 (r2, mag2) is being eclipsed
// Returns 0 = success ; otherwise failure
int findEclipses (double e, double r1, double r2, double mag1, double mag2,
		  double sin_i, double time0, double omega,
		  double quadA1, double quadB1, double quadA2, double quadB2,
		  double *pphase1, double *pdmag1, double *pphase2, double *pdmag2)
{
  double prevFlux, currFlux, nextFlux, plateauFlux, plateauFluxEpsilon ;
  double fluxEclipse1 = 0.0, fluxEclipse2 = 0.0, sumPhase1 = 0.0, sumPhase2 = 0.0 ;
  double E, time, p[DIM] ;
  unsigned long countEclipse1 = 0, countEclipse2 = 0 ;
  int isAfterEclipse1 = 0, isAfterEclipse2 = 0 ;

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

  // Note that I use MAXFLOAT in a double since is may need to be increased by annealingFluctuation()
  if ((p[D_SIN_I] <= 0.0) || (p[D_SIN_I] > 1.0)) return (1) ;
  if ((p[D_ECC] < 0.0) || (p[D_ECC] >= 1.0)) return (2) ; 
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= (1.0 - p[D_ECC]))) return (3) ;
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_ECC] - p[D_R1]))) return (4) ;
  if (p[D_B1] <= 0.0) return (5) ;
  if (p[D_B2] <= 0.0) return (6) ;

  plateauFlux = pow(10.0, -0.4 * mag1) + pow(10.0, -0.4 * mag2) ;
  plateauFluxEpsilon = plateauFlux * FLUX_EPSILON ;

  prevFlux = flux(1.0 - ECLIPSE_SEARCH_STEP, p, MIN_INTEGRATION_STEP, 1) ;
  currFlux = flux(0.0, p, MIN_INTEGRATION_STEP, 0) ;

  for (time = ECLIPSE_SEARCH_STEP ; time < 1.0 - ECLIPSE_SEARCH_STEP ; time += ECLIPSE_SEARCH_STEP)
    {
      nextFlux = flux(time, p, MIN_INTEGRATION_STEP, 0) ;

      if ((plateauFlux - currFlux > plateauFluxEpsilon) && (prevFlux >= currFlux) && (currFlux <= nextFlux))
	{
	  E = eccentricAnomaly(2.0 * M_PI * (time - time0), e) ;

	  if (sin(E + (2.0 * M_PI * omega)) > 0.0)  // star1 is being eclipsed
	    {
	      if (isAfterEclipse1)
		sumPhase1 += (time - 1.0) ;
	      else
		sumPhase1 += time ;

	      countEclipse1++ ;
	      fluxEclipse1 = currFlux ;
	    }
	  else // star2 is being eclipsed
	    {
	      if (isAfterEclipse2)
		sumPhase2 += (time - 1.0) ;
	      else
		sumPhase2 += time ;
	      
	      countEclipse2++ ;
	      fluxEclipse2 = currFlux ;
	    }
	}
      else
	{
	  if (countEclipse1)
	    isAfterEclipse1 = 1 ;

	  if (countEclipse2)
	    isAfterEclipse2 = 1 ;
	}

      prevFlux = currFlux ;
      currFlux = nextFlux ;
    }

  
  if ((countEclipse1 == 0) || (countEclipse2 == 0))  // missing eclipse
    return (-1) ;

  *pphase1 = mod1 (sumPhase1 / countEclipse1) ;
  *pphase2 = mod1 (sumPhase2 / countEclipse2) ;
  
  if ((fluxEclipse1 <= 0.0) || (fluxEclipse2 <= 0.0))  // invalid flux
    return (-2) ;

  *pdmag1 = 2.5 * log10 (plateauFlux / fluxEclipse1) ;
  *pdmag2 = 2.5 * log10 (plateauFlux / fluxEclipse2) ;

  return (0) ;
}


// Finds the phase, and depth of the eclipses
// Assumes a circular orbit (e=0, omega=0)
// Note that the primary eclipse (pphase1, pdmag1) occures when star1 (r1, mag1) is being eclipsed (i.e. is behind),
// and that the secondary eclipse (pphase2, pdmag2) occures when star2 (r2, mag2) is being eclipsed
// Returns 0 = success ; otherwise failure
int findEclipses_noEcc (double r1, double r2, double mag1, double mag2, double sin_i, double time0,
			double quadA1, double quadB1, double quadA2, double quadB2,
			double *pphase1, double *pdmag1, double *pphase2, double *pdmag2)
{
  double p[DIM], plateauFlux, fluxEclipse1 = 0.0, fluxEclipse2 = 0.0 ;

  setLimbDarkening (quadA1, quadB1, quadA2, quadB2) ;
  
  p[D_R1] = r1 ;
  p[D_R2] = r2 ;
  p[D_B1] = pow(10.0, -0.4 * mag1) / integrateWholeDisk1(r1) ;
  p[D_B2] = pow(10.0, -0.4 * mag2) / integrateWholeDisk2(r2) ;
  p[D_SIN_I] = sin_i ;
  p[D_TMO] = time0 ;

  // Note that I use MAXFLOAT in a double since is may need to be increased by annealingFluctuation()
  if ((p[D_SIN_I] <= 0.0) || (p[D_SIN_I] > 1.0)) return (1) ;
  if ((p[D_R1] <= 0.0) || (p[D_R1] >= 1.0)) return (2) ; 
  if ((p[D_R2] <= 0.0) || (p[D_R2] >= (1.0 - p[D_R1]))) return (3) ;   
  if (p[D_B1] <= 0.0) return (4) ;  
  if (p[D_B2] <= 0.0) return (5) ;    

  plateauFlux = pow(10.0, -0.4 * mag1) + pow(10.0, -0.4 * mag2) ;

  *pphase1 = mod1 (time0 + 0.25) ;
  *pphase2 = mod1 (time0 + 0.75) ; 

  fluxEclipse1 = flux_noEcc(*pphase1, p, MIN_INTEGRATION_STEP, 1) ;
  fluxEclipse2 = flux_noEcc(*pphase2, p, MIN_INTEGRATION_STEP, 0) ;

 if ((fluxEclipse1 <= 0.0) || (fluxEclipse2 <= 0.0))  // invalid flux
    return (-1) ;

 *pdmag1 = 2.5 * log10 (plateauFlux / fluxEclipse1) ;
 *pdmag2 = 2.5 * log10 (plateauFlux / fluxEclipse2) ;

  return (0) ;
}
