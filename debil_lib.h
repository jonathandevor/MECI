/************************************************************************
 * DEBiL library of assorted light curve utilities
 *
 * Written by Jonathan Devor (jdevor@cfa.harvard.edu)
 *    Harvard-Smithsonian Center for Astrophysics 
 *    60 Garden St., Cambridge, MA 02138 USA
 *
 ***********************************************************************/

// Prints a given error message (errStr) reported by a given function (functionName) about a
// given light curve file (dataFilename). It is added to a given error repository file (errFilename)
// isFatal: 1 = Error - go to next  ;  0 = Warning - try anyways  ; -1 = Note - just informative
void printError (char *errFilename, int typeID, char *dataFilename, char *functionName, char *errStr) ;


// Allocates and sets: *ptime, *pamp, *perr
// Returns the number of light curve measurements
// Returns a negative number upon error, and does not perform allocation
int loadLC (char *filename, double period, char *errfilename, float **ptime, float **pamp, float **perr) ;


// Same as loadLC(), but it also returns the midpoint of the timeStamps, without modulo (i.e. an epoch)
// Allocates and sets: *ptime, *pamp, *perr
// Returns the number of light curve measurements
// Returns a negative number upon error, and does not perform allocation
int loadLC_midTime (char *filename, double period, char *errfilename, float **ptime, float **pamp, float **perr, double *pmidTime) ;


// Returns 0 = success ; 1 = failure (malloc)
int sortSamples (float *time, float *amp, float *err, int size) ;


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
// Returns the chi2 around the spline (MAXFLOAT in case of error)
double ridOutliers (float *time, float *amp, float *err, int *psize, float varianceLimit) ;


// Returns the reduced chi square of a given model, with given light curve data
double getChi2 (float *time, float *amp, float *err, int size,
		double e, double r1, double r2, double mag1, double mag2,
		double sin_i, double time0, double omega,
		double quadA1, double quadB1, double quadA2, double quadB2) ; 


// Returns the reduced chi square of a given model, with given light curve data
// Assumes a circular orbit (e=0, omega=0)
double getChi2noEcc (float *time, float *amp, float *err, int size,
		     double r1, double r2, double mag1, double mag2, double sin_i, double time0,
		     double quadA1, double quadB1, double quadA2, double quadB2) ;


// Writes the observed and model light curve
void writeLC (FILE *foutFit, FILE *foutDat, float *time, float *amp, float *err, int size,
	      double e, double r1, double r2, double mag1, double mag2,
	      double sin_i, double time0, double omega,
	      double quadA1, double quadB1, double quadA2, double quadB2) ;


// Writes the observed and model light curve
// Assumes a circular orbit (e=0, omega=0)
void writeLC_noEcc (FILE *foutFit, FILE *foutDat, float *time, float *amp, float *err, int size,
		    double r1, double r2, double mag1, double mag2, double sin_i, double time0,
		    double quadA1, double quadB1, double quadA2, double quadB2) ;


// Finds the phase, and depth of the eclipses
// Note that the primary eclipse (pphase1, pdmag1) occures when star1 (r1, mag1) is being eclipsed (i.e. is behind),
// and that the secondary eclipse (pphase2, pdmag2) occures when star2 (r2, mag2) is being eclipsed
// Returns 0 = success ; otherwise failure
int findEclipses (double e, double r1, double r2, double mag1, double mag2,
		  double sin_i, double time0, double omega,
		  double quadA1, double quadB1, double quadA2, double quadB2,
		  double *pphase1, double *pdmag1, double *pphase2, double *pdmag2) ;


// Finds the phase, and depth of the eclipses
// Assumes a circular orbit (e=0, omega=0)
// Note that the primary eclipse (pphase1, pdmag1) occures when star1 (r1, mag1) is being eclipsed (i.e. is behind),
// and that the secondary eclipse (pphase2, pdmag2) occures when star2 (r2, mag2) is being eclipsed
// Returns 0 = success ; otherwise failure
int findEclipses_noEcc (double r1, double r2, double mag1, double mag2, double sin_i, double time0,
			double quadA1, double quadB1, double quadA2, double quadB2,
			double *pphase1, double *pdmag1, double *pphase2, double *pdmag2) ;

