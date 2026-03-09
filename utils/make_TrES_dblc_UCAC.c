/*******************************************************************
 * This program find the colors and RA/DEC (from 2MASS), and proper motion (from UCAC 2.4) of the each LC in a 31-param DEBiL DB.
 * It may also be used to filters through only the light curves that are red  and have proper motions (using UCAC 2.4).
 * Note that this program shouldn't be run in parallel (unless you set different SCAT_TMP_FILENAMEs).
 *
 * compile: gcc make_TrES_dblc_UCAC.c -o make_TrES_dblc_UCAC -Wall -lm -O
 *          gcc make_TrES_dblc_UCAC.c -o filter_TrES_dbl -Wall -lm -O        [after remarking "#define DO_DBLC"]
 *          gcc make_TrES_dblc_UCAC.c -o make_TrES_dblc_UCAC_pm -Wall -lm -O      [with "#define FILTER_OUT_ZERO_PM"]
 *          gcc make_TrES_dblc_UCAC.c -o filter_TrES_dbl_pm -Wall -lm -O      [remarking "#define DO_DBLC" and with "#define FILTER_OUT_ZERO_PM"]
 *
 * run:     make_TrES_dblc_UCAC listAll.edit.dbl and0.2mass listAll.edit.dblc
 *
 *******************************************************************/


#include <math.h>
#include <stdio.h>
#include <string.h>  // strlen(), memcpy()
#include <stdlib.h>  // atol()

#define USE_UCAC_2MASS_COLORS  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#define ERROR_COLOR 0.05

#define ERROR_FILTER (ERROR_COLOR * 0.707106781)   // sqrt(0.5)

#define MAX_STRLEN 256
#define TrES_NUMBER_LENGTH 5   // TrES number inside of filename

#define SCAT_TMP "make_TrES_dblc_tmp"

//#define MAX_RA_DEC_ERR 1.0  // in units of arcsecond
//#define MAX_UCAC_MAG_DIFF 1.0 // in magnitudes
//#define MAX_UCAC_TrES_MAG_DIFF 1.5 // in magnitudes
  
//#define LINUX   // otherwise Solaris
#define DO_DBLC   // otherwise just the 31-param DEBiL format  (for filtering)
#define DO_DBLC_WITH_ABS_MAG

//#define FILTER_OUT_BLUE_ZERO_PM   // outputs only the red binaries with detectable proper motion
//#define MIN_JmH 0.13   // spectral type F0
//#define MIN_HmK 0.03   // spectral type F0


// converts (J-H)_{2mass} into (J-H)_{Bessel & Brett}  using eq. A2 of 
// "Color transformations for the 2mass second incremental data release" by J. Carpenter, 2001
double convert_2mass_JmH (double JmH)
{
  return ((JmH + 0.045) / 0.980) ;
}


// converts (H-Ks)_{2mass} into (H-K)_{Bessel & Brett}  using eq. A4 of 
// "Color transformations for the 2mass second incremental data release" by J. Carpenter, 2001
double convert_2mass_HmK (double HmK)
{
  return ((HmK - 0.028) / 0.996) ;
}


// converts Ks_{2mass} into K_{Bessel & Brett}  using eq. A1 of 
// "Color transformations for the 2mass second incremental data release" by J. Carpenter, 2001
double convert_2mass_K (double Ks)
{
  return (Ks + 0.044) ;
}


void finishLine (FILE *fin)
{
  while (!feof(fin) && (fgetc(fin) != '\n')) ;
}


// returns -1 upon error, otherwise the number of the TrES file
long extractNum (char *nameStr)
{
  char numStr[TrES_NUMBER_LENGTH+1] ;
  int foundNonDigit, i, j, len = strlen(nameStr) ;

  for (i = len - TrES_NUMBER_LENGTH ; i > 0 ; i--)
    {
      foundNonDigit = 0 ;
      
      for (j = 0 ; (!foundNonDigit) && (j < TrES_NUMBER_LENGTH) ; j++)
	if ((nameStr[i+j] < '0') || (nameStr[i+j] > '9'))
	  foundNonDigit = 1 ;

      if (!foundNonDigit)
	{
	  memcpy (numStr, &nameStr[i], TrES_NUMBER_LENGTH * sizeof(char)) ;
	  numStr[TrES_NUMBER_LENGTH] = 0 ;
	  return (atol(numStr)) ;
	}
    }

  return (-1) ;
}


int main(int argc, char **argv)
{ 
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20;
  float pmRA, pmDEC, ra_dec_err ;
  char x1[MAX_STRLEN], str[MAX_STRLEN], strRA[MAX_STRLEN], strDEC[MAX_STRLEN], strScatFilename[MAX_STRLEN] ;
  FILE *finDebil, *finData, *fout, *finUCAC ;
  long nameNum = -1, nameNum2mass, prevNameNum ;
  double magJ, magH, magK, magJ_UCAC, magH_UCAC, magK_UCAC ;  //  all 2MASS

  if (argc != 4)
    {
      printf ("usage %s  <debil file> <2mass data file> <output file>\n", argv[0]) ;
      return (1) ;
    }

  finDebil = fopen (argv[1], "rt") ;
  finData = fopen (argv[2], "rt") ;
  fout = fopen (argv[3], "wt") ;

  if (!finDebil || !finData || !fout)
    {
      printf ("ERROR: couldn't open a file (%d %d %d)\n", !finDebil, !finData, !fout) ;
      return (2) ;
    }

  sprintf (str, "rm -i %s.%s.*", argv[3], SCAT_TMP) ;
  system (str) ;
  printf ("There should be a 'No such file or directory' error above. Please ignore it.\n") ;

  while (31 == fscanf(finDebil, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		      &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {
      printf ("scanning matches for '%s'", x1) ;

      prevNameNum = nameNum ;
      nameNum = extractNum (x1) ;

      printf (" (%lu) -- ", nameNum) ;

      if (nameNum < 0)
	{
	  printf ("ERROR in file name ('%s') - skipping", x1) ;
	  continue ;
	}

      //-------------- look for 2mass color -------------------------------

      if ((prevNameNum >= nameNum) || (prevNameNum < 0))
	{
	  rewind (finData) ;
	  finishLine (finData) ;  // gets rid of header line
	}

      while (6 == fscanf (finData, "%ld %*s  %*s  %s  %s  %lf %lf %lf", &nameNum2mass, strRA, strDEC, &magJ, &magH, &magK))
	{
	  finishLine (finData) ;
	  
	  if (nameNum == nameNum2mass)
	    break ;
	}
            
      if (nameNum != nameNum2mass)
	{
	  printf ("Warning: Couldn't find 2mass color\n") ;
	  continue ;
	}

#ifdef FILTER_OUT_BLUE_ZERO_PM
      if ((magJ - magH < MIN_JmH) || (magH - magK < MIN_HmK))
	{
	  printf ("Warning: One of the color is too blue (J-H=%f   H-K=%f)\n", magJ - magH, magH - magK) ;
	  continue ;
	}
#endif

      //------------------ look for UCAC 2.4 proper motion ------------------

      sprintf (strScatFilename, "%s.%s.%ld", argv[3], SCAT_TMP, nameNum) ;

#ifdef LINUX
      sprintf (str,"/data/oir/wcs/bin.redhat/scat -atc ucac2.4 %s %s J2000 | tail -1 > %s", strRA, strDEC, strScatFilename) ;
#else
      sprintf (str,"/data/oir/wcs/bin.solaris/scat -atc ucac2.4 %s %s J2000 | tail -1 > %s", strRA, strDEC, strScatFilename) ;
#endif
      
      system (str) ;
      
      finUCAC = fopen (strScatFilename, "rt") ;
      if (!finUCAC)
	{
	  printf ("Warnnig: Couln't open the UCAC file ('%s')\n", strScatFilename) ;
	  continue ;
	}

      if (6 != fscanf (finUCAC, "%*f %*d:%*d:%*f  %*d:%*d:%*f   %lf %lf %lf %*f  %f %f %f", &magJ_UCAC, &magH_UCAC, &magK_UCAC, &pmRA, &pmDEC, &ra_dec_err))
	{
	  printf ("Warning: The UCAC map doesn't extend to this location, or invalid UCAC parsing ('%s')\n", strScatFilename) ;
	  fclose (finUCAC) ;
	  sprintf (str, "rm %s", strScatFilename) ;
	  system (str) ;
	  continue ;
	}

      fclose (finUCAC) ;
      
      sprintf (str, "rm %s", strScatFilename) ;
      system (str) ;

#ifdef USE_UCAC_2MASS_COLORS
      magJ = magJ_UCAC ;
      magH = magH_UCAC ;
      magK = magK_UCAC ;
#endif

      /*
	if (ra_dec_err <= MAX_RA_DEC_ERR)
      	{
	printf ("Warning: UCAC closest star is too far (%f)\n", ra_dec_err) ;
	continue ;
      	}
      */
      
#ifdef FILTER_OUT_BLUE_ZERO_PM
      if ((pmRA == 0.0) && (pmDEC == 0.0))
	{
	  printf ("No detected PM\n") ;
	  continue ;
	}
#endif
      
      printf ("OK\n") ;

      //-----------------------------------
      
      fprintf (fout, "%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
	       x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
	       x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
	  	  
#ifdef DO_DBLC
#ifdef DO_DBLC_WITH_ABS_MAG
      fprintf (fout, "R  5   %s  %s   %.3f  %.3f  %.3f\n", strRA, strDEC, pmRA, pmDEC, ra_dec_err) ;
      fprintf (fout, "J-H  %f  %f\n", convert_2mass_JmH(magJ - magH), ERROR_COLOR) ;
      fprintf (fout, "H-K  %f  %f\n", convert_2mass_HmK(magH - magK), ERROR_COLOR) ;
      fprintf (fout, " J   %f  %f\n", convert_2mass_JmH(magJ - magH) + convert_2mass_HmK(magH - magK) + convert_2mass_K(magK), ERROR_FILTER) ;
      fprintf (fout, " H   %f  %f\n", convert_2mass_HmK(magH - magK) + convert_2mass_K(magK), ERROR_FILTER) ;
      fprintf (fout, " K   %f  %f\n", convert_2mass_K(magK), ERROR_FILTER) ;
#else
      fprintf (fout, "R  2   %s  %s   %.3f  %.3f  %.3f\n", strRA, strDEC, pmRA, pmDEC, ra_dec_err) ;
      fprintf (fout, "J-H  %f  %f\n", convert_2mass_JmH (magJ - magH), ERROR_COLOR) ;
      fprintf (fout, "H-K  %f  %f\n", convert_2mass_HmK (magH - magK), ERROR_COLOR) ;
#endif
#endif
    }	     

  fclose (finData) ;
  fclose (finDebil) ;
  fclose (fout) ;
  return (0) ;
}
