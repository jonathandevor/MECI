/*******************************
 * Find the best fit of binary system to a set of given observables
 * Configured for DEBiL (31-param), over OGLE (I-band)
 *
 * compile:  gcc MECIexp.c -o MECIexp  -lm -Wall -O4
 *
 * run:  MECIexp ag.dbl
 *
 * if using contour:
 *
 *       cp /home/jdevor/work/MECIexp/IDL_tt.template .
 *
 *       idl
 *       IDL> .run ag.dbl.contour_idl.scr
 *       IDL> exit
 *
 ********************************/

#include <math.h>
#include <values.h>      // for MAXFLOAT
#include <stdio.h>
#include "starTypes.h"

#define FILTER V   // choose from:  B, V, R, I, K, J 

#define PAIR_ARR_SIZE 2  // shouldn't be larger than (NUM_STAR_TYPES^2)

#define DO_CONTOURS


struct pairStar
{
  int star1, star2 ;  // 0 .. (NUM_STAR_TYPES-1)
  double sigma ;
} pairArr[PAIR_ARR_SIZE] ;


double sqr (double x)
{
  return (x * x) ;
}

double cube (double x)
{
  return (x * x * x) ;
}


// performs one iteration of an insertion sort
void updatePairArr (struct pairStar pair)
{
  int index, j ;

  if (pair.sigma > pairArr[PAIR_ARR_SIZE-1].sigma)
    return ;

  index = PAIR_ARR_SIZE-2 ;
  while ((index >= 0) && (pair.sigma < pairArr[index].sigma))
    index-- ;

  index++ ;

  for (j = (PAIR_ARR_SIZE-1) ; j > index ; j--)
    {
      pairArr[j].star1 = pairArr[j-1].star1 ;
      pairArr[j].star2 = pairArr[j-1].star2 ;
      pairArr[j].sigma = pairArr[j-1].sigma ;
    }

  pairArr[index].star1 = pair.star1 ;
  pairArr[index].star2 = pair.star2 ;
  pairArr[index].sigma = pair.sigma ;
}



FILE* contourHeader (char *nameStr, const char *str, FILE *foutIDL)
{
  char filename[256] ;

  fprintf (foutIDL, "device, file=\"%s_%s.ps\"\n", nameStr, str) ;
  fprintf (foutIDL, "z=read_ascii(\"%s_%s.txt\", template=tt)\n", nameStr, str) ;
  fprintf (foutIDL, "contour, z.field3, z.field1, z.field2, /irregular, /follow") ;
  fprintf (foutIDL, ",xtitle='Star1 type', ytitle='Star2 type', title='%s  %s'", nameStr, str) ;
  fprintf (foutIDL, ",xtickname=['O5','B0','B5','A0','A5','F0','F5','G0','G5','K0','K5','M0','M2','M5'], xticks=13") ;
  fprintf (foutIDL, ",ytickname=['O5','B0','B5','A0','A5','F0','F5','G0','G5','K0','K5','M0','M2','M5'], yticks=13") ;
  fprintf (foutIDL, ",levels = [1,2,5,10,20,50,100,200,500,1000]") ;
  fprintf (foutIDL, ",c_annotation=['1','2','5','10','20','50','100','200','500','1000']") ;   
  fprintf (foutIDL, "\n") ;   
  fprintf (foutIDL, "device, /close\n\n") ;

  sprintf (filename, "%s_%s.txt", nameStr, str) ;

  //  printf ("%s\n", filename) ;

  return (fopen (filename, "wt")) ;
}



int main(int argc, char **argv)
{
  struct starInfo star[NUM_STAR_TYPES] ;
  struct pairStar pair ;
  int i, j ;
  double mag2mag1, mag2mag1err ;
  double tmag2mag1, t_r1, t_r2, inv_a ;
  double sigma2mag2mag1, sigma2r1, sigma2r2 ;
  //  double JmK, JmKerr, tJmK, BmV, BmVerr, tBmV ; 
  FILE *fin ;
  double x_period, x_r1, x_r1err, x_r2, x_r2err, x_mag1, x_mag1err, x_mag2, x_mag2err, x_chi2 ;
  char x_filename[256] ;

#ifdef DO_CONTOURS
  FILE *fout_mag2mag1, *fout_r1, *fout_r2, *fout_comb, *foutIDL ;
  char nameStr[256] ;
#endif 
  
  if (argc != 2)
    {
      printf (" Usage:  %s  <DEBiL database (31-col)>\n", argv[0]) ;
      return (1) ;
    }


  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (2) ;
    }  
  
  initialize (star) ;

#ifdef DO_CONTOURS
  sprintf (nameStr, "%s.contour_idl.scr", argv[1]) ;
  foutIDL = fopen (nameStr, "wt") ;
  
  if (!foutIDL)
    {
      fclose (fin) ;
      printf ("ERROR: couldn't open file '%s'\n", nameStr) ;
      return (3) ;
    }

  fprintf (foutIDL, "restore,'IDL_tt.template'\n") ;
  fprintf (foutIDL, "set_plot, 'ps'\n\n") ;
#endif

  //-----------------------------------------------------


  while (11 == fscanf (fin, "%s %lf %*f %*f %lf %lf %lf %lf %lf %lf %lf %lf %*f %*f %*f %*f %*f %*f %*d %*d %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f\n",
		       x_filename, &x_period, &x_r1, &x_r1err, &x_r2, &x_r2err, &x_mag1, &x_mag1err, &x_mag2, &x_mag2err, &x_chi2))
    {
      mag2mag1 = x_mag2 - x_mag1 ;
      mag2mag1err = sqrt(sqr(x_mag1err) + sqr(x_mag2err)) ;

      for (i = 0 ; i < PAIR_ARR_SIZE ; i++)
	{
	  pairArr[i].star1 = -1 ;
	  pairArr[i].star2 = -1 ;
	  pairArr[i].sigma = MAXFLOAT ;
	}

      //--------------------------------------------------------

#ifdef DO_CONTOURS
      for (i = 0, j = 0 ; x_filename[j] ; j++)  // removes path
	if (x_filename[j] == '/') 
	  i = j+1 ;
      
      sprintf (nameStr, "%s_%s.contour", argv[1], &x_filename[i]) ;

      //    fout_color = contourHeader (nameStr, "color", foutIDL) ;
      fout_mag2mag1 = contourHeader (nameStr, "mag", foutIDL) ;
      fout_r1 = contourHeader (nameStr, "r1", foutIDL) ;
      fout_r2 = contourHeader (nameStr, "r2", foutIDL) ;
      fout_comb = contourHeader (nameStr, "comb", foutIDL) ;
      
      if (!fout_mag2mag1 || !fout_r1 || !fout_r2 || !fout_comb)
	{
	  printf ("ERROR: couldn't open an output file\n") ;
	  fclose (fin) ;
	  //	  if (fout_color) fclose (fout_color) ;
	  if (fout_mag2mag1) fclose (fout_mag2mag1) ;
	  if (fout_r1) fclose (fout_r1) ;
	  if (fout_r2) fclose (fout_r2) ;
	  if (fout_comb) fclose (fout_comb) ;
	  return (4) ;
	}  
#endif


      //--------------------------------------------------------

      for (i = 0 ; i < NUM_STAR_TYPES ; i++)
	for (j = 0 ; j < NUM_STAR_TYPES ; j++)
	  {
	    // values of all possible star pairings:
	    /*
	      tJmK = -2.5 * (log10(pow(10.0, -0.4 * star[i].J) + pow(10.0, -0.4 * star[j].J)) -
	      log10(pow(10.0, -0.4 * star[i].K) + pow(10.0, -0.4 * star[j].K))) ;

	      tBmV = -2.5 * (log10(pow(10.0, -0.4 * star[i].B) + pow(10.0, -0.4 * star[j].B)) -
	      log10(pow(10.0, -0.4 * star[i].V) + pow(10.0, -0.4 * star[j].V))) ;
	    */

	    tmag2mag1 = star[j].FILTER - star[i].FILTER ;

	    // 0.266505 = (G/3pi * (day/sec)^2)^(-1/3)
	    inv_a = 0.266505 * pow(x_period * x_period * ((star[i].rho * cube(star[i].r)) + (star[j].rho * cube(star[j].r))), -1.0 / 3.0) ;
	    t_r1 = star[i].r * inv_a ;
	    t_r2 = star[j].r * inv_a ;

	    //--------------

	    //	    sigma2color = sqr((JmK - tJmK) / JmKerr) + sqr((BmV - tBmV) / BmVerr) ;
	    sigma2mag2mag1 = sqr((mag2mag1 - tmag2mag1) / mag2mag1err) ;
	    sigma2r1 = sqr((x_r1 - t_r1) / x_r1err) ;
	    sigma2r2 = sqr((x_r2 - t_r2) / x_r2err) ;


	    pair.star1 = i ;
	    pair.star2 = j ;
	    //	    pair.sigma = sqrt(sigma2color + sigma2mag2mag1 + sigma2r1 + sigma2r2) ;
	    pair.sigma = sqrt(sigma2mag2mag1 + sigma2r1 + sigma2r2) ;

	    updatePairArr(pair) ;


	    //--------------
#ifdef DO_CONTOURS
	    if ((i < NUM_MAIN_SEQUENCE) && (j < NUM_MAIN_SEQUENCE))
	      {
		//		fprintf (fout_color, "%d %d %f\n", i, j, sqrt(sigma2color)) ;
		fprintf (fout_mag2mag1, "%d %d %f\n", i, j, sqrt(sigma2mag2mag1)) ;
		fprintf (fout_r1, "%d %d %f\n", i, j, sqrt(sigma2r1)) ;
		fprintf (fout_r2, "%d %d %f\n", i, j, sqrt(sigma2r2)) ;
		fprintf (fout_comb, "%d %d %f\n", i, j, pair.sigma) ;
	      }
#endif
	  }

      printf ("%s  P=%.2f chi2=%.2f :   ", x_filename, x_period, x_chi2) ;

      for (i = 0 ; i < PAIR_ARR_SIZE ; i++)
	{
	  if ((pairArr[i].star1 == -1) || (pairArr[i].star2 == -1))
	    printf ("ERROR ") ;
	  else
	    printf ("%c%d;%c - %c%d;%c  %.2f  ", 
		    star[pairArr[i].star1].type, star[pairArr[i].star1].num, star[pairArr[i].star1].branch,
		    star[pairArr[i].star2].type, star[pairArr[i].star2].num, star[pairArr[i].star1].branch,
		    pairArr[i].sigma) ;
	}

      printf ("\n") ;

#ifdef DO_CONTOURS
      //    fclose (fout_color) ;
      fclose (fout_mag2mag1) ;
      fclose (fout_r1) ;
      fclose (fout_r2) ;
      fclose (fout_comb) ;
#endif
    }


  //-------------------------------------------------------------------------------

#ifdef DO_CONTOURS
  fprintf (foutIDL, "end\n") ;
  fclose (foutIDL) ;
#endif

  fclose (fin) ;
  return (0) ;
}
