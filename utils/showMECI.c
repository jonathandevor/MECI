/**********************************************************
 *
 * This program makes a catalogue of the MECI best-fits
 *
 * compile:  gcc showMECI.c -o showMECI -Wall -lm -O4
 *
 * run:      ls *.meci_fit > meci_fit.list
 *           showMECI meci_fit.list meci_fit.sm
 *           
 *           sm
 *           : inp meci_fit.sm
 *           : quit
 *
 *           ps2pdf meci_fit.ps
 *           acroread meci_fit.pdf &
 *
 **********************************************************/


#include <stdio.h>
#include <string.h>   // for strcmp() and strlen()


#define NUM_BOXES_X 4  // for SM file
#define NUM_BOXES_Y 6


//-------------------------------------------------------


void printSM (FILE *foutSM, char *filenameFit, char *filenameDat)
{
  static int i = 0, j = NUM_BOXES_Y ;
  int indexStart, indexEnd = strlen (filenameFit) - 1 ;

  i++ ;

  if (i > NUM_BOXES_X)
    {
      i = 1 ;
      j -= 2 ;
    }

  if (j <= 0)
    {
      j = NUM_BOXES_Y ;
      fprintf (foutSM, "page\n\n") ;
    }

  
  fprintf (foutSM, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j) ;
  
  fprintf (foutSM, "data %s\n", filenameFit) ;
  fprintf (foutSM, "read {phaseFit 1 meci 2}\n") ;
  fprintf (foutSM, "data %s\n", filenameDat) ;
  fprintf (foutSM, "read {phaseDat 1 obs 2 res 4}\n") ;
  
  //------
  while ((indexEnd >= 0) && (filenameFit[indexEnd] != '.'))
    indexEnd-- ;
  
  if (indexEnd >= 0)
    filenameFit[indexEnd] = 0 ;
  
  indexStart = indexEnd - 1 ;
  while ((indexStart >= 0) && (filenameFit[indexStart] != '/'))
    indexStart-- ;
  
  if (indexStart >= 0)
    memmove (filenameFit, &filenameFit[indexStart+1], indexEnd - indexStart) ;
  //------

  fprintf (foutSM, "limits 0 1 obs\n") ;
  fprintf (foutSM, "limits 0 1 $fy2 $fy1\n") ; // invert axis direction
  // add margin:	   limits 0 1 $($fy2+0.2) $($fy1-0.1)

  fprintf (foutSM, "ctype blue\n") ;
  fprintf (foutSM, "connect phaseFit meci\n") ;
  fprintf (foutSM, "ctype black\n") ;
  fprintf (foutSM, "points phaseDat obs\n") ;
  
  fprintf (foutSM, "ticksize 0 0 0 0\n") ;
  fprintf (foutSM, "xlabel %s\n", filenameFit) ;
  fprintf (foutSM, "box\n\n") ;
  //------------------------------
  
  fprintf (foutSM, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j-1) ;
  fprintf (foutSM, "limits 0 1 0.3 -0.3\n") ; // invert axis direction
  fprintf (foutSM, "points phaseDat res\n") ;
  
  fprintf (foutSM, "ticksize 0.05 0.2 0.05 0.1\n") ;
  fprintf (foutSM, "xlabel %s\n", filenameFit) ;
  fprintf (foutSM, "ylabel Residuals\n") ;
  fprintf (foutSM, "box\n\n") ;

  fprintf (foutSM, "set phaseFit = 0\n") ;
  fprintf (foutSM, "set phaseDat = 0\n") ;
  fprintf (foutSM, "set obs = 0\n") ;
  fprintf (foutSM, "set meci = 0\n") ;
  fprintf (foutSM, "set res = 0\n\n") ;
}



int main (int argc, char **argv)
{
  FILE *fin, *foutSM ;
  char filenameFit[256], filenameDat[256] ;

  if (argc != 3)
    {
      printf ("usage: %s <MECI fit list> <output SM file>\n", argv[0]) ;
      printf ("\n\nYou can create a <MECI fit list> by running, for example:\n") ;
      printf ("ls *.meci_fit > meci_fit.list\n") ;
      return (1) ;
    }

  if (!strcmp(argv[1], argv[2]))
    {
      printf ("ERROR: Input and output filenames must be different.\n") ;
      return (2) ;
    }

  if (!(fin = fopen (argv[1], "rt")))
    {
      printf ("ERROR: Couldn't open the MECI fit list ('%s').\n", argv[1]) ;
      return (3) ;
    }

  if (!(foutSM = fopen (argv[2], "wt")))
    {
      printf ("ERROR: Couldn't open the output SM file ('%s').\n", argv[2]) ;
      return (4) ;
    }

  //-------------------------------------

  fprintf (foutSM, "device postportfile %s.ps\n\n", argv[2]) ;
  fprintf (foutSM, "expand 0.4\n") ;
  fprintf (foutSM, "ptype 1 1\n\n") ;

  while (1 == fscanf (fin, "%s\n", filenameFit))
    { 
      strcpy (filenameDat, filenameFit) ;
      filenameDat[strlen(filenameFit) - 3] = 0 ;
      strcat (filenameDat, "dat") ;

      printSM (foutSM, filenameFit, filenameDat) ;
    }
   
  fprintf (foutSM, "hardcopy\n") ;

  if (!feof(fin))
    printf ("\nWarning: Premature end of input file.\n") ;

  printf ("\nCreated '%s'\n", argv[2]) ;
  printf ("\nTo run the SM file, enter: 'sm', then 'inp %s', then 'quit'.\n", argv[2]) ;
  printf ("This will create a postscript file ('%s.ps') of the model fits catalogue.\n\n", argv[2]) ;

  fclose (fin) ;
  fclose (foutSM) ;
  return (0) ;
}
