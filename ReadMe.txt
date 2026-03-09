
(1) Programs:
=============
 meci.c - The MECI source code
 meci_2lc.c - Analyses two light curves at once (same object, different filters)
 debil_lib.c - DEBiL library. Contains MECI's light curve generator
 debil_lib.h - DEBiL library header file



(2)  Tables:
============
 SunIsochrone.txt  - Kim et al.(2002) isochrones (Yonsei-Yale)
 IsochroneBaraffe98_MHsun0.0Lmix1.0.txt - Baraffe et al. (1998) isochrones, Lmix=1.0 (solar-like)
 IsochroneBaraffe98_MHsun0.0Lmix1.9.txt - Baraffe et al. (1998) isochrones, Lmix=1.9 (low mass) 

 LimbDarkTable.txt - Combined limb darkening coefficients from ATLAS (Kurucz 1992) and PHOENIX (Claret 1998, 2000)

 MECI_ReleaseNotes.txt - MECI file formats (isochrone and limb darkening tables, input and output files)
 Readme.txt - This file



(3) MECI utilities:  (All the utility source files are self contained, unless stated otherwise)
===================
 utils/showMECI.c - Creates a SuperMongo script for plotting MECI contours
 utils/filterMECI.c - Filters the results of MECI
 utils/histMECI.c   - Creates histogram plots from the results of MECI
 utils/meci_singleIter.c - Runs MECI on a single pairing iteration

 utils/convert_Debil31.c - Creates a MECI input file from a DEBiL (31-column) database (with no color info)
 utils/convert_OGLE_Debil45.c - Creates a MECI input file from a DEBiL (45-column ; OGLE) database (with color info)
 utils/make_TrES_dblc_UCAC.c - Creates a MECI input file from a DEBiL (31-column) database, TrES 2mass file (with color info), and UCAC database
 utils/make_TrES_dblc_USNO-B.c - Creates a MECI input file from a DEBiL (31-column) database, TrES 2mass file (with color info), and USNO-B database

 utils/meci_sim.c - Creates random-parameter light curves with a MECI input file
 utils/compare_log.c - Compares the MECI output files with the logged "true parameters"



(4) MECI-express:
=================
MECIexp/MECIexp.c - The MECI-express source code (using 31-column DEBiL database)
MECIexp/MECIexp45.c - The MECI-express source code (using 45-column DEBiL database)
MECIexp/starTypes.h - Contains the star characteristics, based on Allen's Astrophysical Quantities (Cox 2000)
