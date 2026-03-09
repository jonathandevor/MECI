
#define NUM_STAR_TYPES 25

#define NUM_MAIN_SEQUENCE 14  // MS stars are assumed to be at the begining of the array


struct starInfo
{
  char branch ;// {m, g, s} - main sequence ; giants ; supergiants
  char type ;  // {O, B, A, F, G, K, M}
  int num ;    // 0 .. 9
  float B ;    // mag
  float V ;    // mag
  float R ;    // mag
  float I ;    // mag
  float J ;    // mag
  float K ;    // mag
  float r ;    // r_sol
  float rho ;  // g/cm^3
};


void initialize (struct starInfo star[NUM_STAR_TYPES])
{
  //  ------- main squence --------

  // B0
  star[0].branch = 'm' ;
  star[0].type = 'B' ;
  star[0].num = 0 ;
  star[0].V = -4.0 ;
  star[0].B = star[0].V - 0.30 ;
  star[0].R = star[0].V + 0.13 ;
  star[0].I = star[0].R + 0.29 ;
  star[0].K = star[0].V + 0.83 ;
  star[0].J = star[0].K - 0.16 ;
  star[0].r = 7.4 ;
  star[0].rho = 1.4084 * pow (10.0, -1.4) ;

  // B5
  star[1].branch = 'm' ;
  star[1].type = 'B' ;
  star[1].num = 5 ;
  star[1].V = -1.2 ;
  star[1].B = star[1].V - 0.17 ;
  star[1].R = star[1].V + 0.06 ;
  star[1].I = star[1].R + 0.16 ;
  star[1].K = star[1].V + 0.42 ;
  star[1].J = star[1].K - 0.07 ;
  star[1].r = 3.9 ;
  star[1].rho = 1.4084 * pow (10.0, -1.0) ;

  // B8
  star[2].branch = 'm' ;
  star[2].type = 'B' ;
  star[2].num = 8 ;
  star[2].V = -0.25 ;
  star[2].B = star[2].V - 0.11 ;
  star[2].R = star[2].V + 0.02 ;
  star[2].I = star[2].R + 0.10 ;
  star[2].K = star[2].V + 0.24 ;
  star[2].J = star[2].K - 0.03 ;
  star[2].r = 3.0 ;
  star[2].rho = 1.4084 * pow (10.0, -0.85) ;

  // A0
  star[3].branch = 'm' ;
  star[3].type = 'A' ;
  star[3].num = 0 ;
  star[3].V = 0.65 ;
  star[3].B = star[3].V - 0.02 ;
  star[3].R = star[3].V - 0.02 ;
  star[3].I = star[3].R + 0.02 ;
  star[3].K = star[3].V ;
  star[3].J = star[3].K ;
  star[3].r = 2.4 ;
  star[3].rho = 1.4084 * pow (10.0, -0.7) ;

  // A5
  star[4].branch = 'm' ;
  star[4].type = 'A' ;
  star[4].num = 5 ;
  star[4].V = 1.95 ;
  star[4].B = star[4].V + 0.15 ;
  star[4].R = star[4].V - 0.16 ;
  star[4].I = star[4].R - 0.06 ;
  star[4].K = star[4].V - 0.38 ;
  star[4].J = star[4].K + 0.08 ;
  star[4].r = 1.7 ;
  star[4].rho = 1.4084 * pow (10.0, -0.4) ;

  // F0
  star[5].branch = 'm' ;
  star[5].type = 'F' ;
  star[5].num = 0 ;
  star[5].V = 2.7 ;
  star[5].B = star[5].V + 0.30 ;
  star[5].R = star[5].V - 0.30 ;
  star[5].I = star[5].R - 0.17 ;
  star[5].K = star[5].V - 0.70 ;
  star[5].J = star[5].K + 0.16 ;
  star[5].r = 1.5 ;
  star[5].rho = 1.4084 * pow (10.0, -0.3) ;

  // F5
  star[6].branch = 'm' ;
  star[6].type = 'F' ;
  star[6].num = 5 ;
  star[6].V = 3.5 ;
  star[6].B = star[6].V + 0.44 ;
  star[6].R = star[6].V - 0.40 ;
  star[6].I = star[6].R - 0.24 ;
  star[6].K = star[6].V - 1.10 ;
  star[6].J = star[6].K + 0.27 ;
  star[6].r = 1.3 ;
  star[6].rho = 1.4084 * pow (10.0, -0.2) ;

  // G0
  star[7].branch = 'm' ;
  star[7].type = 'G' ;
  star[7].num = 0 ;
  star[7].V = 4.4 ;
  star[7].B = star[7].V + 0.58 ; 
  star[7].R = star[7].V - 0.50 ;
  star[7].I = star[7].R - 0.31 ;
  star[7].K = star[7].V - 1.41 ;
  star[7].J = star[7].K + 0.36 ;
  star[7].r = 1.1 ;
  star[7].rho = 1.4084 * pow (10.0, -0.1) ;

  // G5
  star[8].branch = 'm' ;
  star[8].type = 'G' ;
  star[8].num = 5 ;
  star[8].V = 5.1 ;
  star[8].B = star[8].V + 0.68 ; 
  star[8].R = star[8].V - 0.54 ;
  star[8].I = star[8].R - 0.35 ;
  star[8].K = star[8].V - 1.58 ;  // interpolated
  star[8].J = star[8].K + 0.38 ;  // interpolated
  star[8].r = 0.92 ;
  star[8].rho = 1.4084 * pow (10.0, -0.1) ;

  // K0
  star[9].branch = 'm' ;
  star[9].type = 'K' ;
  star[9].num = 0 ;
  star[9].V = 5.9 ;  
  star[9].B = star[9].V + 0.81 ; 
  star[9].R = star[9].V - 0.64 ;
  star[9].I = star[9].R - 0.42 ;
  star[9].K = star[9].V - 1.96 ;
  star[9].J = star[9].K + 0.53 ;
  star[9].r = 0.85 ;
  star[9].rho = 1.4084 * pow (10.0, 0.1) ;

  // K5  
  star[10].branch = 'm' ;
  star[10].type = 'K' ;
  star[10].num = 5 ;
  star[10].V = 7.35 ;
  star[10].B = star[10].V + 1.15 ; 
  star[10].R = star[10].V - 0.99 ;
  star[10].I = star[10].R - 0.63 ;
  star[10].K = star[10].V - 2.85 ;
  star[10].J = star[10].K + 0.72 ;
  star[10].r = 0.72 ;
  star[10].rho = 1.4084 * pow (10.0, 0.25) ;

  // M0
  star[11].branch = 'm' ;
  star[11].type = 'M' ;
  star[11].num = 0 ;
  star[11].V = 8.8 ;
  star[11].B = star[11].V + 1.40 ; 
  star[11].R = star[11].V - 1.28 ;
  star[11].I = star[11].R - 0.91 ;
  star[11].K = star[11].V - 3.65 ;
  star[11].J = star[11].K + 0.84 ;
  star[11].r = 0.6 ;
  star[11].rho = 1.4084 * pow (10.0, 0.35) ;

  // M2
  star[12].branch = 'm' ;
  star[12].type = 'M' ;
  star[12].num = 2 ;
  star[12].V = 9.9 ;
  star[12].B = star[12].V + 1.49 ;
  star[12].R = star[12].V - 1.50 ;
  star[12].I = star[12].R - 1.19 ;
  star[12].K = star[12].V - 4.11 ;
  star[12].J = star[12].K + 0.86 ;
  star[12].r = 0.5 ;
  star[12].rho = 1.4084 * pow (10.0, 0.8) ;

  // M5
  star[13].branch = 'm' ;
  star[13].type = 'M' ;
  star[13].num = 5 ;
  star[13].V = 12.3 ;
  star[13].B = star[13].V + 1.64 ; 
  star[13].R = star[13].V - 1.80 ;
  star[13].I = star[13].R - 1.67 ;
  star[13].K = star[13].V - 6.17 ;
  star[13].J = star[13].K + 0.95 ;
  star[13].r = 0.27 ;
  star[13].rho = 1.4084 * pow (10.0, 1.0) ;


  //  ------- giant branch --------

  // G5
  star[14].branch = 'g' ;
  star[14].type = 'G' ;
  star[14].num = 5 ;
  star[14].V = 0.9 ;
  star[14].B = star[14].V + 0.86 ;
  star[14].R = star[14].V - 0.69 ;
  star[14].I = star[14].R - 0.48 ;
  star[14].K = star[14].V - 2.1 ;   // interpolated
  star[14].J = star[14].K + 0.57 ;  // interpolated
  star[14].r = 10.0 ;
  star[14].rho = 1.4084 * pow (10.0, -3.0) ;
  
  // K0
  star[15].branch = 'g' ;
  star[15].type = 'K' ;
  star[15].num = 0 ;
  star[15].V = 0.7 ;
  star[15].B = star[15].V + 1.00 ; 
  star[15].R = star[15].V - 0.77 ;
  star[15].I = star[15].R - 0.53 ;
  star[15].K = star[15].V - 2.31 ;
  star[15].J = star[15].K + 0.64 ;
  star[15].r = 15.0 ;
  star[15].rho = 1.4084 * pow (10.0, -3.5) ;
  
  // K5
  star[16].branch = 'g' ;
  star[16].type = 'K' ;
  star[16].num = 5 ;
  star[16].V = -0.2 ;
  star[16].B = star[16].V + 1.50 ; 
  star[16].R = star[16].V - 1.20 ;
  star[16].I = star[16].R - 0.90 ;
  star[16].K = star[16].V - 3.60 ;
  star[16].J = star[16].K + 0.96 ;
  star[16].r = 25.0 ;
  star[16].rho = 1.4084 * pow (10.0, -4.1) ;
  
  // M0
  star[17].branch = 'g' ;
  star[17].type = 'M' ;
  star[17].num = 0 ;
  star[17].V = -0.4 ;
  star[17].B = star[17].V + 1.56 ; 
  star[17].R = star[17].V - 1.23 ;
  star[17].I = star[17].R - 0.94 ;
  star[17].K = star[17].V - 3.85 ;
  star[17].J = star[17].K + 1.02 ;
  star[17].r = 40.0 ;
  star[17].rho = 1.4084 * pow (10.0, -4.7) ;

  //  ------- supergiant branch --------

  // B5
  star[18].branch = 's' ;
  star[18].type = 'B' ;
  star[18].num = 5 ;
  star[18].V = -6.2 ;
  star[18].B = star[18].V - 0.10 ;
  star[18].R = star[18].V - 0.02 ;
  star[18].I = star[18].R + 0.07 ;
  star[18].K = star[18].V + 0.13 ;
  star[18].J = star[18].K + 0.01 ;
  star[18].r = 50.0 ;
  star[18].rho = 1.4084 * pow (10.0, -3.8) ;

 // A0
  star[19].branch = 's' ;
  star[19].type = 'A' ;
  star[19].num = 0 ;
  star[19].V = -6.3 ;
  star[19].B = star[19].V - 0.01 ;
  star[19].R = star[19].V - 0.03 ;
  star[19].I = star[19].R - 0.05 ;
  star[19].K = star[19].V - 0.19 ;
  star[19].J = star[19].K + 0.07 ;
  star[19].r = 60.0 ;
  star[19].rho = 1.4084 * pow (10.0, -4.1) ;

 // A5
  star[20].branch = 's' ;
  star[20].type = 'A' ;
  star[20].num = 5 ;
  star[20].V = -6.6 ;
  star[20].B = star[20].V + 0.09 ;
  star[20].R = star[20].V - 0.12 ;
  star[20].I = star[20].R - 0.13 ;
  star[20].K = star[20].V - 0.48 ;
  star[20].J = star[20].K + 0.15 ;
  star[20].r = 60.0 ;
  star[20].rho = 1.4084 * pow (10.0, -4.2) ;

 // F0
  star[21].branch = 's' ;
  star[21].type = 'F' ;
  star[21].num = 0 ;
  star[21].V = -6.6 ;
  star[21].B = star[21].V + 0.17 ; 
  star[21].R = star[21].V - 0.21 ;
  star[21].I = star[21].R - 0.20 ;
  star[21].K = star[21].V - 0.64 ;
  star[21].J = star[21].K + 0.19 ;
  star[21].r = 80.0 ;
  star[21].rho = 1.4084 * pow (10.0, -4.6) ;

 // F5
  star[22].branch = 's' ;
  star[22].type = 'F' ;
  star[22].num = 5 ;
  star[22].V = -6.6 ;
  star[22].B = star[22].V + 0.32 ;
  star[22].R = star[22].V - 0.35 ;
  star[22].I = star[22].R - 0.23 ;
  star[22].K = star[22].V - 0.93 ;
  star[22].J = star[22].K + 0.28 ;
  star[22].r = 100.0 ;
  star[22].rho = 1.4084 * pow (10.0, -5.0) ;

 // G0
  star[23].branch = 's' ;
  star[23].type = 'G' ;
  star[23].num = 0 ;
  star[23].V = -6.4 ;
  star[23].B = star[23].V + 0.76 ;
  star[23].R = star[23].V - 0.51 ;
  star[23].I = star[23].R - 0.33 ;
  star[23].K = star[23].V - 1.44 ;
  star[23].J = star[23].K + 0.41 ;
  star[23].r = 120.0 ;
  star[23].rho = 1.4084 * pow (10.0, -5.2) ;

 // K0
  star[24].branch = 's' ;
  star[24].type = 'K' ;
  star[24].num = 0 ;
  star[24].V = -6.0 ;
  star[24].B = star[24].V + 1.25 ; 
  star[24].R = star[24].V - 0.76 ;
  star[24].I = star[24].R - 0.48 ;
  star[24].K = star[24].V - 2.15 ;
  star[24].J = star[24].K + 0.58 ;
  star[24].r = 200.0 ;
  star[24].rho = 1.4084 * pow (10.0, -5.8) ;
}
