/* ------------------------------------------------------------------------- *
 *                                                                           *
 *                                EZ-Vortex                                  *
 *                    A Code for Simulating Vortex Filament                  *
 *                                                                           *
 *       Copyright (C) 2000 Daniel Margerit                                  *
 *                          margerit@imft.fr                                 *
 *                                                                           *
 *                                Version  1                                 *
 *                                                                           *
 * This code is adapted from  EZ-Scroll a Code for Simulating Scroll Waves   *
 * Copyright (C) 1998  Dwight Barkley                                        *
 * with courtesy of Dwight Barkley (barkley@maths.warwick.ac.uk)             *
 * ------------------------------------------------------------------------- *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 17:49:29 $
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "ezvortex.h"
#include "ezstep3d.h"

/* 
 * Global variables used throughout the EZ-Vortex Code (see ezvortex.h)
 * -------------------------------------------------------------------- */

Real  epsilon, nu_bar_param;
Real  *gamma_param, *delta_0_bar_param, *m_0_param;

Real  *c_v, *c_w, *S, *S_0, *delta_bar, *m_flux, beta;
Real  *int_S_o_S0_tmp, *int_S_o_S0;
Real  *fields, *u, *u_s, *u_ss, *sigma, *u_0, *v, *v_0, *v_m;

Real  *c, *d, *a;

Real  *xm, *xm_2, *trigs; /* Fourier int components */
long int   ifax[13];

Real  ds, dt, error_stop;
int   np, np_fft_o_2, np_fft, nf, n_laguerre_mode, field_size, nsteps, istep, write_filament, n_b, 
      simulating_resolution, rotating_resolution,         
      verbose;
Real  xmin,xmax,ymin,ymax,zmin,zmax;


/* Global variables for this file only
 * ----------------------------------- */

static  Real  length, length_y, ts;
static  int   plot_step, hist_step, hist_x, hist_y, hist_z, 
              ic_type, binary_write, hist_step_tmp;
static  FILE *history_file, *filament_file, *history_file_dat, *fpp; 

/* Private functions 
 * ----------------- */

static void  Initialize      (void);
static void  Allocate_memory (void);
static void  Finish_up       (void);
static void  Write_history   (int wrt_step);
static void  Write_fc        (void);
static void  Write_snapshot  (char *filename);
static void  Write_snapshot_dat (void);
static void  Read_snapshot (void);
static void  Read_ic         (void);
static void  Generate_ic     (void);
static int   factorial       (int);
static Real  EllipticE        (Real);
static Real  EllipticK        (Real);


#define NEXT_LINE(fp) while(getc(fp)!='\n');  /* macro used to skip to end 
						 of input line */

/* ========================================================================= */

void main(void)
{



#if CLOSED
       xmin = -2, xmax = 2.;
       ymin = -2, ymax = 2.;
       zmin =-2, zmax = 2.;
#else
#if 0
       xmin = -0.6250, xmax = 0.6250;
       ymin = -0.6250, ymax = 0.6250;
       zmin = -0.6250, zmax = 0.6250;
#endif
#if 1
       xmin = -10.21/4, xmax =  10.21/4;
       ymin = -10.21/4, ymax =  10.21/4;
       zmin = -10.21/4, zmax =  10.21/4;
#endif
#if 0
       xmin = -7.85/4, xmax =  7.85/4;
       ymin = -7.85/4, ymax =  7.85/4;
       zmin = -7.85/4, zmax =  7.85/4;
#endif
#if 0
       xmin = -0.8976/1.5, xmax = 0.8976/1.5;
       ymin = -0.8976/1.5, ymax = 0.8976/1.5;
       zmin = -0.8976/1.5, zmax = 0.8976/1.5;
#endif

#endif
  
  Initialize();
  
        
#if !COMPUTE 
        hist_step_tmp = hist_step;
	hist_step =1;
        nsteps = (int)((float)nsteps/hist_step_tmp);
        
#endif

  for(istep=0; istep<nsteps; istep++) {
    if( Event_check() )          break;
    if( (istep%plot_step) == 0 ) Draw();


                   
    if( hist_step && (istep%hist_step)==0 ) {
#if COMPUTE
      /*toto*/
                             
                            
                             Write_history(istep);
                            
                          
                           
#endif
#if 1
			      printf("%d\n",istep);
             
#endif
			      
#if MOVIE
    if(1 && (istep%25)==0){ /* time step to do snapshot  */
			      Save_image();
                              printf("isauve=%d\n",istep);
    }
#endif
#if COMPUTE
			     
			      Write_snapshot_dat();
                            
#endif
    }
#if COMPUTE
     
      Step(fpp); 
   
      
#else
      if (istep!=0){
	if( hist_step && (istep%hist_step)==0 ) {
                 Read_snapshot ();
#if MOVE
#if 0
      zmax -=ts*hist_step_tmp/(2*M_PI);
      zmin -=ts*hist_step_tmp/(2*M_PI);
      Mover_continuous(0.,0.,-ts*hist_step_tmp/(2*M_PI));
#endif
#if 1
      zmax -=ts*hist_step_tmp*gamma_param[0]/(2*M_PI)/0.5;
      zmin -=ts*hist_step_tmp*gamma_param[0]/(2*M_PI)/0.5;
      Mover_continuous(0.,0.,-ts*hist_step_tmp*gamma_param[0]/(2*M_PI)/0.5);
#endif
#if 0
      zmax -=ts*hist_step_tmp*0.19545174621606;
      zmin -=ts*hist_step_tmp*0.19545174621606;
      Mover_continuous(0.,0.,-ts*hist_step_tmp*0.19545174621606);
#endif
#endif
	}}
      
#endif

#if COMPUTE
#if MOVE
#if 0
      zmax -=ts/(2*M_PI);
      zmin -=ts/(2*M_PI);
      Mover_continuous(0.,0.,-ts/(2*M_PI));
#endif
#if 1
      zmax -=ts*gamma_param[0]/(2*M_PI)/0.5;
      zmin -=ts*gamma_param[0]/(2*M_PI)/0.5;
      Mover_continuous(0.,0.,-ts*gamma_param[0]/(2*M_PI)/0.5);
#endif
#if 0
      zmax -=ts*0.19545174621606;
      zmin -=ts*0.19545174621606;
      Mover_continuous(0.,0.,-ts*0.19545174621606);
#endif
#endif
#endif
  }
  
  Finish_up();
}
/* ========================================================================= */

void Initialize (void)
{
  /* Read task file, open output files, and initialize graphics and time
   * stepping.  */

  
  double p_in;
  FILE *fp;
  int initial_mode, j, i;
  time_t tt1;

 Real S_tmp;
 Real s_step =  2*M_PI/(NP-1);
 Real point_m[3], point_c[3], point_p[3], vec_tmp[3];
 Real one_o_2ds;
 Real tmp;


 Real S0_o_S, one_nu_tmp, delta_0_bar_param_4, coeff_tmp, S0_o_S_4;
 int in, im;


 /* ----------------------------- 
   * Write informations 
   * ----------------------------- */
 

  printf("\n\n\n\n\n\n\n\n\n\n");
  printf("              EZ_Vortex Simulation\n\n\n\n");
#if CLOSED
  printf("The filaments are closed.\n");
#else
  printf("The filaments are open.\n");
#endif

#if LOCAL_INDUCTION 
    printf("The local induction equation of motion is used.\n");
#endif

#if  CALL_AND_TING
    printf("The Callegari and Ting equation of motion is used.\n");
#endif

#if  DE_SINGU
    printf("The De-singlarization equation of motion is used.\n");
    printf("It requires much more points than in the Callegari and Ting or local induction.\n");
#endif

#if  M1_KNIO_KLEIN
    printf("The Method M1 of Klein/Knio is used.\n");
#endif



#if STRUCTURE
  printf("The inner structure is taken into account\n");
#else
  printf("The inner structure is NOT taken into account.\n");
#endif

#if  SPECTRAL
    printf("The derivatives are computed spectrally.\n");
#else
  printf("The derivatives are computed with finite differences .\n");
#endif

#if  EXPLICIT
    printf("The explicit time-steping is used.\n");
#if GAUSS_SEIDEL
    printf("The Gauss-Seidel method is used.\n");
#else
    printf("The Jacobi method is used.\n");
#endif
#endif

#if  NEWTON
    printf("The implicit time-steping with Newton solver is used.\n");
#endif

#if  ADAMS_BASHFORTH
    printf("The Adams Bashforth time-steping is used.\n");
#if GAUSS_SEIDEL
    printf("The Gauss-Seidel method is used.\n");
#else
    printf("The Jacobi method is used.\n");
#endif
#endif




 /* ----------------------------- 
   * Control of the spatial macros 
   * ----------------------------- */

  if ((((LOCAL_INDUCTION||CALL_AND_TING)||DE_SINGU)||M1_KNIO_KLEIN)==0){
    printf("\n\n\n\n\n");
    printf("Either LOCAL_INDUCTION or CALL_AND_TING or DE_SINGU or M1_KNIO_KLEIN has to be set to 1 in ezvortex.h\n\n");
    exit(0);
  }
  
   if ((((LOCAL_INDUCTION&&CALL_AND_TING)==1)||((LOCAL_INDUCTION&&DE_SINGU)==1)
	                                     ||((CALL_AND_TING&&DE_SINGU)==1))==1){
    printf("\n\n\n\n\n");
    printf("Only one of  LOCAL_INDUCTION or CALL_AND_TING or DE_SINGU or M1_KNIO_KLEIN has to be set to 1 in ezvortex.h\n\n");
   exit(0);
   }

   if ((((M1_KNIO_KLEIN&&CALL_AND_TING)==1)||((M1_KNIO_KLEIN&&DE_SINGU)==1)
	                                     ||((M1_KNIO_KLEIN&&LOCAL_INDUCTION)==1))==1){
    printf("\n\n\n\n\n");
    printf("Only one of  LOCAL_INDUCTION or CALL_AND_TING or DE_SINGU or M1_KNIO_KLEIN has to be set to 1 in ezvortex.h\n\n");
   exit(0);
   }

 /* ----------------------------- 
   * Control of the stepping macros 
   * ----------------------------- */


  if (((EXPLICIT ||NEWTON)||ADAMS_BASHFORTH)==0){
    printf("\n\n\n\n\n");
    printf("Either EXPLICIT or NEWTON or ADAMS_BASHFORTH  has to be set to 1 in ezvortex.h\n\n");
    exit(0);
  }
  
   if ((((EXPLICIT&&NEWTON)==1)||((EXPLICIT&&ADAMS_BASHFORTH)==1)
	                                     ||((NEWTON&&ADAMS_BASHFORTH)==1))==1){
    printf("\n\n\n\n\n");
    printf("Only one of  EXPLICIT or NEWTON or DE_SINGU or ADAMS_BASHFORTH has to be set to 1 in ezvortex.h\n\n");
   exit(0);
   }





#if COMPUTE 
  /* ----------------------------- 
   * Read parameters from task.dat 
   * ----------------------------- */

  if((fp=fopen("task.dat","r"))==NULL) {
    if(verbose) fprintf(stderr,"Cannot open task file: task.dat \n");
    exit(1);
  }
  else {
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); nu_bar_param=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); epsilon=p_in;
                             NEXT_LINE(fp);
    fscanf(fp,"%d", &np);       
    fscanf(fp,",%d",&nf);    NEXT_LINE(fp);
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); ts=p_in;
    fscanf(fp,"%lg",&p_in);  NEXT_LINE(fp); error_stop=p_in;
                             NEXT_LINE(fp);
    fscanf(fp,"%d",&nsteps); NEXT_LINE(fp);
    fscanf(fp,"%d",&plot_step);       NEXT_LINE(fp);
    fscanf(fp,"%d",&write_filament);  NEXT_LINE(fp);
    fscanf(fp,"%d",&hist_step);       NEXT_LINE(fp);
    fscanf(fp,"%d",&n_b);             NEXT_LINE(fp);
    fscanf(fp,"%d,%d,%d",&hist_x,&hist_y,&hist_z); NEXT_LINE(fp);
                                            NEXT_LINE(fp);
    fscanf(fp,"%d",&initial_mode);         NEXT_LINE(fp);
    fscanf(fp,"%d",&ic_type);               NEXT_LINE(fp);
    fscanf(fp,"%d", &simulating_resolution); 
    fscanf(fp,",%d",&rotating_resolution);  NEXT_LINE(fp);
    fscanf(fp,"%d",&binary_write);          NEXT_LINE(fp);
                                            NEXT_LINE(fp);
    fscanf(fp,"%d",&verbose);

    fclose(fp);
  }

#if CONV_ANALYSE
    nsteps = 1;
#endif

#if 1
  if (NP%2==0) {
       printf("The number of points np should be an odd number !!!!\n\n");
       exit(0);
    }
#endif

#if SPECTRAL
 if (NP!=257) {
       printf("\n The number of points np should be 257 !!!!\n\n");
       exit(0);
    }
#endif

#else 
 /* ----------------------------- 
   * Read parameters from history.dat 
   * ----------------------------- */
  if((history_file_dat=fopen("history.dat","r"))==NULL) { 
    printf("No file history.dat !!!!");
    exit(1);
  }
    fscanf(history_file_dat,"%lg",&p_in);  
    NEXT_LINE(history_file_dat); nu_bar_param=p_in;
    fscanf(history_file_dat,"%lg",&p_in);  
    NEXT_LINE(history_file_dat); epsilon=p_in;
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d", &np);       
    fscanf(history_file_dat,",%d",&nf);    
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%lg",&p_in);
    NEXT_LINE(history_file_dat); ts=p_in;
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d",&nsteps); 
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d",&plot_step);      
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d",&hist_step);       
    NEXT_LINE(history_file_dat);
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d",&initial_mode);         
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d", &simulating_resolution); 
    fscanf(history_file_dat,",%d",&rotating_resolution);  
    NEXT_LINE(history_file_dat);
    NEXT_LINE(history_file_dat);
    fscanf(history_file_dat,"%d",&verbose);
  
#endif


#if SPECTRAL
   np_fft_o_2 = (NP-1)/2;
   np_fft = (NP-1);
#endif



    gamma_param = VECTOR(NF);
    m_0_param = VECTOR(NF);
    delta_0_bar_param = VECTOR(NF);
    c_v = VECTOR(NF);
    c_w = VECTOR(NF);
    delta_bar = VECTOR(NF);
    S_0  = VECTOR(NF);
    S = VECTOR(NF);
    m_flux = VECTOR(NF);
    int_S_o_S0_tmp = VECTOR(NF);
    int_S_o_S0  = VECTOR(NF);



  /* ----------------------------------------------------------------- 
   * Process input parameters and write various things to screen 
   *
   * Note:
   * ----------------------------------------------------------------- */

  Allocate_memory();
  Read_ic();


  /* ----------- Non-similar  Core IC   -------------------------------- */
#if NON_SIMIL_PART     

for(j=0;j<=NF-1;j++) {
     DD(0,j) = gamma_param[j]/2/M_PI/delta_0_bar_param[j]/delta_0_bar_param[j];
     CC(0,j) = m_flux[j]/2/M_PI/delta_0_bar_param[j]/delta_0_bar_param[j];

#if 0
     DD(1,j) = 0;
     CC(1,j) = 0;
  
     DD(2,j) = 0.4;
     CC(2,j) = 0.4;
     DD(3,j) = 0.2;
     CC(3,j) = 0.2;
     DD(4,j) = 0.1;
     CC(4,j) = 0.1;
     DD(5,j) = 0.05;
     CC(5,j) = 0.05;
#endif

#if 0 /*  top hat */
     DD(1,j) = 0.1000*DD(0,j);
     CC(1,j) = 0;
     DD(2,j) = -0.1571*DD(0,j);
     CC(2,j) = 0.;
     DD(3,j) = -0.1464*DD(0,j);
     CC(3,j) = 0.;
     DD(4,j) = -0.0588*DD(0,j);
     CC(4,j) = 0.;
     DD(5,j) = 0.0211*DD(0,j);
     CC(5,j) = 0.0;
     DD(6,j) = 0.0658*DD(0,j);
     CC(6,j) = 0;
     DD(7,j) = 0.0755*DD(0,j);
     CC(7,j) = 0.;
     DD(8,j) = 0.0611*DD(0,j);
     CC(8,j) = 0.;
     DD(9,j) = 0.0352*DD(0,j);
     CC(9,j) = 0.;
     DD(10,j) = 0.0080*DD(0,j);
     CC(10,j) = 0.0;
     DD(11,j) = -0.0138*DD(0,j);
     CC(11,j) = 0;
     DD(12,j) = -0.0272*DD(0,j);
     CC(12,j) = 0.;
     DD(13,j) = -0.0318*DD(0,j);
     CC(13,j) = 0.;
     DD(14,j) = -0.0290*DD(0,j);
     CC(14,j) = 0.;
     DD(15,j) = -0.0211*DD(0,j);
     CC(15,j) = 0.0;
     DD(16,j) = -0.0106*DD(0,j);
     CC(16,j) = 0;
     DD(17,j) = 0.0002*DD(0,j);
     CC(17,j) = 0.;
     DD(18,j) = 0.0096*DD(0,j);
     CC(18,j) = 0.;
     DD(19,j) = 0.0165*DD(0,j);
     CC(19,j) = 0.;
     DD(20,j) = 0.0204*DD(0,j);
     CC(20,j) = 0.0;
#endif

#if 0 /*  Rankine */
     DD(1,j) = 0.5000000000*DD(0,j);
     CC(1,j) = 0;
     DD(2,j) = 0.1666666667*DD(0,j);
     CC(2,j) = 0.;
     DD(3,j) = -0.04166666667*DD(0,j);
     CC(3,j) = 0.;
     DD(4,j) = -0.1583333333*DD(0,j);
     CC(4,j) = 0.;
     DD(5,j) = -0.2097222222*DD(0,j);
     CC(5,j) = 0.0;
     DD(6,j) = -0.2164682540*DD(0,j);
     CC(6,j) = 0;
     DD(7,j) = -0.1944692460*DD(0,j);
     CC(7,j) = 0.;
     DD(8,j) = -0.1557512125*DD(0,j);
     CC(8,j) = 0.;
     DD(9,j) = -0.1092016645*DD(0,j);
     CC(9,j) = 0.;
     DD(10,j) = -0.06118824655*DD(0,j);
     CC(10,j) = 0.0;
     DD(11,j) = -0.01607804442*DD(0,j);
     CC(11,j) = 0;
     DD(12,j) = 0.02332889927*DD(0,j);
     CC(12,j) = 0.;
     DD(13,j) = 0.05543992963*DD(0,j);
     CC(13,j) = 0.;
     DD(14,j) = 0.07957349396*DD(0,j);
     CC(14,j) = 0.;
     DD(15,j) = 0.09571701938*DD(0,j);
     CC(15,j) = 0.0;
     DD(16,j) = 0.1043308936*DD(0,j);
     CC(16,j) = 0;
     DD(17,j) = 0.1061915099*DD(0,j);
     CC(17,j) = 0.;
     DD(18,j) = 0.1022672451*DD(0,j);
     CC(18,j) = 0.;
     DD(19,j) = 0.09362204447*DD(0,j);
     CC(19,j) = 0.;
     DD(20,j) = 0.08134200369*DD(0,j);
     CC(20,j) = 0.0;
#endif


}
               

 for(j=0;j<=N_LAGUERRE_MODE;j++) {
      for(i=0;i<N_LAGUERRE_MODE;i++) {
      AA(i,j) = factorial(i+j)/(factorial(i)*factorial(j)*exp((i+j+1)*log(2)));
      }
  }
#endif

 /* ----------------------------------------------------------------- */

#if SPECTRAL
/* ----------- Spectral int coefficients   -------------------------------- */

      for(i=1;i<=np_fft_o_2;i++) {
      j = 2 * i - 1;
      tmp = (Real)(i-1);
      xm[j-1] = tmp;
      xm[j] = tmp;
      }

      for(i=1;i<=NP-1;i++) {
      xm_2[i-1] = xm[i-1] * xm[i-1];
      }
     
     
      set99_(trigs, ifax, &np_fft);
   
      

/* ----------------------------------------------------------------- */
#endif


  Write_snapshot("ic_dat.m");

 
  dt         = ts;
  ds         = 2.*M_PI/(NP-1.);
  one_o_2ds = 1. / (ds * 2);
 

  for(j=0;j<=NF-1;j++) {

      /* Initialyse  Ux_s and SIGMA(i,j) */

     for(i=1;i<=NP;i++) {
                  point_m[X] = Ux(i-1,j); 
	     point_m[Y] = Uy(i-1,j); 
	     point_m[Z] = Uz(i-1,j);

             point_c[X] = Ux(i,j); 
	     point_c[Y] = Uy(i,j); 
	     point_c[Z] = Uz(i,j);

             point_p[X] = Ux(i+1,j); 
	     point_p[Y] = Uy(i+1,j); 
	     point_p[Z] = Uz(i+1,j);

	     /*  X_s derivative and find sigma */
             SUB(vec_tmp,point_p,point_m); 
             Ux_s(i,j) =  vec_tmp[X] * one_o_2ds;      
             Uy_s(i,j) =  vec_tmp[Y] * one_o_2ds;  
             Uz_s(i,j) =  vec_tmp[Z] * one_o_2ds;        
	     SIGMA(i,j) = NORMV(vec_tmp) * one_o_2ds;

     }



      /* Find the length S0(t) */
         S_tmp =0.;
         for(i=1;i<=NP-1;i++) {
         S_tmp +=SIGMA(i,j);
	 }
         S_0[j] = S_tmp * ds;
        
      /* Initialyse other parameters */

  c_v[j]         = C_V(delta_0_bar_param[j]);
  c_w[j]         = C_W(delta_0_bar_param[j],m_0_param[j],gamma_param[j]); 
  delta_bar[j]   = delta_0_bar_param[j]; 
  S0_o_S = 1.;

#if NON_SIMIL_PART
        FIND_NON_SIMIL_PART
#endif

#if 0  
 printf("%g\n%g\n",c_v[j], c_w[j]);
 exit(0);
#endif

  m_flux[j] = m_0_param[j];
  int_S_o_S0[j] = 0;
  S[j] = S_0[j];
 
  }

 

  if(verbose) {
    printf("\n\nModel Parameters: \n");
    printf("epsilon   = %g\n", epsilon);
    printf("nu_bar     = %g\n", nu_bar_param);
   for(j=0;j<=NF-1;j++) {
    printf("gamma[%d]     = %g\n", j , gamma_param[j]);
    printf("delta_0_bar[%d]   = %g\n",j , delta_0_bar_param[j]);
    printf("m_0[%d]       = %g\n", j, m_0_param[j]);
   }
    printf("\nNumerical Parameters: \n");
    printf("NP = %-6d NF = %-6d  N_b= %-6d ts  = %-10g\n", 
	   NP, NF, n_b, ts);
    printf("dt   = %-10g ds= %-10g\n", 
	   dt, ds);
#if BIOT_ON       
    printf("; Biot Integral --\n");
#else
    printf(" --\n");
#endif


#if EXPLICIT
    printf("Explicit ");
#endif
#if NEWTON 
    printf("Implicit with Newton solver");
#endif
#if  ADAMS_BASHFORTH
    printf("Adams-Bashforth stepping");
#endif

    printf("\nNumber of time steps = %d\n", nsteps);
    printf("  time steps per plot (and/or filament computation) = %d\n", 
	   plot_step);
    if(write_filament) printf("  writing filament data\n");
    else printf("  not writing filament data\n");
    if(hist_step) printf("  writing history data\n");
    else printf("  not writing history data\n");
    printf("\n\n");
  }
#if 0
  /* ------------------ 
   * Perform some tests 
   * ------------------ */

  if( V_DIFF_ON && Dv==0.) {
    fprintf(stderr,"***** V_DIFF_ON is 1 and Dv == 0. ******\n");
    exit(1);
  }
  if(!V_DIFF_ON && Dv!=0.) {
    fprintf(stderr,"***** V_DIFF_ON is 0 and Dv != 0. ******\n");
    exit(1);
  }
  if(hist_step && (
     hist_x<1 || hist_x>NX || 
     hist_y<1 || hist_y>NY || 
     hist_z<1 || hist_z>NZ ) ) {
    fprintf(stderr,"***** history point out of range ******\n");
    exit(1);
  }
  if(ts > 1.0 ) {
    fprintf(stderr,"***** ts > 1 (the diffusion stability limit) ******\n");
    exit(1);
  }
#endif

  /* ------------ 
   * Final things
   * ------------ */

 
  Step_ini();

#if COMPUTE

  /* toto */
  
  fpp = fopen("look", "w");
   
  if(hist_step){

  history_file = fopen("history_dat.m", "w");
 

  fprintf(history_file,"epsilon=%g;\n",epsilon);
  fprintf(history_file,"nu_bar=%g;\n",nu_bar_param);
  for(j=0;j<=NF-1;j++) {
    fprintf(history_file,"gamma(%d)=%g;\n", j+1, gamma_param[j]);
    fprintf(history_file,"m0(%d)=%g;\n", j+1, m_0_param[j]);
  }

  fprintf(history_file,"error_stop=%g;\n",error_stop);


  fprintf(history_file,"np=%d;\n",NP);
  fprintf(history_file,"nf=%d;\n",NF);
  fprintf(history_file,"nb=%d;\n",n_b);
  fprintf(history_file,"ts=%g;\n",ts);
  fprintf(history_file,"nsteps=%d;\n",nsteps);
  fprintf(history_file,"hist_step=%d;\n",hist_step);}
#endif


#if COMPUTE 
 if(hist_step){


  history_file_dat = fopen("history.dat", "w");

    fprintf(history_file_dat,"%g \t\t parameter nu_bar\n",nu_bar_param); 
    fprintf(history_file_dat,"%g \t\t parameter epsilon\n\n",epsilon);
    fprintf(history_file_dat,"%d,%d \t\t np,nf (ODD Numb of pts per filament(<= NP if ic.dat), Num of filaments : this number must be coherent with the initial condition type set below) \n",np,nf);
    fprintf(history_file_dat,"%g \t ts (Timestep)\n\n",ts);
    fprintf(history_file_dat,"%d \t\t number of time steps to take\n",nsteps);
    fprintf(history_file_dat,"%d \t\t time steps per plot\n",plot_step);
    fprintf(history_file_dat,"%d \t\t time steps per write to history file\n\n",hist_step);
    fprintf(history_file_dat,"%d \t\t number of periodic boxes (open filament) \n",n_b);
    fprintf(history_file_dat,"%d \t\t initial display mode of filament : 0 for curve, 1 for worm\n",initial_mode);
    fprintf(history_file_dat,"%d,%d \t\t simulation resolution, rotation resolution\n\n",simulating_resolution,rotating_resolution);
    fprintf(history_file_dat,"%d \t\t verbose: Level of diagnostics\n",verbose);
  

   for(j=0;j<=NF-1;j++) {
   fprintf(history_file_dat,"gamma_param[%d],delta_bar[%d],m_flux[%d] =%g, %g, %g\n",j,j,j,gamma_param[j],delta_bar[j],m_flux[j]);
   }
  

  time(&tt1); 
  fprintf(history_file_dat,"File written: %s",ctime(&tt1));
  fprintf(history_file_dat,"Comments: 3D \n");
  fprintf(history_file_dat,"\n");
  fprintf(history_file_dat,"\n");
  fprintf(history_file_dat,"\n");
  fprintf(history_file_dat,"\n");
  fprintf(history_file_dat,"\n");}
#endif

  if(write_filament) filament_file = fopen("filament.dat", "w"); 
  Draw_ini(initial_mode);
}
/* ========================================================================= */

void Allocate_memory (void)
{
  /* -----------------------------------------------------------------------
   * There are NX, NY and NZ actual grid points in each direction.
   * field_size = (NX+2)*(NY+2)*(NZ+2) because of the way boundaries are
   * treated.  See ezscroll.h for definition of NX, NY, NZ.
   *                                                                         
   * I try to allocate one big block of memory for everything (u, v,
   * sigma_u, and if needed sigma_v) and then set appropriate pointers.  If
   * this fails, try for blocks of size 2*field_size, which is the smallest
   * block I can accept.  The handling of memory can be changed with
   * appropriate changes to the macros INDEX etc in ezscroll.h.
   * ----------------------------------------------------------------------- */

  field_size = (NP+2)*3*NF;

#if EXPLICITE 
#if GAUSS_SEIDEL
  fields = (Real *) calloc(4*field_size + (NP+2)*NF, sizeof(Real));
#else
   fields = (Real *) calloc(5*field_size + (NP+2)*NF, sizeof(Real));
#endif
#endif
#if NEWTON
  fields = (Real *) calloc(5*field_size + (NP+2)*NF, sizeof(Real));
#endif
#if ADAMS_BASHFORTH
#if GAUSS_SEIDEL
  fields = (Real *) calloc(7*field_size + (NP+2)*NF, sizeof(Real));
#else
  fields = (Real *) calloc(8*field_size + (NP+2)*NF, sizeof(Real));
#endif
#endif

 /* -In some cases sigma or u_s or u_ss is not needed and this place can be shorter */



  if (fields != NULL) {
    /* --------------------------------------------------- 
     * We have one big block of memory.  Set some pointers 
     * --------------------------------------------------- */
    u         = fields;
    u_s       = u      +   field_size;
    u_ss      = u_s    +   field_size;
    sigma     = u_ss   +   field_size;
#if NEWTON
    u_0       = sigma  +   (NP+2)*NF;
#endif

#if EXPLICIT
#if !GAUSS_SEIDEL
    u_0       = sigma  +   (NP+2)*NF;
#endif
#endif



#if ADAMS_BASHFORTH
    v         = sigma  +   (NP+2)*NF;
    v_0       = v      +   field_size;
    v_m       = v_0    +   field_size;
#if !GAUSS_SEIDEL
    u_0       = v_m    +   field_size;
#endif
#endif

    if(verbose>2) printf("Memory allocated in one block.\n");
  }
  else {
    /* ------------------------------------------------------------------- 
     * Did not get one big block. Try to find blocks of size field_size.
     * ------------------------------------------------------------------- */

#define CALL_CALLOC     (Real *) calloc(field_size, sizeof(Real))
#define CALL_CALLOC_SCAL (Real *) calloc((NP+2)*NF, sizeof(Real))

    if(    ((fields  = CALL_CALLOC ) != NULL) 
	&& ((u_s = CALL_CALLOC ) != NULL)
	&& ((u_ss = CALL_CALLOC ) != NULL)
 	&& ((sigma = CALL_CALLOC_SCAL ) != NULL)
#if  NEWTON
	&& ((u_0 = CALL_CALLOC ) != NULL)
#endif
#if EXPLICIT
#if !GAUSS_SEIDEL
        && ((u_0 = CALL_CALLOC ) != NULL)
#endif
#endif

#if ADAMS_BASHFORTH
        && ((v = CALL_CALLOC ) != NULL)
 	&& ((v_0 = CALL_CALLOC ) != NULL)
       	&& ((v_m = CALL_CALLOC ) != NULL)
#if !GAUSS_SEIDEL
        && ((u_0 = CALL_CALLOC ) != NULL)
#endif
#endif

       ) {
      u         = fields;
      /* Found memory
       * ------------ */
      if(verbose>2) printf("Using multiple memory blocks.\n");
    }
    else {
      /* Could not find enough memory.
       * ----------------------------- */
      fprintf(stderr, "\n ****Error: Not enough memory available\n");
      exit(-1);
    }
  }

#if NON_SIMIL_PART  
  N_LAGUERRE_MODE = 20;
  d = (Real *) calloc((N_LAGUERRE_MODE+1)*NF, sizeof(Real));
  c = (Real *) calloc((N_LAGUERRE_MODE+1)*NF, sizeof(Real));
  a = (Real *) calloc((N_LAGUERRE_MODE+1)*(N_LAGUERRE_MODE+1), sizeof(Real));
#endif

#if SPECTRAL
  xm = (Real *) calloc(NP+1,sizeof(Real));
  xm_2 = (Real *) calloc(NP+1, sizeof(Real));
  trigs = (Real *) calloc(3*((NP-1)/2)+1, sizeof(Real));
#endif





}
/* ========================================================================= */

void Finish_up (void)
{ 
  QuitX(); 
  Write_fc(); 
  Write_snapshot("fc_dat.m");
  if(hist_step) fclose(history_file);
  if(hist_step) fclose(history_file_dat);

  free(u);
  free(u_s);
  free(u_ss);

#if SPECTRAL
  free(xm);
  free(xm_2);
  free(trigs);
#endif


}
/* ========================================================================= */

void Write_filament_data (float x0, float y0, float z0, 
			  float x1, float y1, float z1) 
{
}
/* ========================================================================= */

void Write_history (int wrt_step)
{
    /* 
      *  Write  conditions to a file history_dat.m
      */
  int i,j;
  int wrt_inc = wrt_step/hist_step * NF;

  
  for(i=1;i<=NF;i++) {
    fprintf (history_file, "X(%d,:)=[\n",i+wrt_inc);
    for(j=1;j<=NP;j++){
      fprintf (history_file, "%.17e\n", Ux(j,i-1));
    }
   fprintf (history_file, "]';\n");
  }

 for(i=1;i<=NF;i++) {
    fprintf (history_file, "Y(%d,:)=[\n",i+wrt_inc);
    for(j=1;j<=NP;j++){
      fprintf (history_file, "%.17e\n", Uy(j,i-1));
    }
   fprintf (history_file, "]';\n");
  }

 for(i=1;i<=NF;i++) {
    fprintf (history_file, "Z(%d,:)=[\n",i+wrt_inc);
    for(j=1;j<=NP;j++){
      fprintf (history_file, "%.17e\n", Uz(j,i-1));
    }
   fprintf (history_file, "]';\n");
  }

 


}
/* ========================================================================= */


void Write_fc (void)
{
/* Write final condition file */

  FILE *fp;
  time_t tt1;
  int i,j,k;

  time(&tt1); 
  fp = fopen("fc.dat", "w");

  fprintf(fp,"Model Parameters: nu_bar, epsilon = ");
  fprintf(fp,"%g, %g\n", nu_bar_param, epsilon);
   for(j=0;j<=NF-1;j++) {
   fprintf(fp,"gamma_param[%d],delta_bar[%d],m_flux[%d] =%g, %g, %g\n",j,j,j,gamma_param[j],delta_bar[j],m_flux[j]);
   }
  fprintf(fp,"Numerical Parameters: NP, NF,N_b, ts = ");
  fprintf(fp,"%d, %d, %d, %g \n", NP, NF, n_b, ts);
  fprintf(fp,"File written: %s",ctime(&tt1));
  fprintf(fp,"Comments: 3D \n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");

  if(binary_write) {
    /* Write binary data */
    fprintf(fp,"Binary values of ux, uy and uz follow\n");

    for(j=1;j<=NF;j++) {
	fwrite(&Ux(0,j-1), sizeof(Real), (NP+2)*3, fp);
    }    
 
  }
  else {
    /* Write ascii data */
    fprintf(fp,"Ascii values of ux, uy and uz follow\n");
    
      for(j=1;j<=NF;j++) {
	for(i=0;i<=NP+1;i++) {
	  fprintf (fp, "%g %g %g\n", (float)Ux(i,j-1), (float)Uy(i,j-1), (float)Uz(i,j-1));
	}
      }
  }  
  fclose(fp); 
}
/*---------------------------------------------------------------------------*/

void Write_snapshot_dat (void)

     /* 
      *  Write  conditions to a file history.dat 
      */
{
  time_t tt1;
  int i,j;


  time(&tt1); 
  
    if(binary_write) {
    /* Write binary data */
    fprintf(history_file_dat,"Binary values of ux, uy and uz follow\n");

    for(j=1;j<=NF;j++) {
	fwrite(&Ux(0,j-1), sizeof(Real), (NP+2)*3, history_file_dat);
    }    
 
  }
  else {
    /* Write ascii data */
    fprintf(history_file_dat,"Ascii values of ux, uy and uz follow\n");
    
      for(j=1;j<=NF;j++) {
	for(i=0;i<=NP+1;i++) {
	  fprintf (history_file_dat, "%g %g %g\n", (float)Ux(i,j-1), (float)Uy(i,j-1), (float)Uz(i,j-1));
	}
      }
  }  


}



/*---------------------------------------------------------------------------*/
/* ========================================================================= */

void Read_snapshot (void)
{ 

  /* Reads file (history.dat) 
   *  */

  double ux_in, uy_in, uz_in;
  double gamma_in, delta_in, m0_in;

  Real   *u_tmp, rp;
  int    np_ic, nf_ic, npts_ic, npts_ic3, index, i_tmp, j_tmp, 
         i, j;
  char   f_type, dummy;
 

 
  /* Skip to next line and read first character to determine type 
     B(inary) or A(scii) */
 /*  NEXT_LINE(history_file_dat); */
 
  f_type = getc(history_file_dat); NEXT_LINE(history_file_dat); 
  
  if ( (f_type !='B') && (f_type !='A') ) {
    if(verbose) 
      printf("\n history.dat exists but of unrecognized type Binary or Ascii \n"); 
    exit(1);
  }
 
  np_ic = np;
  nf_ic = nf;

#if 0
  if(verbose) printf("\nReading history.dat with np, nf = %d, %d... \n\n", 
		     np_ic, nf_ic);
#endif

  npts_ic  = (np_ic+2) * nf_ic;
  npts_ic3  = (np_ic+2) * 3 *  nf_ic;  
 

  /* Allocate temporary memory and read from file */


  u_tmp =(Real *) malloc((unsigned)(npts_ic3)*sizeof(Real));
 

  if(f_type =='B') {    
    /* Binary data file */
    fread(u_tmp, sizeof(Real), npts_ic3, history_file_dat);
  }
  else {          
    /* Ascii data file */
    for(index=0;index<npts_ic;index++) {
      fscanf (history_file_dat, "%lg %lg %lg\n", &ux_in, &uy_in, &uz_in);
     
      /*  printf("%lg %lg %lg\n", ux_in, uy_in, uz_in);*/
      
  
      u_tmp[index*3] = ux_in;
      u_tmp[index*3+1] = uy_in;
      u_tmp[index*3+2] = uz_in;
    }
  }

  /* Copy into u and v */
  rp = (np_ic+2.-1.0)/(NP+2.-1.0);
 
      

  for(j=1;j<=NF;j++) {
        j_tmp = 3*(NP+2)*(j-1);
      for(i=1;i<=NP+2;i++) {
	i_tmp = (int)(rp*(i-1));
	Ux(i-1,j-1) = u_tmp[i_tmp * 3 + j_tmp];
        Uy(i-1,j-1) = u_tmp[i_tmp * 3 + j_tmp + 1];
	Uz(i-1,j-1) = u_tmp[i_tmp * 3 + j_tmp + 2];
     /*  printf("%lg %lg %lg\n", Ux(i,j), Uy(i,j), Uz(i,j)); */
      }
  }



  free(u_tmp);



}
/* ========================================================================= */



void Write_snapshot (char *filename)

     /* 
      *  Write final conditions to a matlab file fc_dat.m. 
      */
{
  FILE *fp;
  time_t tt1;
  int i,j;

  time(&tt1); 
  fp = fopen(filename, "w");

  fprintf(fp,"nu_bar=%g;\n",nu_bar_param);
   fprintf(fp,"epsilon=%g;\n",epsilon);
  for(j=0;j<=NF-1;j++) {
  fprintf(fp,"gamma(%d)=%g;\n",j+1,gamma_param[j]);
  fprintf(fp,"m0(%d)=%g;\n",j+1,m_0_param[j]);
  }

 
  fprintf(fp,"np=%d;\n",NP);
  fprintf(fp,"nf=%d;\n",NF);
  fprintf(fp,"nb=%d;\n",n_b);
  fprintf(fp,"ts=%g;\n",ts);
  fprintf(fp,"error_stop=%g;\n",error_stop);
  fprintf(fp,"nsteps=%d;\n",nsteps);
  fprintf(fp,"hist_step=%d;\n",hist_step);

  for(i=1;i<=NF;i++) {
    fprintf (fp, "Ux(%d,:)=[\n",i);
    for(j=1;j<=NP;j++){
      fprintf (fp, "%g\n", Ux(j,i-1));
    }
   fprintf (fp, "]';\n");
  }

 for(i=1;i<=NF;i++) {
    fprintf (fp, "Uy(%d,:)=[\n",i);
    for(j=1;j<=NP;j++){
      fprintf (fp, "%g\n", Uy(j,i-1));
    }
   fprintf (fp, "]';\n");
  }

 for(i=1;i<=NF;i++) {
    fprintf (fp, "Uz(%d,:)=[\n",i);
    for(j=1;j<=NP;j++){
      fprintf (fp, "%g\n", Uz(j,i-1));
    }
   fprintf (fp, "]';\n");
  }

 

  fclose(fp);
}

/* ========================================================================= */

void Read_ic (void)
{ 



  /* Reads initial condition file (ic.dat) if it exists, otherwise calls
   * Generate_ic() to generate new initial condition. */

  double ux_in, uy_in, uz_in;
  double gamma_in, delta_in, m0_in;

  Real   *u_tmp, rp;
  int    np_ic, nf_ic, npts_ic, npts_ic3, index, i_tmp, j_tmp, 
         i, j;
  char   f_type, dummy;
  FILE  *fp;

 

#if !COMPUTE 
   
  for(j=0;j<=NF-1;j++) {
   while( (dummy=getc(history_file_dat)) != '='); 
   fscanf(history_file_dat,"%lg, %lg, %lg",&gamma_in,&delta_in,&m0_in);
   gamma_param[j] = gamma_in;
   delta_0_bar_param[j] = delta_in;
   m_0_param[j] = m0_in;
   NEXT_LINE(history_file_dat);
   }    
                 
   /* Skip to 10th line  */
   for(i=0;i<7;i++) NEXT_LINE(history_file_dat); 
   Read_snapshot ();

#else

  if((fp=fopen("ic.dat","r"))==NULL) { 
    Generate_ic();
    return;
  }

  /* Read np_ic etc following = sign on second line of file */
  NEXT_LINE(fp);  

  for(j=0;j<=NF-1;j++) {
   while( (dummy=getc(fp)) != '='); 
   fscanf(fp,"%lg, %lg, %lg",&gamma_in,&delta_in,&m0_in);
   gamma_param[j] = gamma_in;
   delta_0_bar_param[j] = delta_in;
   m_0_param[j] = m0_in;
   NEXT_LINE(fp);
   }    
                 
  while( (dummy=getc(fp)) != '=');                   
  fscanf(fp, "%d, %d", &np_ic, &nf_ic);  

   if(nf_ic != NF)  printf(" Warning : nf_ic != NF");  

  /* Skip to 10th line and read first character to determine type 
     B(inary) or A(scii) */
  for(i=0;i<8;i++) NEXT_LINE(fp); 
  f_type = getc(fp); NEXT_LINE(fp); 
  
  if ( (f_type !='B') && (f_type !='A') ) {
    if(verbose) 
      printf("\n ic.dat exists but of unrecognized type Binary or Ascii \n"); 
    exit(1);
  }

  if(verbose) printf("\nReading ic.dat with np, nf = %d, %d... \n\n", 
		     np_ic, nf_ic);

  npts_ic  = (np_ic+2) * nf_ic;
  npts_ic3  = (np_ic+2) * 3 *  nf_ic;  
 

  /* Allocate temporary memory and read from file */


  u_tmp =(Real *) malloc((unsigned)(npts_ic3)*sizeof(Real));
 

  if(f_type =='B') {    
    /* Binary data file */
    fread(u_tmp, sizeof(Real), npts_ic3, fp);
  }
  else {          
    /* Ascii data file */
    for(index=0;index<npts_ic;index++) {
      fscanf (fp, "%lg %lg %lg\n", &ux_in, &uy_in, &uz_in);
     
      /*  printf("%lg %lg %lg\n", ux_in, uy_in, uz_in);*/
      
  
      u_tmp[index*3] = ux_in;
      u_tmp[index*3+1] = uy_in;
      u_tmp[index*3+2] = uz_in;
    }
  }

  /* Copy into u and v */
  rp = (np_ic+2.-1.0)/(NP+2.-1.0);
 
      

  for(j=1;j<=NF;j++) {
        j_tmp = 3*(NP+2)*(j-1);
      for(i=1;i<=NP+2;i++) {
	i_tmp = (int)(rp*(i-1));
	Ux(i-1,j-1) = u_tmp[i_tmp * 3 + j_tmp];
        Uy(i-1,j-1) = u_tmp[i_tmp * 3 + j_tmp + 1];
	Uz(i-1,j-1) = u_tmp[i_tmp * 3 + j_tmp + 2];
     /*  printf("%lg %lg %lg\n", Ux(i,j), Uy(i,j), Uz(i,j)); */
      }
  }


  free(u_tmp);
#endif

}
/* ========================================================================= */

void Generate_ic (void)
{
  /* Generate new initial condition of type ic_type.  Some cases are given,
   * you might want to add more. */

  Real radius, amplitude, ri_o_ro, k_tmp;
  Real s_step =  2*M_PI/(NP-1);   /*   space angular step  */
  int i, mode;
  Real a_axe, b_axe, center_x, center_y, center_z, length_x, length_y;
  Real wave_length, theta;

  if(verbose) 
    printf("\nGenerating initial condition of ic_type %d. \n\n", ic_type);


   

  if(ic_type>0) {

    /* ------------------------------------------ 
     *                   Compute ic 
     * ------------------------------------------ */

	  /* Compute different Initial Conditions*/
	  switch (ic_type) {
	  case 1:  
	    /* ellipseyz */

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

            a_axe =1., b_axe = 0.75, center_x = 0, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
          }

          if (NF!=1){
              printf("In this case 1 the number of filaments nf has to be 1 !!!\n");
              exit(0);
          }

	    break;

	  case 2:  
	    /* standing waves on a vortex ring */

                radius =(Real)1., amplitude =(Real)0.25, mode =4;

                gamma_param[0] = 1.;
                m_0_param[0] = 0.;
                delta_0_bar_param[0] = 1.;


          for(i=0;i<=NP+1;i++) {
          

 	  Ux(i,0) = 0;
	  Uy(i,0) = (1 + amplitude * cos(mode * (i-1) * s_step)) * cos((i-1) * s_step);
	  Uy(i,0) = radius * Uy(i,0); 
	  Uz(i,0) = (1 + amplitude * cos(mode * (i-1) * s_step)) * sin((i-1) * s_step);
          Uz(i,0) = radius * Uz(i,0);

 	  Ux(i,0) += radius * amplitude * sin(mode * (i-1) * s_step);
	 
          }

   
          if (NF!=1){
              printf("In this case 2 the number of filaments nf has to be 1 !!!\n");
              exit(0);
          }


	    break;

	  case 3:
	    /*  triangle   */

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

	    radius = 0.5;

	  for(i=0;i<=NP+1;i++) {
          
	  Ux(i,0) = 0;
	  Uy(i,0) = 3 * cos((i-1) * s_step) + cos(2 * (i-1) * s_step);
	  Uy(i,0) = radius * Uy(i,0); 
	  Uz(i,0) = 3 * sin((i-1) * s_step) - sin(2 * (i-1) * s_step);
          Uz(i,0) = radius * Uz(i,0); 
          }

          if (NF!=1){
              printf("In this case 3 the number of filaments nf has to be 1 !!!\n");
              exit(0);
          }

	    break;

	  case 4:
	    /* lissa*/

	  radius = 0.5;

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;
	    
	  for(i=0;i<=NP+1;i++) {
          
	  Ux(i,0) = sin(2 * (i-1) * s_step - M_PI/4);
	  Ux(i,0) = radius * Ux(i,0);
	  Uy(i,0) = 2. * cos((i-1) * s_step);
	  Uy(i,0) = radius * Uy(i,0); 
	  Uz(i,0) = 1.5 * sin((i-1) * s_step);
          Uz(i,0) = radius * Uz(i,0); 
          }

          if (NF!=1){
              printf("In this case 4 the number of filaments nf has to be 1 !!!\n");
              exit(0);
          }


	    break;


	  case 5:
	    /* side_b_side */

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

            a_axe =0.5, b_axe = 0.5, center_x = 0, center_y = -1, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
          }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;

           a_axe = 0.5, b_axe = 0.5, center_x = 0, center_y = 1, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Uz(i,1) = sin((i-1) * s_step);
          Uz(i,1) = b_axe* Uz(i,1) + center_z; 
          }

          if (NF!=2){
              printf("In this case 5 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }


	    break;


	  case 6:
	    /* jump */

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          a_axe =1, b_axe = 1, center_x = 0.25, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
          }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


           a_axe = 1, b_axe = 1, center_x = -0.25, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Uz(i,1) = sin((i-1) * s_step);
          Uz(i,1) = b_axe* Uz(i,1) + center_z; 
          }

          if (NF!=2){
              printf("In this case 6 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	    break;

	  case 7:
	    /* face_to_face*/

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

#if 1
      a_axe =1, b_axe = 1, center_x = -0.15, center_y = -0.9, center_z = 0;
#endif
#if 0
      a_axe =1, b_axe = 1, center_x = -0.15, center_y = -0.25, center_z = 0;
#endif

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
          }

          gamma_param[1] = -1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;

#if 1
           a_axe =1, b_axe = 1, center_x = 0.15, center_y = 0.9, center_z = 0;
#endif
#if 0
      a_axe =1, b_axe = 1, center_x = -0.15, center_y = -0.25, center_z = 0;
#endif

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Uz(i,1) = sin((i-1) * s_step);
          Uz(i,1) = b_axe* Uz(i,1) + center_z; 
          }

          if (NF!=2){
              printf("In this case 7 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;



         case 8:
	    /* face to face */

          gamma_param[0] = -1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          a_axe =1, b_axe = 1, center_x = 0.5, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
          }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


           a_axe = 1, b_axe = 1, center_x = -0.5, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Uz(i,1) = sin((i-1) * s_step);
          Uz(i,1) = b_axe* Uz(i,1) + center_z; 
          }

          if (NF!=2){
              printf("In this case 8 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	    break;

	  case 9:
	    /* link*/

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

           a_axe =1., b_axe = 1.;

           center_x = 0, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
	  }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;

          center_x = 0, center_y = -1., center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Uz(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Ux(i,1) = sin((i-1) * s_step);
          Ux(i,1) = b_axe* Ux(i,1) + center_z; 
	    }

          if (NF!=2){
              printf("In this case 9 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;

	  case 10:
	    /* inside*/

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

           a_axe =1.5, b_axe = 0.5;

           center_x = 0, center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
	  }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


	  a_axe =0.5, b_axe = 1.5;
          for(i=0;i<=NP+1;i++) {
         
	  Uz(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Ux(i,1) = sin((i-1) * s_step);
          Ux(i,1) = b_axe* Ux(i,1) + center_z; 
	  }

          if (NF!=2){
              printf("In this case 10 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;

          case 11:
	    /* One straight filament */

          wave_length = 1.25;
          amplitude = (Real)0.01;
          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          length_x = 1.25, length_y = 0.;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += amplitude * sin(Ux(i,0)*2*M_PI/wave_length);
	  Uz(i,0) = 0.;
	  }

     

          if (NF!=1){
              printf("In this case 11 the number of filaments nf has to be 1 !!!\n");
              exit(0);
          }

	    break;

          case 12:
	    /* One helical straight filament */

          wave_length = 1.25;
          amplitude = (Real)0.1;
          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          length_x = 7.5;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
          Uy(i,0) = amplitude * sin(Ux(i,0)*2*M_PI/wave_length);
	  Uz(i,0) = amplitude * cos(Ux(i,0)*2*M_PI/wave_length);
	  }

          if (NF!=1){
              printf("In this case 12 the number of filaments nf has to be 1 !!!\n");
              exit(0);
          }

	    break;

     
	  case 13:
	    /* Oscillation of two trailling vortices */

	       wave_length = 1.25; 
	     /*  wave_length = 7.5;*/
	      /*    wave_length = 7.0;*/

          amplitude = (Real)0.001;
          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

	  length_x = wave_length, length_y = 0.5;

	    /*    length_x = 42., length_y = 0.5;*/

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += amplitude * sin(Ux(i,0)*2*M_PI/wave_length);
	  Uz(i,0) = 0.;
	  }

          gamma_param[1] = -1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


	  for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,1) = -length_y ;
          Uy(i,1) -= amplitude * sin(Ux(i,1)*2*M_PI/wave_length);
	  Uz(i,1) = 0.;
	  }
	
          if (NF!=2){
              printf("In this case 13 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;

	  case 14:
	    /* Crow instability  of two trailling vortices*/

	      wave_length = 1.25; 
	      /*	       wave_length = 7.5;*/
	      /*    wave_length = 7.0;*/

          amplitude = (Real)0.05;
          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

	  length_x = 7.5, length_y = 0.5;
	    /*    length_x = 42., length_y = 0.5;*/

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += amplitude * sin(Ux(i,0)*2*M_PI/wave_length);
	  Uz(i,0) = 0.;
	  }

          gamma_param[1] = -1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


	  for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,1) = -length_y ;
          Uy(i,1) -= amplitude * sin(Ux(i,1)*2*M_PI/wave_length);
	  Uz(i,1) = 0.;
	  }
	
          if (NF!=2){
              printf("In this case 14 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;
	  case 15:
	    /* Crow instability  of two trailling vortices*/
	    /* Mode in a plane*/


                 /* conf1*/
	         /* wave_length = 4.8336;*/
                 /* theta = 47.61/180*M_PI;*/

                 /* conf2*/
	         /*   wave_length = 4.7653;*/
                 /*   theta = 47.50/180*M_PI;*/

                 /* conf3*/
	         /* wave_length = 4.7673;*/
                 /* theta = 47.61/180*M_PI;*/

                 /* conf6*/
	          wave_length = 5.8482;
                  theta = 47.39/180*M_PI;



	        /* epsilon = 0.02*/
	        /*  wave_length = 10.21;*/
                /*  theta = 47.40/180*M_PI;*/

	        /* epsilon = 0.04*/
	        /*   wave_length = 9.3;*/
                /*   theta = 47.36/180*M_PI;*/

	        /* epsilon = 0.07*/
	        /*  wave_length = 8.6;*/
                /*  theta = 47.6/180*M_PI;*/

	        /* epsilon = 0.1 */
	        /*   wave_length = 7.95; */
                /*   theta = 47.72/180*M_PI; */

	        /* epsilon = 0.12 */
	        /*   wave_length = 7.68;*/
                /*   theta = 47.63/180*M_PI;*/

	        /* epsilon = 0.15 */
	        /*  wave_length = 7.28;*/
                /*  theta = 47.63/180*M_PI;*/
  
                /* epsilon = 0.2 */
	        /*  wave_length = 6.684; */
                /*  theta = 47.74186/180*M_PI;*/
  
 
	         /*  theta = 30./180*M_PI;*/
	         /*  theta = 0.;*/


          amplitude = (Real)0.01;


           /* gamma_param[0] = 1.;*/
            m_0_param[0] = 0.;

	   /* conf1*/
           /*gamma_param[0] = 4.15;*/
           /*m_0_param[0] = 0.;*/

	   /* conf2*/
           /*gamma_param[0] = 4.23;*/
           /*m_0_param[0] = 0.;*/

	   /* conf3*/
           /*gamma_param[0] = 5.47;*/
           /*m_0_param[0] = 0.;*/


	   /* conf6*/
           gamma_param[0] = 2.29;
           /*m_0_param[0] = 0.;*/


          delta_0_bar_param[0] = 1.;


	   /* conf1*/
           /* length_y = 0.527/2;*/

	   /* conf2*/
           /*length_y = 0.516/2;*/

	   /* conf3*/
           /* length_y = 0.490/2;*/

	   /* conf6*/
            length_y = 0.540/2;


          /* length_x =7.95;*/



	       /*  length_x = 10.21;*/
	       /*  length_x = 9.3;*/
	       /*  length_x = 8.6;*/
	       /*  length_x = 7.95;*/
	      /*   length_x = 7.68;*/
	         length_x = wave_length;


          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,0) = 0.+amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*sin(theta);
	  }

	   /* conf1*/
           /* gamma_param[1] = -4.15;*/


	   /* conf2*/
           /* gamma_param[1] = -4.23;*/

	   /* conf3*/
           /* gamma_param[1] = -5.47;*/

	   /* conf6*/
            gamma_param[1] = -2.29;



          m_0_param[1] = -0.;
          delta_0_bar_param[1] = 1.;


	  for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,1) = -length_y ;
          Uy(i,1) -= amplitude * sin(Ux(i,1)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,1) = 0.+amplitude * sin(Ux(i,1)*2*M_PI/wave_length)*sin(theta);
	  }
	
          if (NF!=2){
              printf("In this case 15 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;

	  case 16:
	    /* Vortex ring pairs */

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          c_v[0] = C_V(delta_0_bar_param[0]);
          c_w[0] = 0.; /* it must be m_0_param[0] = 0. */
          delta_bar[0] = delta_0_bar_param[0]; 

          a_axe =1., b_axe = 1., center_x = 0., center_y = 0, center_z = 0;

          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = 0 +  center_x;
	  Uy(i,0) = cos((i-1) * s_step);
	  Uy(i,0) = a_axe * Uy(i,0) + center_y; 
	  Uz(i,0) = sin((i-1) * s_step);
          Uz(i,0) = b_axe* Uz(i,0) + center_z; 
          }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;

          c_v[1] = C_V(delta_0_bar_param[1]);
          c_w[1] = 0.; /* it must be m_0_param[0] = 0. */
          delta_bar[1] = delta_0_bar_param[1]; 

           a_axe = 2., b_axe = 2., center_x = -0., center_y = 0, center_z = 0;

           ri_o_ro = 1. / a_axe; 
           k_tmp = 2.*sqrt(ri_o_ro)/(1+ri_o_ro);
           gamma_param[1] = - gamma_param[0]* (EllipticE(k_tmp)/(1-ri_o_ro)-EllipticK(k_tmp)/(1+ri_o_ro)+0.5/ri_o_ro*(log(8.*1./delta_bar[0]/epsilon)-1+c_v[0]+c_w[0]))/(EllipticE(k_tmp)/(1-ri_o_ro)+EllipticK(k_tmp)/(1+ri_o_ro)-0.5*(log(8.*a_axe/delta_bar[1]/epsilon)-1+c_v[1]+c_w[1]));

 

    


          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = 0 +  center_x;
	  Uy(i,1) = cos((i-1) * s_step);
	  Uy(i,1) = a_axe * Uy(i,1) + center_y; 
	  Uz(i,1) = sin((i-1) * s_step);
          Uz(i,1) = b_axe* Uz(i,1) + center_z; 
          }

          if (NF!=2){
              printf("In this case 16 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	    break;

	  case 17:
	    /* Crow instability  of a trailling vortex pairs*/
	    /* Mode in a plane*/



	        /* epsilon = 0.02*/
	        /* wave_length = 10.21;*/
                /* theta = 47.40/180*M_PI;*/

	        /* epsilon = 0.04*/
	        /* wave_length = 9.3;*/
                /* theta = 47.36/180*M_PI;*/

	        /* epsilon = 0.07*/
	        /*  wave_length = 8.6;*/
                /*  theta = 47.6/180*M_PI;*/

	        /* epsilon = 0.1 */
	        /*  wave_length = 7.95;*/
                /*  theta = 47.72/180*M_PI;*/

	        /* epsilon = 0.12 */
	        /*  wave_length = 7.68;*/
                /*  theta = 47.63/180*M_PI;*/

	        /* epsilon = 0.15 */
	        /*   wave_length = 7.85;*/


	         /* case a) */
                 /*  wave_length = 0.89759790102566; */
	         /* case b) */
                 /*  wave_length = 7.85; */
	         /* case c) */
                 /*  wave_length = 7.85;*/
	         /* case a1) */
                 /*  wave_length = 1.25663706143592; */
	         /* case b1) */
                 /*  wave_length = 7.85398163397448;*/
	         /* case c1) */
                 /*  wave_length = 7.85398163397448;*/

                 /* case b2) */
                   wave_length = 10.21017;
                
  
		   /* case a) */
		   /* theta = (105.8607737185099)/180*M_PI; */
		    /*  theta = (111.2362107577465)/180*M_PI;*/
                  
	           /* case b) */
		   /*   theta = (145.4555560967843)/180*M_PI; */
                   /*     theta = (145.6842984862458)/180*M_PI;  */

	           /* case c) */
		    /* theta = (116.9034734253750)/180*M_PI; */
                    /* theta = (118.6527720619110)/180*M_PI; */ 

		   /* case a1) */
		   /* theta = (262.8168308070007-180)/180*M_PI; */
		    /*  theta = (83.07994602308850)/180*M_PI;*/
                  
	           /* case b1) */
		   /*   theta = (140.3662170412589)/180*M_PI; */
                     /*  theta = (145.6842984862458)/180*M_PI; */ 

	           /* case c1) */
		   /*  theta = (110.1343841500866)/180*M_PI;  */
                    /* theta = (118.6527720619110)/180*M_PI; */ 

		   /* case b2) */
		   /*   theta = (145.4555560967843)/180*M_PI; */
                     /*  theta = (146.4349)/180*M_PI;   */
                       theta = (147.1176)/180*M_PI;  

	   /* amplitude = (Real)0.001;*/
              amplitude = (Real)0.02;
    	   /* case a) */
	   /* amplitude = amplitude/57.45853169821793;*/
            /*  amplitude = amplitude/52.74925546273909;*/


           /* case b) */
	    /* amplitude = amplitude/9.72309249711754; */ 
            /*    amplitude = amplitude/9.79650359927358;   */   

   	   /* case c) */
	   /*  amplitude = amplitude/9.58319756964557;*/ 
            /*   amplitude = amplitude/9.73257164722301;*/ 

   	   /* case a1) */
	   /*  amplitude = amplitude/48.50135675572913;*/
            /*   amplitude = amplitude/48.51724900265352;*/


           /* case b1) */
	    /* amplitude = amplitude/10.00604144059512; */ 
            /*    amplitude = amplitude/10.08441274204574; */ 

   	   /* case c1) */
	    /* amplitude = amplitude/9.35090400316568;*/
            /*   amplitude = amplitude/9.73257164722301;*/ 

             /* case b2) */
	    /* amplitude = amplitude/3.7; */ 
                amplitude = amplitude/3.7699;  

          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          length_y = 0.5;
          /* length_x =7.95;*/



	       /*  length_x = 10.21;*/
	       /*  length_x = 9.3;*/
	       /*  length_x = 8.6;*/
	       /*  length_x = 7.95;*/
	      /*   length_x = 7.68;*/
	         length_x = wave_length;

#if 1  /* sym */
          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += +amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,0) = 0.+amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*sin(theta);
	  }
#endif 
#if 0  /* antisym */
          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += -amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,0) = 0.-amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*sin(theta);
	  }
#endif 


          gamma_param[1] = -1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


	  for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,1) = -length_y ;
          Uy(i,1) -= amplitude * sin(Ux(i,1)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,1) = 0.+amplitude * sin(Ux(i,1)*2*M_PI/wave_length)*sin(theta);
	  }

          
          /* length_y = length_y * 0.14;*/
          /* ri_o_ro = 0.14;*/

           length_y = length_y * 0.3;
           ri_o_ro = 0.3;

	   /* case a) */
	   /*    theta = (131.2450250513009)/180*M_PI;  */ 
           /*  theta = (130.2267477637575)/180*M_PI;   */  


	   /* case b) */
	   /* theta = (103.8504006964758)/180*M_PI;   */
           /*   theta = (103.7325921533532)/180*M_PI;  */   


           /* case c) */
	   /* theta = (167.0280924219493)/180*M_PI;  */
            /*   theta = (166.3874383186925)/180*M_PI;   */ 

	   /* case a1) */
	     /*  theta = (132.5335848090105)/180*M_PI;  */ 
            /* theta = (131.8618979964265)/180*M_PI;     */ 


	   /* case b1) */
	    /* theta = (104.3541493925483)/180*M_PI;   */
             /* theta = (104.3541493925483)/180*M_PI;   */  


           /* case c1) */
	   /* theta = (167.5446761870687)/180*M_PI;  */
            /*   theta = (166.3874383186925)/180*M_PI;   */ 

           /* case b2) */
	   /* theta = (103.8504006964758)/180*M_PI;   */
            /*   theta = (132.3262)/180*M_PI;     */
                  theta = (134.8293)/180*M_PI;

	      /* case a) */
	       /*    amplitude = amplitude*57.45853169821793;*/
	       /* amplitude = amplitude*52.74925546273909;*/


	      /* case b) */
	      /*  amplitude = amplitude*9.72309249711754;*/
	      /*  amplitude = amplitude*9.79650359927358;*/


	      /* case c) */
	      /*  amplitude = amplitude*9.58319756964557;*/
	      /* amplitude = amplitude*9.73257164722301;*/


	      /* case a1) */
	        /*   amplitude = amplitude*48.50135675572913;*/
	        /*  amplitude = amplitude*48.51724900265352;*/


	      /* case b1) */
	      /* amplitude = amplitude*10.00604144059512;*/
	     /*  amplitude = amplitude*10.08441274204574;*/


	      /* case c1) */
	      /*  amplitude = amplitude*9.35090400316568; */
	      /* amplitude = amplitude*9.73257164722301;*/

	      /* case b) */
	      /*  amplitude = amplitude*3.7;*/
	        amplitude = amplitude*3.7699;
 

          gamma_param[2] = - gamma_param[0] /(1+3*ri_o_ro*ri_o_ro)
                                   *(3*ri_o_ro+ri_o_ro*ri_o_ro*ri_o_ro);
          gamma_param[2] = - 0.3;
          
          m_0_param[2] = 0.;
          delta_0_bar_param[2] = 0.5;

#if 1  /* sym */
          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,2) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,2) = length_y ;
          Uy(i,2) += +amplitude * sin(Ux(i,2)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,2) = 0.+amplitude * sin(Ux(i,2)*2*M_PI/wave_length)*sin(theta);
	  }
#endif
#if 0  /* antisym */
          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,2) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,2) = length_y ;
          Uy(i,2) += -amplitude * sin(Ux(i,2)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,2) = 0.-amplitude * sin(Ux(i,2)*2*M_PI/wave_length)*sin(theta);
	  }
#endif


          gamma_param[3] = -gamma_param[2];
          m_0_param[3] = 0.;
          delta_0_bar_param[3] = delta_0_bar_param[2];


	  for(i=0;i<=NP+1;i++) {
         
	  Ux(i,3) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,3) = -length_y ;
          Uy(i,3) -= amplitude * sin(Ux(i,3)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,3) = 0.+amplitude * sin(Ux(i,3)*2*M_PI/wave_length)*sin(theta);
	  }
		
          if (NF!=4){
              printf("In this case 17 the number of filaments nf has to be 4 !!!\n");
              exit(0);
          }

	  break;
 

	  case 18:
	    /* Rotating  trailling vortices*/
	    /* Mode in a plane*/



	        /* epsilon = 0.02*/
	        /* wave_length = 10.21;*/
                /* theta = 47.40/180*M_PI;*/

	        /* epsilon = 0.04*/
	        /* wave_length = 9.3;*/
                /* theta = 47.36/180*M_PI;*/

	        /* epsilon = 0.07*/
	        /*  wave_length = 8.6;*/
                /*  theta = 47.6/180*M_PI;*/

	        /* epsilon = 0.1 */
	        /*  wave_length = 7.95;*/
                /*  theta = 47.72/180*M_PI;*/

	        /* epsilon = 0.12 */
	        /*  wave_length = 7.68;*/
                /*  theta = 47.63/180*M_PI;*/

	        /* epsilon = 0.15 */
	         wave_length = 7.28;
                 theta = 47.63/180*M_PI;
  
 
	         /*  theta = 30./180*M_PI;*/
	         /*  theta = 0.;*/


          amplitude = (Real)0.001;
          gamma_param[0] = 1.;
          m_0_param[0] = 0.;
          delta_0_bar_param[0] = 1.;

          length_y = 0.5;
          /* length_x =7.95;*/



	       /*  length_x = 10.21;*/
	       /*  length_x = 9.3;*/
	       /*  length_x = 8.6;*/
	       /*  length_x = 7.95;*/
	      /*   length_x = 7.68;*/
	         length_x = wave_length;


          for(i=0;i<=NP+1;i++) {
         
	  Ux(i,0) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,0) = length_y ;
          Uy(i,0) += amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,0) = 0.+amplitude * sin(Ux(i,0)*2*M_PI/wave_length)*sin(theta);
	  }

          gamma_param[1] = 1.;
          m_0_param[1] = 0.;
          delta_0_bar_param[1] = 1.;


	  for(i=0;i<=NP+1;i++) {
         
	  Ux(i,1) = (i-(NP+1)/2.) * length_x /(NP-1);
	  Uy(i,1) = -length_y ;
          Uy(i,1) -= amplitude * sin(Ux(i,1)*2*M_PI/wave_length)*cos(theta);
	  Uz(i,1) = 0.+amplitude * sin(Ux(i,1)*2*M_PI/wave_length)*sin(theta);
	  }
	
          if (NF!=2){
              printf("In this case 15 the number of filaments nf has to be 2 !!!\n");
              exit(0);
          }

	  break;


	  default:
	    /* good idea to set some default for out-of-range ic_type */
	    if(verbose) printf("ic_type out of range\n");
	    exit(1);
	  }

	  /* Map phi: C^1 -> R^k, i.e. from (p_i, p_r) to (u,v).  Only
           * the phase, theta, is needed.  */


	
      
    
  }

  else {

    /* ------------------------------------------ 
     * Compute ic for a perturbed vortex ring
     *  
     * ------------------------------------------ */

       radius =(Real)1., amplitude =(Real)0.01, mode =3;

    amplitude =(Real)0.1;

    gamma_param[0] = 1.;
    m_0_param[0] = 0.;
    delta_0_bar_param[0] = 1.;


    for(i=0;i<=NP+1;i++) {
          
	  Ux(i,0) = 0;
	  Uy(i,0) = (1 + amplitude * cos(mode * (i-1) * s_step)) * cos((i-1) * s_step);
	  Uy(i,0) = radius * Uy(i,0); 
	  Uz(i,0) = (1 + amplitude * cos(mode * (i-1) * s_step)) * sin((i-1) * s_step);
          Uz(i,0) = radius * Uz(i,0); 
     }

    if (NF!=1){
       printf("In this case 0 the number of filaments nf has to be 1 !!!\n");
       exit(0);
    }


  }

 

}


/* ========================================================================= */

static int   factorial(int n){
  int val_tmp, ii;

    val_tmp = 1;
    for(ii=1;ii<=n;ii++) {
      val_tmp = val_tmp * ii;
    }
   return val_tmp;
}
/* ========================================================================= */

/* ========================================================================= */

static Real   EllipticE(Real k){
  /* Complete elliptic integrals of the second kind */
  /* Algorithm by J. R. Herndon */
  /* Source www.netlib.org */

  Real val_tmp, t;

   t = 1-k * k;
   val_tmp = ((0.040905094*t+0.085099193)*t+0.44479204)*t+1.0-(((0.01382999*t+0.08150224)*t+0.24969795)*t)*log(t);    
   return val_tmp;
}
/* ========================================================================= */

/* ========================================================================= */

static Real   EllipticK(Real k){
  /*  Complete elliptic integral of the first kind */
  /* Algorithm by J. R. Herndon */
  /* Source www.netlib.org */
  Real val_tmp, t;

   t = 1-k * k;
   val_tmp =   ((0.032024666*t+0.054555509)*t+0.097932891)*t+1.3862944-(((0.010944912*t+0.060118519)*t+0.12475074)*t+0.5)*log(t); 
   return val_tmp;
}
/* ========================================================================= */
