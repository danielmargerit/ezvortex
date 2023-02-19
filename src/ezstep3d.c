/* ------------------------------------------------------------------------- 
 * ezstep3d.c -- Time Stepping Routines 
 *
 * Copyright (C) 2000 Daniel Margerit
 *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 15:03:55 $
 * ------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ezvortex.h"
#include "ezstep3d.h"





/* -------------------------------------------------------------------------
 * I use many pre-processor macros to make the time-stepping flexible.
 * See ezstep3d.h for the definitions.
 * ------------------------------------------------------------------------- */

static void  Average_periodic_directions (void);
static void  Impose_boundary_conditions (void);
static int   s1, s2;

const Real zero       = 0.;
const Real half       = 0.5;
const Real third      = 1./3.;
const Real one        = 1.;
const Real two        = 2.;
const Real six        = 6.;
const Real twentyfour = 24.;

static  Real one_o_ds2, one_o_2ds;

#if M1_KNIO_KLEIN
  const Real c_ttm = -0.4202;/* to be use for the M1 Knio/Klein method*/
#endif

#if   ADAMS_BASHFORTH
   const Real adam_a = 23./12.;
   const Real adam_b = -16./12.;
   const Real adam_c = 5./12.;
#endif


/* ========================================================================= */

void Step (FILE *fpp)
{
 
 


 Real omega = 400., error = 100., error_tmp, S_tmp; 
 Real radius, amplitude;
 int mode, counter = 0;
 Real s_step =  2*M_PI/(NP-1);
 Real point_m[3], point_c[3], point_c0[3], point_p[3];
 Real vec_tmp[3], vec_tmp1[3], vec_tmp2[3];
 Real vec_kb[3], vec_local[3];  
 Real vec_inter[3], vec_auto[3], vec_add[3];
 Real sigma_tmp, one_o_sigma3, one_o_norm3;
 Real S0_o_S, gamma_tmp;
 Real delta_0_bar_param_4, coeff_tmp, one_nu_tmp, S0_o_S_4;
 Real cst, cst1;
 int in, im;
 


#if SPECTRAL
   Real *work, sum_tmp; /* Fourier int components */
   int inc, jump, lot, isgn;
   Real *u_spect, *u_s_spect, *u_ss_spect;
   register int k,jj;
#endif
 

 


#if CALL_AND_TING 
 Real *lambda_u, *ux_tmp, *uy_tmp,*uz_tmp, *ux_s, *uy_s, *uz_s, *sigma_u;
 Real  lambda_tmp;
#endif

#if DE_SINGU
 Real cut_off_length, cut_off_length_tmp, cut_off_length_tmp2;
 Real *ux_tmp, *uy_tmp,  *uz_tmp, *ux_s, *uy_s, *uz_s, *sigma_u;/* if open */
#endif

#if M1_KNIO_KLEIN
  Real vec_auto1[3], vec_auto2[3]; 
  Real sigma_max, delta_ttm, sigma_1, sigma_2;
  Real delta_ttm_tmp, ln_o_ln_coeff;
  Real kernel_tanh_coeff_1, kernel_tanh_coeff_2;
  Real sigma3_1, sigma3_2, norm3;
  Real *ux_tmp, *uy_tmp,  *uz_tmp, *ux_s, *uy_s, *uz_s, *sigma_u;/* if open */                             
#endif
                                 
#if !CLOSED
  Real left_right[3], left_right_b[3];
  int flag_overflow;
#endif

 register int i, j, j_k, i_k, i_j, i_b;


#if SPECTRAL
   work = (Real *) calloc(3*(NP+1), sizeof(Real));
   u_spect = (Real *) calloc(3*(NP+2), sizeof(Real));
   u_s_spect =  (Real *) calloc(3*(NP+2), sizeof(Real));
   u_ss_spect =  (Real *) calloc(3*(NP+2), sizeof(Real));
#endif
 
   
#if 1   /* real physical step */
          Impose_boundary_conditions ();

	  /* store U in U_0 to operate the Newton iterative stepping */

   

/* 
 * Beginning of explicite
 * ------------------------------------------------------- */
#if EXPLICIT   

           for(j=0;j<=NF-1;j++) {
   
             /* Find  beta(t) */ 
             beta       = gamma_param[j]/4./M_PI; 
#if 0    
	     beta       =1.;
#endif

         
#if STRUCTURE
#if  CALL_AND_TING 
	     beta       =  beta * (-log(epsilon)+log(S[j])-1+c_v[j]+c_w[j]);
#endif
#if LOCAL_INDUCTION
#if CLOSED
	     beta       =  beta * (-log(epsilon)+log(S[j])-1+c_v[j]+c_w[j]);
#else
	     beta       =  beta * (-log(epsilon)+log(2)-1+c_v[j]+c_w[j]);
#endif
#endif
#endif


#if  DE_SINGU
      cut_off_length = epsilon;
#if STRUCTURE
       cut_off_length = cut_off_length * exp(-c_v[j]-c_w[j]);
#endif
#endif

#if  M1_KNIO_KLEIN
       delta_ttm = epsilon;
       delta_ttm = delta_ttm * exp(c_ttm-c_v[j]+1-c_w[j]);
       sigma_max = 0.;
           for(i=1;i<=NP;i++) {
            sigma_max = max(sigma_max,SIGMA(i,j));
	   }
       sigma_max = ds * sigma_max;
       sigma_1 = 3 * sigma_max;  /* you can choose something else than 3 */
       sigma_2 = 2 * sigma_1;    /* you can choose something else than 2 */
       sigma3_1 = sigma_1 * sigma_1 * sigma_1;
       sigma3_2 = sigma_2 * sigma_2 * sigma_2;
       /* printf("%g ",sigma_max);*/ 
 
#endif
    

              /* geometry derivatives */              
#if !SPECTRAL

               FIND_DERIVATIVES_FINITE_DIFF
#else
             /*  find X_s, sigma and K b by spectral method */

#if !CLOSED
          left_right[0] = Ux(NP,j)-Ux(1,j);
          left_right[1] = Uy(NP,j)-Uy(1,j);
          left_right[2] = Uz(NP,j)-Uz(1,j);
           for(i=1;i<=NP-1;i++) {
             Ux(i,j) =  Ux(i,j)-(i-1)*ds*left_right[0]/2/M_PI;      
             Uy(i,j) =  Uy(i,j)-(i-1)*ds*left_right[1]/2/M_PI;  
             Uz(i,j) =  Uz(i,j)-(i-1)*ds*left_right[2]/2/M_PI; 
	     } 
#endif
	 
          
	   
               FIND_DERIVATIVES_SPECTRAL
              
#if !CLOSED
           for(i=1;i<=NP-1;i++) {
             Ux(i,j) =  Ux(i,j)+(i-1)*ds*left_right[0]/2/M_PI;      
             Uy(i,j) =  Uy(i,j)+(i-1)*ds*left_right[1]/2/M_PI;  
             Uz(i,j) =  Uz(i,j)+(i-1)*ds*left_right[2]/2/M_PI; 
             Ux_s(i,j) =  Ux_s(i,j)+left_right[0]/2/M_PI;      
             Uy_s(i,j) =  Uy_s(i,j)+left_right[1]/2/M_PI;  
             Uz_s(i,j) =  Uz_s(i,j)+left_right[2]/2/M_PI;
	     } 

#endif    
               FIND_KB_SPECTRAL  
#endif
  

       Impose_boundary_conditions (); /* added */

  

      
	    /* loops on a filament */

#if CLOSED
               for(i=1;i<=NP-1;i++) {
#else
               for(i=1;i<=NP;i++) {
#endif
                
             point_c0[X] = Ux(i,j); 
	     point_c0[Y] = Uy(i,j); 
	     point_c0[Z] = Uz(i,j);            
	     
             point_c[X] = Ux(i,j);
	     point_c[Y] = Uy(i,j);
	     point_c[Z] = Uz(i,j);           	             
          
             vec_kb[X] = Ux_ss(i,j);
             vec_kb[Y] = Uy_ss(i,j);
             vec_kb[Z] = Uz_ss(i,j);

	     /* Multiply by beta */
	     vec_local[X] = beta * vec_kb[X];
	     vec_local[Y] = beta * vec_kb[Y];
	     vec_local[Z] = beta * vec_kb[Z];

            /* Find  Q2(t) : global mutual induction term*/ 

             FIND_GLOBAL_MUTUAL_INDUCTION

          
  
            /* Find  Q(t): global auto induction term */ 

             
             FIND_GLOBAL_AUTO_INDUCTION
          
 
            /* Add different contributions */
             vec_add[X] = 0.;
	     vec_add[Y] = 0.;
	     vec_add[Z] = 0.;

	     ADD_TO_VEC;
             
    

#if 1        /* Add also Mutual interaction by Q2 */
	     
              vec_add[X] += vec_inter[X];
              vec_add[Y] += vec_inter[Y];
              vec_add[Z] += vec_inter[Z];
#endif
 
#if CONV_ANALYSE
        /* Analyse of convergence at t=0: velocity fiel in file 'look' */
        fprintf(fpp,"%g\t%g\t%g\n",vec_add[X],vec_add[Y],vec_add[Z]);
#endif

 

          /* printf("%g,%g,%g\n ",vec_add[X],vec_add[Y],vec_add[Z]);*/

             /* Multiply by dt */
	     vec_add[X] = dt * vec_add[X];
	     vec_add[Y] = dt * vec_add[Y];
	     vec_add[Z] = dt * vec_add[Z];
    
                       
             /* add  to previous time step value */
             ADD(vec_tmp,point_c0,vec_add);
              
#if GAUSS_SEIDEL
             Ux(i,j) =  vec_tmp[X];      
             Uy(i,j) =  vec_tmp[Y];  
             Uz(i,j) =  vec_tmp[Z]; 
                 
#else
               
             Ux_0(i,j) =  vec_tmp[X];      
             Uy_0(i,j) =  vec_tmp[Y]; 
             Uz_0(i,j) =  vec_tmp[Z]; 
                      
 
               
#endif          
 
	     }  /* end of loops on a filament */
        
    
#if !GAUSS_SEIDEL
	    /* loops on a filament */
#if CLOSED
             for(i=1;i<=NP-1;i++) {
#else
             for(i=1;i<=NP;i++) {
#endif
             Ux(i,j) =  Ux_0(i,j);      
             Uy(i,j) =  Uy_0(i,j);  
             Uz(i,j) =  Uz_0(i,j); 
	     } 
#endif

     

        /* Find the new length S(t) */
         S_tmp =0.;
         for(i=1;i<=NP-1;i++) {
         S_tmp +=SIGMA(i,j);
	 }
         S[j] = S_tmp * ds;

        /* Find the new stetching terms */
           S0_o_S = S_0[j]/S[j];

         int_S_o_S0[j] = int_S_o_S0[j] + dt * (1/ S0_o_S);  
        

           /* Find the new core thickness delta_bar(t) */
           delta_bar[j]      = delta_0_bar_param[j]*delta_0_bar_param[j]; 
           delta_bar[j]      += 4 * nu_bar_param * int_S_o_S0[j];
           delta_bar[j]  = sqrt(delta_bar[j] * S0_o_S);
  
           m_flux[j] = m_0_param[j] * S0_o_S * S0_o_S;

           /* Find  the new Cv(t) and Cw(t) */
           c_v[j]         = C_V(delta_bar[j]);
           c_w[j]         = C_W(delta_bar[j], m_flux[j],gamma_param[j]); 

#if NON_SIMIL_PART
      FIND_NON_SIMIL_PART
#endif




	   }    /* end of loops over filaments */
#if SPECTRAL
   free(work);
   free(u_spect);
   free(u_s_spect);
   free(u_ss_spect);
#endif
 
#endif   /* end of explicite */

/* 
 * End of explicite
 * ------------------------------------------------------- */

/* 
 * Beginning of Newton
 * ------------------------------------------------------- */
#if   NEWTON

           for(j=0;j<=NF-1;j++) {
             for(i=0;i<=NP+1;i++) {
             Ux_0(i,j) = Ux(i,j); 
	     Uy_0(i,j) = Uy(i,j); 
	     Uz_0(i,j) = Uz(i,j);

	     }
            int_S_o_S0_tmp[j] = int_S_o_S0[j];
	   }


	   while ((error > error_stop) && (counter <= 11)){
#if 0
       
       printf("counter = %d\n", counter);
#endif
           if (counter == 11) { 
                   printf("Level 11 reached in Newton solver\n"); 
                   exit(0);
	   }
         
	   error = 0.;


           for(j=0;j<=NF-1;j++) {
   
             /* Find  beta(t) */ 
             beta       = gamma_param[j]/4./M_PI; 
#if 0    
	     beta       =1.;
#endif

         
#if STRUCTURE
#if  CALL_AND_TING 
	     beta       =  beta * (-log(epsilon)+log(S[j])-1+c_v[j]+c_w[j]);
            
#endif
#if LOCAL_INDUCTION
#if CLOSED
	     beta       =  beta * (-log(epsilon)+log(S[j])-1+c_v[j]+c_w[j]);
#else
	     beta       =  beta * (-log(epsilon)+log(2)-1+c_v[j]+c_w[j]);
#endif
#endif
#endif


#if  DE_SINGU
      cut_off_length = epsilon;
#if STRUCTURE
       cut_off_length = cut_off_length * exp(-c_v[j]-c_w[j]);
#endif
#endif

#if  M1_KNIO_KLEIN
       delta_ttm = epsilon;
       delta_ttm = delta_ttm * exp(c_ttm-c_v[j]+1-c_w[j]);
       sigma_max = 0.;
           for(i=1;i<=NP;i++) {
            sigma_max = max(sigma_max,SIGMA(i,j));
	   }
       sigma_max = ds * sigma_max;
       sigma_1 = 3 * sigma_max;  /* you can choose something else than 3 */
       sigma_2 = 2 * sigma_1;    /* you can choose something else than 2 */
       sigma3_1 = sigma_1 * sigma_1 * sigma_1;
       sigma3_2 = sigma_2 * sigma_2 * sigma_2;
       /*printf("%g ",delta_ttm_tmp);*/
 
#endif


             
              /* geometry derivatives */              
#if !SPECTRAL

               FIND_DERIVATIVES_FINITE_DIFF 
      
#else
             /*  find X_s, sigma and K b by spectral method */
             
#if !CLOSED
          left_right[0] = Ux(NP,j)-Ux(1,j);
          left_right[1] = Uy(NP,j)-Uy(1,j);
          left_right[2] = Uz(NP,j)-Uz(1,j);
           for(i=1;i<=NP-1;i++) {
             Ux(i,j) =  Ux(i,j)-(i-1)*ds*left_right[0]/2/M_PI;      
             Uy(i,j) =  Uy(i,j)-(i-1)*ds*left_right[1]/2/M_PI;  
             Uz(i,j) =  Uz(i,j)-(i-1)*ds*left_right[2]/2/M_PI; 
	     } 
#endif
               FIND_DERIVATIVES_SPECTRAL
#if !CLOSED
           for(i=1;i<=NP-1;i++) {
             Ux(i,j) =  Ux(i,j)+(i-1)*ds*left_right[0]/2/M_PI;      
             Uy(i,j) =  Uy(i,j)+(i-1)*ds*left_right[1]/2/M_PI;  
             Uz(i,j) =  Uz(i,j)+(i-1)*ds*left_right[2]/2/M_PI; 
             Ux_s(i,j) =  Ux_s(i,j)+left_right[0]/2/M_PI;      
             Uy_s(i,j) =  Uy_s(i,j)+left_right[1]/2/M_PI;  
             Uz_s(i,j) =  Uz_s(i,j)+left_right[2]/2/M_PI;
	     } 

#endif    
               FIND_KB_SPECTRAL  
#endif

                  
	    /* loops on a filament */
#if CLOSED
             for(i=1;i<=NP-1;i++) {
             /* for(i=48;i<=48;i++) {*/
#else
             for(i=1;i<=NP;i++) {
#endif
             point_c0[X] = Ux_0(i,j); 
	     point_c0[Y] = Uy_0(i,j); 
	     point_c0[Z] = Uz_0(i,j);

             point_c[X] = Ux(i,j);
	     point_c[Y] = Uy(i,j);
	     point_c[Z] = Uz(i,j);

          
             vec_kb[X] = Ux_ss(i,j);
             vec_kb[Y] = Uy_ss(i,j);
             vec_kb[Z] = Uz_ss(i,j);

        


	     /* Multiply by beta */
   
	     vec_local[X] = beta * vec_kb[X];
	     vec_local[Y] = beta * vec_kb[Y];
	     vec_local[Z] = beta * vec_kb[Z];

            /* Find  Q2(t) : global mutual induction term*/ 

             FIND_GLOBAL_MUTUAL_INDUCTION



            /* Find  Q(t): global auto induction term */ 

             
             FIND_GLOBAL_AUTO_INDUCTION
          

            /* Add different contributions */
             vec_add[X] = 0.;
	     vec_add[Y] = 0.;
	     vec_add[Z] = 0.;

	     ADD_TO_VEC;
             
       

#if 1        /* Add also Mutual interaction by Q2 */
	     
              vec_add[X] += vec_inter[X];
              vec_add[Y] += vec_inter[Y];
              vec_add[Z] += vec_inter[Z];
#endif

#if CONV_ANALYSE
        /* Analyse of convergence at t=0: velocity fiel in file 'look' */
        fprintf(fpp,"%.17e\t%.17e\t%.17e\n",vec_add[X],vec_add[Y],vec_add[Z]);
#endif

 

          /* printf("%g,%g,%g\n ",vec_add[X],vec_add[Y],vec_add[Z]);*/

             /* Multiply by dt */
	     vec_add[X] = dt * vec_add[X];
	     vec_add[Y] = dt * vec_add[Y];
	     vec_add[Z] = dt * vec_add[Z];
    

             /* add  to previous time step value */
             ADD(vec_tmp,point_c0,vec_add);
             Ux(i,j) =  vec_tmp[X];      
             Uy(i,j) =  vec_tmp[Y];  
             Uz(i,j) =  vec_tmp[Z]; 

    


          /* find error for the implicite iteration */        
             SUB(vec_tmp2, vec_tmp, point_c);
             error_tmp = NORMV(vec_tmp2);
	     /*   printf("error = %g, i = %d\n",  error_tmp, i); */
             if (error_tmp > error)  {
                              error = error_tmp;
	     }

	     }  /* end of loops on a filament */

     

        /* Find the new length S(t) */
         S_tmp =0.;
         for(i=1;i<=NP-1;i++) {
         S_tmp +=SIGMA(i,j);
	 }
         S[j] = S_tmp * ds;

        /* Find the new stetching terms */
           S0_o_S = S_0[j]/S[j];

         int_S_o_S0[j] = int_S_o_S0_tmp[j] + dt* (1/ S0_o_S);
    

           /* Find the new core thickness delta_bar(t) */
           delta_bar[j]      = delta_0_bar_param[j]*delta_0_bar_param[j]; 
           delta_bar[j]      += 4 * nu_bar_param * int_S_o_S0[j];
           delta_bar[j]  = sqrt(delta_bar[j] * S0_o_S);
  
           m_flux[j] = m_0_param[j] * S0_o_S * S0_o_S;

           /* Find  the new Cv(t) and Cw(t) */
           c_v[j]         = C_V(delta_bar[j]);
           c_w[j]         = C_W(delta_bar[j], m_flux[j],gamma_param[j]); 

#if NON_SIMIL_PART
      FIND_NON_SIMIL_PART
#endif




	   }    /* end of loops over filaments */
         
           counter++;
           
           Impose_boundary_conditions ();

    }

#if SPECTRAL
   free(work);
   free(u_spect);
   free(u_s_spect);
   free(u_ss_spect);
#endif
 
#endif /* end of Newton */
/* 
 * End of Newton
 * ------------------------------------------------------- */


/* 
 * Beginning of Adam-Bashforth
 * ------------------------------------------------------- */
#if   ADAMS_BASHFORTH
            
           for(j=0;j<=NF-1;j++) {
   
             /* Find  beta(t) */ 
             beta       = gamma_param[j]/4./M_PI; 
#if 0    
	     beta       =1.;
#endif

         
#if STRUCTURE
#if  CALL_AND_TING 
	     beta       =  beta * (-log(epsilon)+log(S[j])-1+c_v[j]+c_w[j]);
#endif
#if LOCAL_INDUCTION
#if CLOSED
	     beta       =  beta * (-log(epsilon)+log(S[j])-1+c_v[j]+c_w[j]);
#else
	     beta       =  beta * (-log(epsilon)+log(2)-1+c_v[j]+c_w[j]);
#endif
#endif
#endif


#if  DE_SINGU
      cut_off_length = epsilon;
#if STRUCTURE
       cut_off_length = cut_off_length * exp(-c_v[j]-c_w[j]);
#endif
#endif

#if  M1_KNIO_KLEIN
       delta_ttm = epsilon;
       delta_ttm = delta_ttm * exp(c_ttm-c_v[j]+1-c_w[j]);
       sigma_max = 0.;
           for(i=1;i<=NP;i++) {
            sigma_max = max(sigma_max,SIGMA(i,j));
	   }
       sigma_max = ds * sigma_max;
       sigma_1 = 3 * sigma_max;  /* you can choose something else than 3 */
       sigma_2 = 2 * sigma_1;    /* you can choose something else than 2 */
       sigma3_1 = sigma_1 * sigma_1 * sigma_1;
       sigma3_2 = sigma_2 * sigma_2 * sigma_2;
       /*printf("tttttttttttttttttttttttttttttttt%g ",sigma_max);*/
       /*exit(0);*/
 
#endif

          


              /* geometry derivatives */              
#if !SPECTRAL

               FIND_DERIVATIVES_FINITE_DIFF  
#else
#if !CLOSED
          left_right[0] = Ux(NP,j)-Ux(1,j);
          left_right[1] = Uy(NP,j)-Uy(1,j);
          left_right[2] = Uz(NP,j)-Uz(1,j);
           for(i=1;i<=NP-1;i++) {
             Ux(i,j) =  Ux(i,j)-(i-1)*ds*left_right[0]/2/M_PI;      
             Uy(i,j) =  Uy(i,j)-(i-1)*ds*left_right[1]/2/M_PI;  
             Uz(i,j) =  Uz(i,j)-(i-1)*ds*left_right[2]/2/M_PI; 
	     } 
#endif
	   
               FIND_DERIVATIVES_SPECTRAL
        
           
#if !CLOSED
           for(i=1;i<=NP-1;i++) {
             Ux(i,j) =  Ux(i,j)+(i-1)*ds*left_right[0]/2/M_PI;      
             Uy(i,j) =  Uy(i,j)+(i-1)*ds*left_right[1]/2/M_PI;  
             Uz(i,j) =  Uz(i,j)+(i-1)*ds*left_right[2]/2/M_PI; 
             Ux_s(i,j) =  Ux_s(i,j)+left_right[0]/2/M_PI;      
             Uy_s(i,j) =  Uy_s(i,j)+left_right[1]/2/M_PI;  
             Uz_s(i,j) =  Uz_s(i,j)+left_right[2]/2/M_PI;
	     } 


#endif          
              FIND_KB_SPECTRAL  
#endif
          
          Impose_boundary_conditions (); /* added */
#if 0
                for(i=0;i<=NP+1;i++) {
       fprintf(fpp,"%d\t%g\t%g\t%g\n",i,Ux_ss(i,j),Uy_ss(i,j),Uz_ss(i,j));
	       }
#endif
#if 0
           for(i=1;i<=NP;i++) {
                  printf("%g\t",sqrt(Ux_ss(i,j)* Ux_ss(i,j)+Uy_ss(i,j)*Uy_ss(i,j)+Uz_ss(i,j)*Uz_ss(i,j)));
		}
                exit(0);
#endif

	    /* loops on a filament */
#if CLOSED
             for(i=1;i<=NP-1;i++) {
#else
               for(i=1;i<=NP;i++) {

             /*for(i=1;i<=1;i++) {*/
#endif

  
             point_c0[X] = Ux(i,j); 
	     point_c0[Y] = Uy(i,j); 
	     point_c0[Z] = Uz(i,j);

    	      /* printf("%g\t%g\t%g\n",point_c0[X],point_c0[Y],point_c0[Z]);*/
	     
             point_c[X] = Ux(i,j);
	     point_c[Y] = Uy(i,j);
	     point_c[Z] = Uz(i,j);
            	             


     
             vec_kb[X] = Ux_ss(i,j);
             vec_kb[Y] = Uy_ss(i,j);
             vec_kb[Z] = Uz_ss(i,j);  
               

	     /* Multiply by beta */
	     vec_local[X] = beta * vec_kb[X];
	     vec_local[Y] = beta * vec_kb[Y];
	     vec_local[Z] = beta * vec_kb[Z];

            /* Find  Q2(t) : global mutual induction term*/ 

             FIND_GLOBAL_MUTUAL_INDUCTION



            /* Find  Q(t): global auto induction term */ 

             
             FIND_GLOBAL_AUTO_INDUCTION
          

            /* Add different contributions */
             vec_add[X] = 0.;
	     vec_add[Y] = 0.;
	     vec_add[Z] = 0.;

	     ADD_TO_VEC;
             
       

#if 1        /* Add also Mutual interaction by Q2 */
	     
              vec_add[X] += vec_inter[X];
              vec_add[Y] += vec_inter[Y];
              vec_add[Z] += vec_inter[Z];
#endif

        
#if CONV_ANALYSE
        /* Analyse of convergence at t=0: velocity field in file 'look' */
        fprintf(fpp,"%.17e\t%.17e\t%.17e\n",vec_add[X],vec_add[Y],vec_add[Z]);
#endif

 

         /* store the velocity field */
             
             Vx(i,j) =  vec_add[X];      
             Vy(i,j) =  vec_add[Y];  
             Vz(i,j) =  vec_add[Z];


             if (istep==0)  {
             Vx_0(i,j) =  vec_add[X];      
             Vy_0(i,j) =  vec_add[Y];  
             Vz_0(i,j) =  vec_add[Z];
             Vx_m(i,j) =  vec_add[X];      
             Vy_m(i,j) =  vec_add[Y];  
             Vz_m(i,j) =  vec_add[Z];                       
	     }


             
             if (istep==1)  {
             /* Multiply by dt */
	     vec_add[X] = dt * (1.5 * Vx(i,j) - 0.5 * Vx_0(i,j));
	     vec_add[Y] = dt * (1.5 * Vy(i,j) - 0.5 * Vy_0(i,j));
	     vec_add[Z] = dt * (1.5 * Vz(i,j) - 0.5 * Vz_0(i,j));    
             /* add  to previous time step value */
             ADD(vec_tmp,point_c0,vec_add);
#if GAUSS_SEIDEL
             Ux(i,j) =  vec_tmp[X];      
             Uy(i,j) =  vec_tmp[Y];  
             Uz(i,j) =  vec_tmp[Z]; 
#else
             Ux_0(i,j) =  vec_tmp[X];      
             Uy_0(i,j) =  vec_tmp[Y];  
             Uz_0(i,j) =  vec_tmp[Z]; 
#endif
	     }
             else{
              /* Multiply by dt */
	     vec_add[X] = dt * (adam_a * Vx(i,j) + adam_b * Vx_0(i,j) + adam_c * Vx_m(i,j));
	     vec_add[Y] = dt * (adam_a * Vy(i,j) + adam_b * Vy_0(i,j) + adam_c * Vy_m(i,j));
	     vec_add[Z] = dt * (adam_a * Vz(i,j) + adam_b * Vz_0(i,j) + adam_c * Vz_m(i,j));

             /* add  to previous time step value */
             ADD(vec_tmp,point_c0,vec_add);
  
#if GAUSS_SEIDEL
             Ux(i,j) =  vec_tmp[X];      
             Uy(i,j) =  vec_tmp[Y];  
             Uz(i,j) =  vec_tmp[Z]; 
#else
             Ux_0(i,j) =  vec_tmp[X];      
             Uy_0(i,j) =  vec_tmp[Y];  
             Uz_0(i,j) =  vec_tmp[Z]; 
#endif
	     }

	     Vx_m(i,j) = Vx_0(i,j);
	     Vy_m(i,j) = Vy_0(i,j);
	     Vz_m(i,j) = Vz_0(i,j);

	     Vx_0(i,j) = Vx(i,j);
	     Vy_0(i,j) = Vy(i,j);
	     Vz_0(i,j) = Vz(i,j);

	     }  /* end of loops on a filament */

#if !GAUSS_SEIDEL
	    /* loops on a filament */
#if CLOSED
             for(i=1;i<=NP-1;i++) {
#else
             for(i=1;i<=NP;i++) {
#endif
             Ux(i,j) =  Ux_0(i,j);      
             Uy(i,j) =  Uy_0(i,j);  
             Uz(i,j) =  Uz_0(i,j); 
	     } 
#endif

        /* Find the new length S(t) */
         S_tmp =0.;
         for(i=1;i<=NP-1;i++) {
         S_tmp +=SIGMA(i,j);
	 }
         S[j] = S_tmp * ds;

        /* Find the new stetching terms */
           S0_o_S = S_0[j]/S[j];

         int_S_o_S0[j] = int_S_o_S0[j] + dt * (1/ S0_o_S);  
        

           /* Find the new core thickness delta_bar(t) */
           delta_bar[j]      = delta_0_bar_param[j]*delta_0_bar_param[j]; 
           delta_bar[j]      += 4 * nu_bar_param * int_S_o_S0[j];
           delta_bar[j]  = sqrt(delta_bar[j] * S0_o_S);
  
           m_flux[j] = m_0_param[j] * S0_o_S * S0_o_S;

           /* Find  the new Cv(t) and Cw(t) */
           c_v[j]         = C_V(delta_bar[j]);
           c_w[j]         = C_W(delta_bar[j], m_flux[j],gamma_param[j]); 

#if NON_SIMIL_PART
      FIND_NON_SIMIL_PART
#endif




	   }    /* end of loops over filaments */
#if SPECTRAL
   free(work);
   free(u_spect);
   free(u_s_spect);
   free(u_ss_spect);
#endif
 
#endif /* end of Adam-Bashforth */
/* 
 * End of Adam-Bashforth
 * ------------------------------------------------------- */

#endif  /* end of real physical step */

}    
/* ========================================================================= */

void Step_ini (void)
{
   one_o_ds2 = 1. / (ds * ds);
   one_o_2ds = 1. / (ds * 2);
}
/* ========================================================================= */

void Average_periodic_directions (void)
{
}
/* ========================================================================= */

/* Define W and Sigma_w macros.  w is generic field (u or v) */
#define W(i,j,k)          w[INDEX(i,j,k)]
#define Sigma_w(s,i,j,k)  sigma_w[(s)*FIELD_SIZE + INDEX(i,j,k)]

static void Impose_boundary_conditions (void)
{ 
          int j;
          Real left_right[3];
          

            for(j=0;j<=NF-1;j++) {
#if CLOSED         
	  Ux(0,j) = Ux(NP-1,j);   /* set fictive points */
	  Uy(0,j) = Uy(NP-1,j);
	  Uz(0,j) = Uz(NP-1,j);
 

          Ux(NP+1,j) = Ux(2,j);
	  Uy(NP+1,j) = Uy(2,j);
	  Uz(NP+1,j) = Uz(2,j);

    /*         printf("%g,%g,%g\n",Ux(1,j),Uy(1,j),Uz(1,j));*/
    /*       printf("%g,%g,%g\n",Ux(NP,j),Uy(NP,j),Uz(NP,j));*/

#if 1
	  Ux_s(0,j) = Ux_s(NP-1,j);  
	  Uy_s(0,j) = Uy_s(NP-1,j);
	  Uz_s(0,j) = Uz_s(NP-1,j);
          SIGMA(0,j) = SIGMA(NP-1,j);

          Ux_s(NP+1,j) = Ux_s(2,j);
	  Uy_s(NP+1,j) = Uy_s(2,j);
	  Uz_s(NP+1,j) = Uz_s(2,j);
	  SIGMA(NP+1,j) = SIGMA(2,j);

	  Ux_ss(0,j) = Ux_ss(NP-1,j);  
	  Uy_ss(0,j) = Uy_ss(NP-1,j);
	  Uz_ss(0,j) = Uz_ss(NP-1,j);

          Ux_ss(NP+1,j) = Ux_ss(2,j);
	  Uy_ss(NP+1,j) = Uy_ss(2,j);
	  Uz_ss(NP+1,j) = Uz_ss(2,j);
#endif

#if 1
	  Ux(NP,j) = Ux(1,j);   /* set periodic point */
	  Uy(NP,j) = Uy(1,j);
	  Uz(NP,j) = Uz(1,j);
 /*	  Uz(NP,j) = Uz(1,j)+1e-16;*/

#if 1
	  Ux_s(NP,j) = Ux_s(1,j);   
	  Uy_s(NP,j) = Uy_s(1,j);
	  Uz_s(NP,j) = Uz_s(1,j);
          SIGMA(NP,j)  = SIGMA(1,j);

	  Ux_ss(NP,j) = Ux_ss(1,j);   
	  Uy_ss(NP,j) = Uy_ss(1,j);
	  Uz_ss(NP,j) = Uz_ss(1,j);
#endif

#endif

      /*     printf("%g,%g,%g\n",Ux(1,j),Uy(1,j),Uz(1,j));*/
       /*    printf("%g,%g,%g\n",Ux(NP,j),Uy(NP,j),Uz(NP,j));*/

#endif
#if !CLOSED
          left_right[0] = Ux(NP,j)-Ux(1,j);
          left_right[1] = Uy(NP,j)-Uy(1,j);
          left_right[2] = Uz(NP,j)-Uz(1,j);
         /* set fictive points */

         /* find the left fictive point */

	  Ux(0,j) = Ux(NP-1,j)-left_right[0];   
	  Uy(0,j) = Uy(NP-1,j)-left_right[1];
	  Uz(0,j) = Uz(NP-1,j)-left_right[2];

         /* find the right fictive point */

          Ux(NP+1,j) = Ux(2,j)+left_right[0];
	  Uy(NP+1,j) = Uy(2,j)+left_right[1];
	  Uz(NP+1,j) = Uz(2,j)+left_right[2];

#if 1
	  Ux_s(0,j) = Ux_s(NP-1,j);  
	  Uy_s(0,j) = Uy_s(NP-1,j);
	  Uz_s(0,j) = Uz_s(NP-1,j);
          SIGMA(0,j)  = SIGMA(NP-1,j);

          Ux_s(NP+1,j) = Ux_s(2,j);
	  Uy_s(NP+1,j) = Uy_s(2,j);
	  Uz_s(NP+1,j) = Uz_s(2,j);
          SIGMA(NP+1,j)  = SIGMA(2,j);

	  Ux_ss(0,j) = Ux_ss(NP-1,j);  
	  Uy_ss(0,j) = Uy_ss(NP-1,j);
	  Uz_ss(0,j) = Uz_ss(NP-1,j);

          Ux_ss(NP+1,j) = Ux_ss(2,j);
	  Uy_ss(NP+1,j) = Uy_ss(2,j);
	  Uz_ss(NP+1,j) = Uz_ss(2,j);
#endif

#if 1
	  Ux(NP,j) = Ux(1,j)+left_right[0];   /* set periodic point */
	  Uy(NP,j) = Uy(1,j)+left_right[1];
	  Uz(NP,j) = Uz(1,j)+left_right[2];
#endif

#if 1
	  Ux_s(NP,j) = Ux_s(1,j);   
	  Uy_s(NP,j) = Uy_s(1,j);
	  Uz_s(NP,j) = Uz_s(1,j);
          SIGMA(NP,j)  = SIGMA(1,j);

	  Ux_ss(NP,j) = Ux_ss(1,j);   
	  Uy_ss(NP,j) = Uy_ss(1,j);
	  Uz_ss(NP,j) = Uz_ss(1,j);
#endif

#endif

	    }
}
/* ========================================================================= */

#undef W
#undef Sigma_w



