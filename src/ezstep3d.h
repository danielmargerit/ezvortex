/* ------------------------------------------------------------------------- 
 * ezstep3d.h -- Macros for EZSTEP
 *
 * Copyright (C) 2000 Daniel Margerit
 *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 16:43:23 $
 * ------------------------------------------------------------------------- */

#include "ezvortex.h"
#ifndef _EZSTEP3D_
#define _EZSTEP3D_

/* ---------------------------------------------------------------------
 * Define the macros C_V(delta_bar) and C_W(delta_bar, m_flux,gamma_param)  according to 
 * the core model  you wish to simulate. Gaussian and uniform core are given 
 * but you can add other core.    
 * --------------------------------------------------------------------- */

#if SIMIL_PART /* Gaussian core */
  #define C_V(delta_bar)  ( (1+M_GAMMA-log(2.))/2.-log((delta_bar)) )
  #define C_W(delta_bar, m_flux,gamma_param) ( -2. * ((m_flux)/(delta_bar)/(gamma_param)) * ((m_flux)/(delta_bar)/(gamma_param)) )
#endif

#if UNIFORM_CORE  /* Uniform core : to be coherent you have to choose nu_bar=0 in task.dat*/  
  #define C_V(delta_bar)  ( 1./2+1./4-log((delta_bar)) )
  #define C_W(delta_bar, m_flux,gamma_param) ( -4. * ((m_flux)/(delta_bar)/(gamma_param)) * ((m_flux)/(delta_bar)/(gamma_param)) )
#endif

#if NON_SIMIL_PART 
#    define FIND_NON_SIMIL_PART \
           delta_0_bar_param_4 = delta_0_bar_param[j] * delta_0_bar_param[j];\
             one_nu_tmp = delta_bar[j] * delta_bar[j]/delta_0_bar_param_4;\
             one_nu_tmp = one_nu_tmp/S0_o_S;\
             delta_0_bar_param_4 = delta_0_bar_param_4 * delta_0_bar_param_4;\
             coeff_tmp = 4* M_PI * M_PI/gamma_param[j]/gamma_param[j];\
             coeff_tmp = coeff_tmp * delta_0_bar_param_4;\
             S0_o_S_4 = S0_o_S * S0_o_S;\
	     S0_o_S_4 = S0_o_S_4 * S0_o_S_4;\
                                            \
                                            \
             for(in=0;in<3;in++) {\
	      /*  if (in||im){N_LAGUERRE_MODE*/\
	     for(im=0;im<3;im++) {\
               if ((in==0)&&(im==0)){;}\
               else{\
             c_v[j] += coeff_tmp * DD(in,j) * DD(im,j) * AA(in,im) /(in+im)*exp(-(in+im)*log(one_nu_tmp));\
             c_w[j] -= 4* coeff_tmp / (delta_bar[j] * delta_bar[j]) * S0_o_S_4 * CC(in,j) * CC(im,j) * AA(in,im) /(in+im)*exp(-(in+im)*log(one_nu_tmp));\
	     }\
	     }\
	     }
#endif


/* ------------------------------------------------------------------------- 
 *          In general you should not need to change anything below          
 * ------------------------------------------------------------------------- */

/* I rely on the compiler to take care of most optimization.  The only place
 * where I give the compiler help is with the Laplacian formulas because
 * none of the compiles do the obvious thing.  */



#define X 0                           /* Indices for vector components */
#define Y 1
#define Z 2


#define  VECTOR_INT(n)  (int *)malloc((unsigned)(n)*sizeof(int))


/* ------------------------------------- 
 * Defining the method used:         
 * Local induction, Callegari and Ting or De-singularization,... 
 * ------------------------------------- */

#if CLOSED
#    define FIND_GLOBAL_MUTUAL_INDUCTION \
             vec_inter[X] = 0.;\
	     vec_inter[Y] = 0.;\
	     vec_inter[Z] = 0.;\
          for(j_k=0;j_k<=NF-1;j_k++) {\
            if (j_k!=j){\
                  \
                gamma_tmp = gamma_param[j_k]/4./M_PI;\
                \
               \
                for(i_k=1;i_k<=NP-1;i_k++) {\
                  vec_tmp[X] = Ux_s(i_k,j_k);\
                  vec_tmp[Y] = Uy_s(i_k,j_k);\
                  vec_tmp[Z] = Uz_s(i_k,j_k);\
                \
                  vec_tmp1[X] = Ux(i_k,j_k);\
                  vec_tmp1[Y] = Uy(i_k,j_k);\
                  vec_tmp1[Z] = Uz(i_k,j_k);\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_inter[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_inter[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_inter[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
	     }\
	     }
#endif
#if !CLOSED
#    define FIND_GLOBAL_MUTUAL_INDUCTION \
             vec_inter[X] = 0.;\
	     vec_inter[Y] = 0.;\
	     vec_inter[Z] = 0.;\
          for(j_k=0;j_k<=NF-1;j_k++) {\
            if (j_k!=j){\
                  \
                gamma_tmp = gamma_param[j_k]/4./M_PI;\
                \
                left_right[0] = Ux(NP,j_k)-Ux(1,j_k);\
                left_right[1] = Uy(NP,j_k)-Uy(1,j_k);\
                left_right[2] = Uz(NP,j_k)-Uz(1,j_k);\
               \
                for(i_k=1;i_k<=NP;i_k++) {\
                 \
                  vec_tmp[X] = Ux_s(i_k,j_k);\
                  vec_tmp[Y] = Uy_s(i_k,j_k);\
                  vec_tmp[Z] = Uz_s(i_k,j_k);\
                \
                  vec_tmp1[X] = Ux(i_k,j_k);\
                  vec_tmp1[Y] = Uy(i_k,j_k);\
                  vec_tmp1[Z] = Uz(i_k,j_k);\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
               \
               /* proceed by the crossproduct*/\
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_inter[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_inter[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_inter[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
             \
               for(i_b=1;i_b<=n_b;i_b++) {\
                   left_right_b[0] = i_b * left_right[0];\
                   left_right_b[1] = i_b * left_right[1];\
                   left_right_b[2] = i_b * left_right[2];\
           /* find the the mutual induction of the left part*/\
               for(i_k=1;i_k<=NP-1;i_k++) {\
                 \
                  vec_tmp[X] = Ux_s(i_k,j_k);\
                  vec_tmp[Y] = Uy_s(i_k,j_k);\
                  vec_tmp[Z] = Uz_s(i_k,j_k);\
                \
                  vec_tmp1[X] = Ux(i_k,j_k)-left_right_b[0];\
                  vec_tmp1[Y] = Uy(i_k,j_k)-left_right_b[1];\
                  vec_tmp1[Z] = Uz(i_k,j_k)-left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
               \
               /* proceed by the crossproduct*/\
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_inter[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_inter[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_inter[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
             \
                /* find the the mutual induction of the right part*/\
             for(i_k=2;i_k<=NP;i_k++) {\
                 \
                  vec_tmp[X] = Ux_s(i_k,j_k);\
                  vec_tmp[Y] = Uy_s(i_k,j_k);\
                  vec_tmp[Z] = Uz_s(i_k,j_k);\
                \
                  vec_tmp1[X] = Ux(i_k,j_k)+left_right_b[0];\
                  vec_tmp1[Y] = Uy(i_k,j_k)+left_right_b[1];\
                  vec_tmp1[Z] = Uz(i_k,j_k)+left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
               /* proceed by the crossproduct*/\
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_inter[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_inter[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_inter[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
          	}\
	     }\
	     }
#endif   /* End of FIND_GLOBAL_MUTUAL_INDUCTION */


#if LOCAL_INDUCTION
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION \
             vec_auto[X] = 0.;\
	     vec_auto[Y] = 0.;\
	     vec_auto[Z] = 0.;
#endif
#endif /* End of if LOCAL_INDUCTION */

#if CALL_AND_TING 
#if CLOSED
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION  \
             vec_auto[X] = 0.;\
	     vec_auto[Y] = 0.;\
	     vec_auto[Z] = 0.;\
                              \
             ux_tmp  = VECTOR(NP+2);\
             uy_tmp  = VECTOR(NP+2);\
             uz_tmp  = VECTOR(NP+2);\
             ux_s = VECTOR(NP+2); \
             uy_s = VECTOR(NP+2); \
             uz_s = VECTOR(NP+2); \
             sigma_u = VECTOR(NP+2);\
             lambda_u = VECTOR(NP+2);\
                                        \
             for(i_j=1;i_j<=(NP+1)/2;i_j++) {\
             int ll = i_j+i-1;\
               if (ll>(NP-1)) {\
                     ll = ll-(NP-1);\
               }\
             ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j);\
             uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j);\
             uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j);\
             ux_s[(NP+1)/2+i_j-1] = Ux_s(ll,j);\
             uy_s[(NP+1)/2+i_j-1] = Uy_s(ll,j);\
             uz_s[(NP+1)/2+i_j-1] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2+i_j-1] =  SIGMA(ll,j);\
	     }\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
             int ll = i-i_j;\
               if (ll<1) {\
                     ll = (NP-1)+ll;\
               }\
             ux_tmp[(NP+1)/2-i_j] = Ux(ll,j);\
             uy_tmp[(NP+1)/2-i_j] = Uy(ll,j);\
             uz_tmp[(NP+1)/2-i_j] = Uz(ll,j);\
             ux_s[(NP+1)/2-i_j] = Ux_s(ll,j);\
             uy_s[(NP+1)/2-i_j] = Uy_s(ll,j);\
             uz_s[(NP+1)/2-i_j] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2-i_j] =  SIGMA(ll,j);\
	     }\
           /* find the left and right distance of points to point i  */\
             lambda_tmp = 0.;\
             cst = sigma_u[(NP+1)/2];\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
              cst1 = sigma_u[(NP+1)/2+i_j];\
              lambda_tmp += (cst+cst1)/2;\
              cst = cst1;\
              lambda_u[(NP+1)/2+i_j] = ds * lambda_tmp;\
		       }\
             lambda_tmp = 0.;\
             cst = sigma_u[(NP+1)/2];\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
              cst1 = sigma_u[(NP+1)/2-i_j];\
              lambda_tmp += (cst+cst1)/2;\
              cst = cst1;\
              lambda_u[(NP+1)/2-i_j] = ds * lambda_tmp;\
	      }\
          /* then proceed to the auto_induction calculus  */\
             gamma_tmp = gamma_param[j]/4./M_PI;\
                for(i_j=1;i_j<=NP-1;i_j++) {\
                 if (i_j!=(NP+1)/2){\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j];\
                  vec_tmp1[Y] = uy_tmp[i_j];\
                  vec_tmp1[Z] = uz_tmp[i_j];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_tmp1[X] = vec_kb[X]/2 * sigma_u[i_j]/lambda_u[i_j];\
                  vec_tmp1[Y] = vec_kb[Y]/2 * sigma_u[i_j]/lambda_u[i_j];\
                  vec_tmp1[Z] = vec_kb[Z]/2 * sigma_u[i_j]/lambda_u[i_j];\
                     \
                  SUB(vec_tmp2,vec_tmp2,vec_tmp1);\
		  /*if(i==i){fprintf(fpp,"%d\t%g\t%g\t%g\n",i_j,vec_tmp2[X],vec_tmp2[Y],vec_tmp2[Z]);}*/\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		 }\
		}\
             free(ux_tmp);\
             free(uy_tmp);\
             free(uz_tmp);\
             free(ux_s);\
             free(uy_s);\
             free(uz_s);\
             free(sigma_u);\
             free(lambda_u);
#endif
#endif

#if !CLOSED
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION  \
             vec_auto[X] = 0.;\
	     vec_auto[Y] = 0.;\
	     vec_auto[Z] = 0.;\
                              \
             ux_tmp  = VECTOR(NP+2);\
             uy_tmp  = VECTOR(NP+2);\
             uz_tmp  = VECTOR(NP+2);\
             ux_s = VECTOR(NP+2); \
             uy_s = VECTOR(NP+2); \
             uz_s = VECTOR(NP+2); \
             sigma_u = VECTOR(NP+2);\
             lambda_u = VECTOR(NP+2);\
                                      \
           /* find the translation vector*/\
                \
                left_right[0] = Ux(NP,j)-Ux(1,j);\
                left_right[1] = Uy(NP,j)-Uy(1,j);\
                left_right[2] = Uz(NP,j)-Uz(1,j);\
                flag_overflow = 0;\
               \
           /* construct the u_tmp vector */\
             for(i_j=1;i_j<=(NP+1)/2;i_j++) {\
             int ll = i_j+i-1;\
               if (ll>NP) {\
                     ll = ll-(NP-1);\
                     flag_overflow = 1;\
               }\
	       if (flag_overflow==0){\
                  ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j);\
                  uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j);\
                  uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j);\
               }\
               else{\
                  ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j)+left_right[0];\
                  uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j)+left_right[1];\
                  uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j)+left_right[2];\
	       }\
             ux_s[(NP+1)/2+i_j-1] = Ux_s(ll,j);\
             uy_s[(NP+1)/2+i_j-1] = Uy_s(ll,j);\
             uz_s[(NP+1)/2+i_j-1] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2+i_j-1] =  SIGMA(ll,j);\
	     }\
                                 \
             flag_overflow = 0;\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
             int ll = i-i_j;\
               if (ll<1) {\
                     ll = (NP-1)+ll;\
                   flag_overflow = 1;\
               }\
	       if (flag_overflow==0){\
                  ux_tmp[(NP+1)/2-i_j] = Ux(ll,j);\
                  uy_tmp[(NP+1)/2-i_j] = Uy(ll,j);\
                  uz_tmp[(NP+1)/2-i_j] = Uz(ll,j);\
               }\
               else{\
                  ux_tmp[(NP+1)/2-i_j] = Ux(ll,j)-left_right[0];\
                  uy_tmp[(NP+1)/2-i_j] = Uy(ll,j)-left_right[1];\
                  uz_tmp[(NP+1)/2-i_j] = Uz(ll,j)-left_right[2];\
	       }\
             ux_s[(NP+1)/2-i_j] = Ux_s(ll,j);\
             uy_s[(NP+1)/2-i_j] = Uy_s(ll,j);\
             uz_s[(NP+1)/2-i_j] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2-i_j] =  SIGMA(ll,j);\
	     }\
           /* find the left and right distance of points to point i  */\
             lambda_tmp = 0.;\
             cst = sigma_u[(NP+1)/2];\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
              cst1 = sigma_u[(NP+1)/2+i_j];\
              lambda_tmp += (cst+cst1)/2;\
              cst = cst1;\
              lambda_u[(NP+1)/2+i_j] = ds * lambda_tmp;\
		       }\
             lambda_tmp = 0.;\
             cst = sigma_u[(NP+1)/2];\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
              cst1 = sigma_u[(NP+1)/2-i_j];\
              lambda_tmp += (cst+cst1)/2;\
              lambda_u[(NP+1)/2-i_j] = ds * lambda_tmp;\
	      }\
              /* then proceed to the auto_induction calculus  */\
             gamma_tmp = gamma_param[j]/4./M_PI;\
               /* find the translation vector*/\
                \
                left_right[0] = ux_tmp[NP]-ux_tmp[1];\
                left_right[1] = uy_tmp[NP]-uy_tmp[1];\
                left_right[2] = uz_tmp[NP]-uz_tmp[1];\
               \
                /* find the autoinduction of the central part  */\
                for(i_j=1;i_j<=NP;i_j++) {\
                 if (i_j!=(NP+1)/2){\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j];\
                  vec_tmp1[Y] = uy_tmp[i_j];\
                  vec_tmp1[Z] = uz_tmp[i_j];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_tmp1[X] = vec_kb[X]/2 * sigma_u[i_j]/lambda_u[i_j];\
                  vec_tmp1[Y] = vec_kb[Y]/2 * sigma_u[i_j]/lambda_u[i_j];\
                  vec_tmp1[Z] = vec_kb[Z]/2 * sigma_u[i_j]/lambda_u[i_j];\
                  SUB(vec_tmp2,vec_tmp2,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		 }\
		}\
                for(i_b=1;i_b<=n_b;i_b++) {\
                   left_right_b[0] = i_b * left_right[0];\
                   left_right_b[1] = i_b * left_right[1];\
                   left_right_b[2] = i_b * left_right[2];\
                /* find the autoinduction of the left part  */\
                for(i_j=1;i_j<=NP-1;i_j++) {\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j]-left_right_b[0];\
                  vec_tmp1[Y] = uy_tmp[i_j]-left_right_b[1];\
                  vec_tmp1[Z] = uz_tmp[i_j]-left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
                /* find the autoinduction of the right  part  */\
                for(i_j=2;i_j<=NP;i_j++) {\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j]+left_right_b[0];\
                  vec_tmp1[Y] = uy_tmp[i_j]+left_right_b[1];\
                  vec_tmp1[Z] = uz_tmp[i_j]+left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
                }\
             /* free pointers  */\
             free(ux_tmp);\
             free(uy_tmp);\
             free(uz_tmp);\
             free(ux_s);\
             free(uy_s);\
             free(uz_s);\
             free(sigma_u);\
             free(lambda_u);
#endif
#endif
#endif /* End of if CALL_AND_TING */

#if DE_SINGU
#if CLOSED
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION  \
                      vec_auto[X] = 0.;\
	              vec_auto[Y] = 0.;\
	              vec_auto[Z] = 0.;\
                              \
                  cut_off_length_tmp = cut_off_length;\
                  cut_off_length_tmp2 = cut_off_length_tmp * cut_off_length_tmp;\
                  /*printf("%g ",cut_off_length_tmp);*/\
                  gamma_tmp = gamma_param[j]/4./M_PI;\
                 for(i_j=1;i_j<=NP-1;i_j++) {\
                  if (i_j!=i){\
                  vec_tmp[X] = Ux_s(i_j,j);\
                  vec_tmp[Y] = Uy_s(i_j,j);\
                  vec_tmp[Z] = Uz_s(i_j,j);\
                \
                  vec_tmp1[X] = Ux(i_j,j);\
                  vec_tmp1[Y] = Uy(i_j,j);\
                  vec_tmp1[Z] = Uz(i_j,j);\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1) * NORMV(vec_tmp1);\
                  one_o_norm3 = sqrt(one_o_norm3+cut_off_length_tmp2);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
                  }\
		 }
#endif
#endif

#if !CLOSED
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION  \
                      vec_auto[X] = 0.;\
	              vec_auto[Y] = 0.;\
	              vec_auto[Z] = 0.;\
                              \
             ux_tmp  = VECTOR(NP+1);\
             uy_tmp  = VECTOR(NP+1);\
             uz_tmp  = VECTOR(NP+1);\
             ux_s = VECTOR(NP+1); \
             uy_s = VECTOR(NP+1); \
             uz_s = VECTOR(NP+1); \
             sigma_u = VECTOR(NP+1);\
                                      \
           /* find the translation vector*/\
                \
                left_right[0] = Ux(NP,j)-Ux(1,j);\
                left_right[1] = Uy(NP,j)-Uy(1,j);\
                left_right[2] = Uz(NP,j)-Uz(1,j);\
                flag_overflow = 0;\
               \
           /* construct the u_tmp vector */\
             for(i_j=1;i_j<=(NP+1)/2;i_j++) {\
             int ll = i_j+i-1;\
               if (ll>NP) {\
                     ll = ll-(NP-1);\
                     flag_overflow = 1;\
               }\
	       if (flag_overflow==0){\
                  ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j);\
                  uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j);\
                  uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j);\
                /*printf("%d\t%d\t%d\n",flag_overflow,(NP+1)/2+i_j-1,ll);*/\
               }\
               else{\
                  ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j)+left_right[0];\
                  uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j)+left_right[1];\
                  uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j)+left_right[2];\
                /*printf("%d\t%d\t%d\n",flag_overflow,(NP+1)/2+i_j-1,ll);*/\
	       }\
             ux_s[(NP+1)/2+i_j-1] = Ux_s(ll,j);\
             uy_s[(NP+1)/2+i_j-1] = Uy_s(ll,j);\
             uz_s[(NP+1)/2+i_j-1] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2+i_j-1] =  SIGMA(ll,j);\
	     }\
                                 \
             flag_overflow = 0;\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
             int ll = i-i_j;\
               if (ll<1) {\
                     ll = (NP-1)+ll;\
                   flag_overflow = 1;\
               }\
	       if (flag_overflow==0){\
                  ux_tmp[(NP+1)/2-i_j] = Ux(ll,j);\
                  uy_tmp[(NP+1)/2-i_j] = Uy(ll,j);\
                  uz_tmp[(NP+1)/2-i_j] = Uz(ll,j);\
	      /*printf("%d\t%d\t%d\n",flag_overflow,(NP-1)/2-i_j+1,ll);*/\
               }\
               else{\
                  ux_tmp[(NP+1)/2-i_j] = Ux(ll,j)-left_right[0];\
                  uy_tmp[(NP+1)/2-i_j] = Uy(ll,j)-left_right[1];\
                  uz_tmp[(NP+1)/2-i_j] = Uz(ll,j)-left_right[2];\
                /*printf("%d\t%d\t%d\n",flag_overflow,(NP+1)/2-i_j,ll);*/\
	       }\
             ux_s[(NP+1)/2-i_j] = Ux_s(ll,j);\
             uy_s[(NP+1)/2-i_j] = Uy_s(ll,j);\
             uz_s[(NP+1)/2-i_j] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2-i_j] =  SIGMA(ll,j);\
	     }\
              /* then proceed to the auto_induction calculus  */\
                  cut_off_length_tmp = cut_off_length;\
                  cut_off_length_tmp2 = cut_off_length_tmp * cut_off_length_tmp;\
                  /*printf("%g ",cut_off_length_tmp);*/\
                  gamma_tmp = gamma_param[j]/4./M_PI;\
               /* find the translation vector*/\
                \
                left_right[0] = ux_tmp[NP]-ux_tmp[1];\
                left_right[1] = uy_tmp[NP]-uy_tmp[1];\
                left_right[2] = uz_tmp[NP]-uz_tmp[1];\
               \
                /* find the autoinduction of the central part  */\
                for(i_j=1;i_j<=NP;i_j++) {\
                 /* if (i_j!=(NP+1)/2){ */\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j];\
                  vec_tmp1[Y] = uy_tmp[i_j];\
                  vec_tmp1[Z] = uz_tmp[i_j];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = sqrt(one_o_norm3 *one_o_norm3 +cut_off_length_tmp2);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		 /*}*/\
		}\
                for(i_b=1;i_b<=n_b;i_b++) {\
                   left_right_b[0] = i_b * left_right[0];\
                   left_right_b[1] = i_b * left_right[1];\
                   left_right_b[2] = i_b * left_right[2];\
                /* find the autoinduction of the left part  */\
                for(i_j=1;i_j<=NP-1;i_j++) {\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j]-left_right_b[0];\
                  vec_tmp1[Y] = uy_tmp[i_j]-left_right_b[1];\
                  vec_tmp1[Z] = uz_tmp[i_j]-left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = sqrt(one_o_norm3 *one_o_norm3 +cut_off_length_tmp2);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
                /* find the autoinduction of the right  part  */\
                for(i_j=2;i_j<=NP;i_j++) {\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j]+left_right_b[0];\
                  vec_tmp1[Y] = uy_tmp[i_j]+left_right_b[1];\
                  vec_tmp1[Z] = uz_tmp[i_j]+left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = sqrt(one_o_norm3 *one_o_norm3 +cut_off_length_tmp2);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
                }\
             /* free pointers  */\
             free(ux_tmp);\
             free(uy_tmp);\
             free(uz_tmp);\
             free(ux_s);\
             free(uy_s);\
             free(uz_s);\
             free(sigma_u);
#endif
#endif
#endif /* End of if DE_SINGU */

#if M1_KNIO_KLEIN
#if CLOSED
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION  \
                      vec_auto[X] = 0.;\
	              vec_auto[Y] = 0.;\
	              vec_auto[Z] = 0.;\
                              \
                  delta_ttm_tmp = delta_ttm;\
                  ln_o_ln_coeff = log(sigma_1/delta_ttm_tmp)/log(sigma_2/sigma_1);\
                  /*printf("%g ",delta_ttm_tmp);*/\
                  gamma_tmp = gamma_param[j]/4./M_PI;\
                  for(i_j=1;i_j<=NP-1;i_j++) {\
                  if (i_j!=i){\
                  vec_tmp[X] = Ux_s(i_j,j);\
                  vec_tmp[Y] = Uy_s(i_j,j);\
                  vec_tmp[Z] = Uz_s(i_j,j);\
                \
                  vec_tmp1[X] = Ux(i_j,j);\
                  vec_tmp1[Y] = Uy(i_j,j);\
                  vec_tmp1[Z] = Uz(i_j,j);\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1) * NORMV(vec_tmp1);\
                  one_o_norm3 = sqrt(one_o_norm3);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  norm3 = one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                  kernel_tanh_coeff_1 = tanh(norm3/sigma3_1);\
                  kernel_tanh_coeff_2 = tanh(norm3/sigma3_2);\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                \
                  vec_auto1[X] = kernel_tanh_coeff_1 * vec_tmp2[X];\
	          vec_auto1[Y] = kernel_tanh_coeff_1 * vec_tmp2[Y];\
	          vec_auto1[Z] = kernel_tanh_coeff_1 * vec_tmp2[Z];\
                \
                  vec_auto2[X] = kernel_tanh_coeff_2 * vec_tmp2[X];\
	          vec_auto2[Y] = kernel_tanh_coeff_2 * vec_tmp2[Y];\
	          vec_auto2[Z] = kernel_tanh_coeff_2 * vec_tmp2[Z];\
                \
                  vec_tmp2[X] = vec_auto1[X] + (vec_auto1[X]-vec_auto2[X]) * ln_o_ln_coeff;\
	          vec_tmp2[Y] = vec_auto1[Y] + (vec_auto1[Y]-vec_auto2[Y]) * ln_o_ln_coeff;\
	          vec_tmp2[Z] = vec_auto1[Z] + (vec_auto1[Z]-vec_auto2[Z]) * ln_o_ln_coeff;\
                \
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
                  }\
		 }	
#endif
#endif

#if !CLOSED
#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION  \
             vec_auto[X] = 0.;\
	     vec_auto[Y] = 0.;\
	     vec_auto[Z] = 0.;\
                              \
             ux_tmp  = VECTOR(NP+1);\
             uy_tmp  = VECTOR(NP+1);\
             uz_tmp  = VECTOR(NP+1);\
             ux_s = VECTOR(NP+1); \
             uy_s = VECTOR(NP+1); \
             uz_s = VECTOR(NP+1); \
             sigma_u = VECTOR(NP+1);\
                                      \
           /* find the translation vector*/\
                \
                left_right[0] = Ux(NP,j)-Ux(1,j);\
                left_right[1] = Uy(NP,j)-Uy(1,j);\
                left_right[2] = Uz(NP,j)-Uz(1,j);\
                flag_overflow = 0;\
               \
           /* construct the u_tmp vector */\
             for(i_j=1;i_j<=(NP+1)/2;i_j++) {\
             int ll = i_j+i-1;\
               if (ll>NP) {\
                     ll = ll-(NP-1);\
                     flag_overflow = 1;\
               }\
	       if (flag_overflow==0){\
                  ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j);\
                  uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j);\
                  uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j);\
                /*printf("%d\t%d\t%d\n",flag_overflow,(NP+1)/2+i_j-1,ll);*/\
               }\
               else{\
                  ux_tmp[(NP+1)/2+i_j-1] = Ux(ll,j)+left_right[0];\
                  uy_tmp[(NP+1)/2+i_j-1] = Uy(ll,j)+left_right[1];\
                  uz_tmp[(NP+1)/2+i_j-1] = Uz(ll,j)+left_right[2];\
                /*printf("%d\t%d\t%d\n",flag_overflow,(NP+1)/2+i_j-1,ll);*/\
	       }\
             ux_s[(NP+1)/2+i_j-1] = Ux_s(ll,j);\
             uy_s[(NP+1)/2+i_j-1] = Uy_s(ll,j);\
             uz_s[(NP+1)/2+i_j-1] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2+i_j-1] =  SIGMA(ll,j);\
	     }\
                                 \
             flag_overflow = 0;\
             for(i_j=1;i_j<=(NP+1)/2-1;i_j++) {\
             int ll = i-i_j;\
               if (ll<1) {\
                     ll = (NP-1)+ll;\
                   flag_overflow = 1;\
               }\
	       if (flag_overflow==0){\
                  ux_tmp[(NP+1)/2-i_j] = Ux(ll,j);\
                  uy_tmp[(NP+1)/2-i_j] = Uy(ll,j);\
                  uz_tmp[(NP+1)/2-i_j] = Uz(ll,j);\
	      /*printf("%d\t%d\t%d\n",flag_overflow,(NP-1)/2-i_j+1,ll);*/\
               }\
               else{\
                  ux_tmp[(NP+1)/2-i_j] = Ux(ll,j)-left_right[0];\
                  uy_tmp[(NP+1)/2-i_j] = Uy(ll,j)-left_right[1];\
                  uz_tmp[(NP+1)/2-i_j] = Uz(ll,j)-left_right[2];\
                /*printf("%d\t%d\t%d\n",flag_overflow,(NP+1)/2-i_j,ll);*/\
	       }\
             ux_s[(NP+1)/2-i_j] = Ux_s(ll,j);\
             uy_s[(NP+1)/2-i_j] = Uy_s(ll,j);\
             uz_s[(NP+1)/2-i_j] = Uz_s(ll,j);\
             sigma_u[(NP+1)/2-i_j] =  SIGMA(ll,j);\
	     }\
              /* then proceed to the auto_induction calculus  */\
                              \
                  delta_ttm_tmp = delta_ttm;\
                  ln_o_ln_coeff = log(sigma_1/delta_ttm_tmp)/log(sigma_2/sigma_1);\
                 /*printf("%g ",delta_ttm_tmp);*/\
             gamma_tmp = gamma_param[j]/4./M_PI;\
               /* find the translation vector*/\
                \
                left_right[0] = ux_tmp[NP]-ux_tmp[1];\
                left_right[1] = uy_tmp[NP]-uy_tmp[1];\
                left_right[2] = uz_tmp[NP]-uz_tmp[1];\
               \
                /* find the autoinduction of the central part  */\
                for(i_j=1;i_j<=NP;i_j++) {\
                 if (i_j!=(NP+1)/2){\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j];\
                  vec_tmp1[Y] = uy_tmp[i_j];\
                  vec_tmp1[Z] = uz_tmp[i_j];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  norm3 = one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                  kernel_tanh_coeff_1 = tanh(norm3/sigma3_1);\
                  kernel_tanh_coeff_2 = tanh(norm3/sigma3_2);\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
               \
                  vec_auto1[X] = kernel_tanh_coeff_1 * vec_tmp2[X];\
	          vec_auto1[Y] = kernel_tanh_coeff_1 * vec_tmp2[Y];\
	          vec_auto1[Z] = kernel_tanh_coeff_1 * vec_tmp2[Z];\
                \
                  vec_auto2[X] = kernel_tanh_coeff_2 * vec_tmp2[X];\
	          vec_auto2[Y] = kernel_tanh_coeff_2 * vec_tmp2[Y];\
	          vec_auto2[Z] = kernel_tanh_coeff_2 * vec_tmp2[Z];\
                \
                  vec_tmp2[X] = vec_auto1[X] + (vec_auto1[X]-vec_auto2[X]) * ln_o_ln_coeff;\
	          vec_tmp2[Y] = vec_auto1[Y] + (vec_auto1[Y]-vec_auto2[Y]) * ln_o_ln_coeff;\
	          vec_tmp2[Z] = vec_auto1[Z] + (vec_auto1[Z]-vec_auto2[Z]) * ln_o_ln_coeff;\
                \
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		 }\
		}\
                for(i_b=1;i_b<=n_b;i_b++) {\
                   left_right_b[0] = i_b * left_right[0];\
                   left_right_b[1] = i_b * left_right[1];\
                   left_right_b[2] = i_b * left_right[2];\
                /* find the autoinduction of the left part  */\
                for(i_j=1;i_j<=NP-1;i_j++) {\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j]-left_right_b[0];\
                  vec_tmp1[Y] = uy_tmp[i_j]-left_right_b[1];\
                  vec_tmp1[Z] = uz_tmp[i_j]-left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
                /* find the autoinduction of the right  part  */\
                for(i_j=2;i_j<=NP;i_j++) {\
                  vec_tmp[X] = ux_s[i_j];\
                  vec_tmp[Y] = uy_s[i_j];\
                  vec_tmp[Z] = uz_s[i_j];\
                \
                  vec_tmp1[X] = ux_tmp[i_j]+left_right_b[0];\
                  vec_tmp1[Y] = uy_tmp[i_j]+left_right_b[1];\
                  vec_tmp1[Z] = uz_tmp[i_j]+left_right_b[2];\
                \
                  SUB(vec_tmp1,point_c,vec_tmp1);\
                  one_o_norm3 = NORMV(vec_tmp1);\
                  one_o_norm3 = one_o_norm3*one_o_norm3*one_o_norm3;\
                  one_o_norm3 = 1. / one_o_norm3;\
                \
                  vec_tmp1[X] = one_o_norm3 * vec_tmp1[X];\
                  vec_tmp1[Y] = one_o_norm3 * vec_tmp1[Y];\
                  vec_tmp1[Z] = one_o_norm3 * vec_tmp1[Z];\
                \
                  CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                  vec_auto[X] += gamma_tmp * ds * vec_tmp2[X];\
	          vec_auto[Y] += gamma_tmp * ds * vec_tmp2[Y];\
	          vec_auto[Z] += gamma_tmp * ds * vec_tmp2[Z];\
		}\
                }\
             /* free pointers  */\
             free(ux_tmp);\
             free(uy_tmp);\
             free(uz_tmp);\
             free(ux_s);\
             free(uy_s);\
             free(uz_s);\
             free(sigma_u);
#endif
#endif
#endif  /* End of if M1_KNIO_KLEIN */



    /* Add different contributions */

#if LOCAL_INDUCTION 
#ifndef ADD_TO_VEC
#   define ADD_TO_VEC  \
                        vec_add[X] += vec_local[X];\
                        vec_add[Y] += vec_local[Y];\
                        vec_add[Z] += vec_local[Z]
#endif 
#endif    

#if CALL_AND_TING 
#ifndef ADD_TO_VEC 
#   define ADD_TO_VEC  \
              vec_add[X] += vec_local[X];\
              vec_add[Y] += vec_local[Y];\
              vec_add[Z] += vec_local[Z];\
                                         \
              vec_add[X] += vec_auto[X];\
              vec_add[Y] += vec_auto[Y];\
              vec_add[Z] += vec_auto[Z];\
/*fprintf(fpp,"%g\t%g\t%g\n",vec_add[X],vec_add[Y],vec_add[Z]);*/
#endif  
#endif
             

#if DE_SINGU
#ifndef ADD_TO_VEC
#   define ADD_TO_VEC  \
              vec_add[X] += vec_auto[X];\
              vec_add[Y] += vec_auto[Y];\
              vec_add[Z] += vec_auto[Z]
#endif
#endif

#if M1_KNIO_KLEIN
#ifndef ADD_TO_VEC
#   define ADD_TO_VEC  \
              vec_add[X] += vec_auto[X];\
              vec_add[Y] += vec_auto[Y];\
              vec_add[Z] += vec_auto[Z];\
/*  fprintf(fpp,"%g\t%g\t%g\n",vec_add[X],vec_add[Y],vec_add[Z]);*/
#endif
#endif


#ifndef ADD_TO_VEC
#define ADD_TO_VEC
#endif

#ifndef FIND_GLOBAL_AUTO_INDUCTION
#define FIND_GLOBAL_AUTO_INDUCTION
#endif

    
/* -------------------------------------         
 * Derivative of the geometry
 * ------------------------------------- */

#ifndef FIND_DERIVATIVES_FINITE_DIFF
#define FIND_DERIVATIVES_FINITE_DIFF  \
             for(i=1;i<=NP;i++) {\
             point_m[X] = Ux(i-1,j);\
	     point_m[Y] = Uy(i-1,j);\
	     point_m[Z] = Uz(i-1,j);\
                                    \
             point_c[X] = Ux(i,j);\
	     point_c[Y] = Uy(i,j);\
	     point_c[Z] = Uz(i,j);\
                                  \
             point_p[X] = Ux(i+1,j);\
	     point_p[Y] = Uy(i+1,j);\
	     point_p[Z] = Uz(i+1,j);\
                                    \
	     /*  X_s derivative and find sigma */\
                                             \
             SUB(vec_tmp,point_p,point_m);\
             Ux_s(i,j) =  vec_tmp[X] * one_o_2ds;\
             Uy_s(i,j) =  vec_tmp[Y] * one_o_2ds;\
             Uz_s(i,j) =  vec_tmp[Z] * one_o_2ds;\
	     SIGMA(i,j) = sigma_tmp = NORMV(vec_tmp) * one_o_2ds;\
                                                  \
                                                  \
             /*  find X_ss ds^2 */\
             one_o_sigma3 = 1./(sigma_tmp * sigma_tmp * sigma_tmp);\
             ADD(vec_tmp1,point_p,point_m);\
             SUB(vec_tmp1,vec_tmp1,point_c);\
             SUB(vec_tmp1,vec_tmp1,point_c);\
                                            \
           /*  find K b */\
                                            \
             CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                                                   \
             vec_tmp2[X] = vec_tmp2[X] * one_o_2ds * one_o_ds2;\
	     vec_kb[X] = vec_tmp2[X] * one_o_sigma3;\
                                                    \
	     vec_tmp2[Y] = vec_tmp2[Y] * one_o_2ds * one_o_ds2;\
	     vec_kb[Y] = vec_tmp2[Y] * one_o_sigma3;\
                                                    \
             vec_tmp2[Z] = vec_tmp2[Z] * one_o_2ds * one_o_ds2;\
	     vec_kb[Z] = vec_tmp2[Z] * one_o_sigma3;\
                                                    \
             Ux_ss(i,j) = vec_kb[X];\
             Uy_ss(i,j) = vec_kb[Y];\
             Uz_ss(i,j) = vec_kb[Z];\
	     }
#endif

#ifndef FIND_DERIVATIVES_SPECTRAL
#define FIND_DERIVATIVES_SPECTRAL \
             /*  find X_s, sigma and K b by spectral method */\
                                    \
             /* fft parameters */   \
                                    \
            for(i=1;i<=NP-1;i++) {\
            U_spect(i-1,X) = Ux(i,j);\
	    U_spect(i-1,Y) = Uy(i,j);\
	    U_spect(i-1,Z) = Uz(i,j);\
	    }\
             \
            /* zero Nyquist mode and smoothing */\
            for(i=NP;i<=NP+1;i++) {\
            U_spect(i-1,X) = 0.;\
	    U_spect(i-1,Y) = 0.;\
	    U_spect(i-1,Z) = 0.;\
	    }\
                                \
    inc = 1;\
    jump = NP+1;\
    lot = 3;\
                                \
    isgn = -1;\
     fft991_(u_spect, work, trigs, ifax, &inc, &jump, &np_fft, &lot, &isgn);\
             \
                                                 \
            /* zero Nyquist mode and smoothing */\
            for(i=NP;i<=NP+1;i++) {\
            U_spect(i-1,X) = 0.;\
	    U_spect(i-1,Y) = 0.;\
	    U_spect(i-1,Z) = 0.;\
	    }\
                                \
          for(i=1;i<= np_fft_o_2;i++) {\
          k = 2*i-1;\
          sum_tmp = 0.;\
	  for(jj=0;jj<= 2;jj++){\
             sum_tmp = sum_tmp + U_spect(k-1,jj) * U_spect(k-1,jj);\
             sum_tmp = sum_tmp + U_spect(k,jj) * U_spect(k,jj);\
	  }\
          sum_tmp = sqrt(sum_tmp);\
	  if (sum_tmp<=1e-12){\
	      for(jj=0;jj<= 2;jj++){\
              U_spect(k-1,jj)=0.;\
              U_spect(k,jj)=0.;\
	      }\
	  }\
	  }\
                               \
                               \
         /* Derivative  */\
         /*  Fourier X_s derivative */\
         for(jj=0;jj<= 2;jj++){\
          for(i=1;i<= np_fft_o_2;i++){\
          k = 2*i-1;\
          U_s_spect(k-1,jj) = - xm[k-1]* U_spect(k,jj);\
          U_s_spect(k,jj) = xm[k-1]* U_spect(k-1,jj);\
	  }\
	 }\
                                 \
                                 \
 	                         \
        /*  Fourier X_ss derivative */\
         for(jj=0;jj<= 2;jj++){\
          for(i=1;i<= np_fft_o_2;i++) {\
          U_ss_spect(i-1,jj) = - xm_2[i-1]* U_spect(i-1,jj);\
	 }\
	  }\
           \
           \
           \
        /* zero Nyquist mode */\
        for(jj=0;jj<= 2;jj++){\
          U_s_spect(NP-1,jj) =0.;\
          U_s_spect(NP,jj) =0.;\
          U_ss_spect(NP-1,jj) =0.;\
          U_ss_spect(NP,jj) =0.;\
	  }\
                                  \
        /* Inverse Fourier Transform */\
        isgn = 1;\
          fft991_(&U_spect(0,X), work, trigs, ifax, &inc, &jump, &np_fft, &lot, &isgn);\
                               \
           for(i=1;i<=NP-1;i++) {\
            Ux(i,j) = U_spect(i-1,X);\
	    Uy(i,j) = U_spect(i-1,Y);\
	    Uz(i,j) = U_spect(i-1,Z);\
	    }\
            \
         /*  X_s derivative and sigma*/\
        fft991_(&U_s_spect(0,X), work, trigs, ifax, &inc, &jump, &np_fft, &lot, &isgn);\
                           \
        for(i=1;i<=NP-1;i++) {\
            Ux_s(i,j) = U_s_spect(i-1,X);\
	    Uy_s(i,j) = U_s_spect(i-1,Y);\
	    Uz_s(i,j) = U_s_spect(i-1,Z);\
	}\
            Ux_s(NP,j) = U_s_spect(0,X);\
	    Uy_s(NP,j) = U_s_spect(0,Y);\
	    Uz_s(NP,j) = U_s_spect(0,Z);\
                                         \
          /*  X_ss derivative */\
        fft991_(&U_ss_spect(0,X), work, trigs, ifax, &inc, &jump, &np_fft, &lot, &isgn);\
               \
            for(i=1;i<=NP-1;i++) {\
            Ux_ss(i,j) = U_ss_spect(i-1,X);\
	    Uy_ss(i,j) = U_ss_spect(i-1,Y);\
	    Uz_ss(i,j) = U_ss_spect(i-1,Z);\
	    }\
            Ux_ss(NP,j) = U_ss_spect(0,X);\
	    Uy_ss(NP,j) = U_ss_spect(0,Y);\
	    Uz_ss(NP,j) = U_ss_spect(0,Z);
#endif

#ifndef FIND_KB_SPECTRAL
#define FIND_KB_SPECTRAL \
             /*  find K b */\
        for(i=1;i<=NP-1;i++) {\
            vec_tmp[X] = Ux_s(i,j);\
	    vec_tmp[Y] = Uy_s(i,j);\
	    vec_tmp[Z] = Uz_s(i,j);\
            SIGMA(i,j) = sigma_tmp = NORMV(vec_tmp);\
            one_o_sigma3 = 1./(sigma_tmp * sigma_tmp * sigma_tmp);\
                                     \
            vec_tmp1[X] = Ux_ss(i,j);\
	    vec_tmp1[Y] = Uy_ss(i,j);\
	    vec_tmp1[Z] = Uz_ss(i,j);\
                                     \
            CROSS_PROD(vec_tmp2,vec_tmp,vec_tmp1);\
                                                   \
	     vec_kb[X] = vec_tmp2[X] * one_o_sigma3;\
                                                    \
	     vec_kb[Y] = vec_tmp2[Y] * one_o_sigma3;\
                                                    \
	     vec_kb[Z] = vec_tmp2[Z] * one_o_sigma3;\
                                                    \
             Ux_ss(i,j) = vec_kb[X];\
             Uy_ss(i,j) = vec_kb[Y];\
             Uz_ss(i,j) = vec_kb[Z];\
	}\
            Ux_ss(NP,j) = Ux_ss(1,j);\
	    Uy_ss(NP,j) = Uy_ss(1,j);\
	    Uz_ss(NP,j) = Uz_ss(1,j);\
            SIGMA(NP,j) = SIGMA(1,j);
#endif 
#endif /*  _EZSTEP3D_  */
