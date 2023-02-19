/* ------------------------------------------------------------------------- *
 * ezvortex.h -- header file for EZ-Vortex
 *
 * Copyright (C) 2000 Daniel Margerit
 *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 16:43:23 $
 * ------------------------------------------------------------------------- */

#ifndef _EZVORTEX_
#define _EZVORTEX_

/* -------------------------------------------------------------------------  */
/* These are the main parameters affecting the compilation of EZ-Filament.    */
/* Define them according to your needs.                                       */
/*                                                                            */
typedef double  Real;      /* precision of Real variables (float or double)   */
                           /*                                             */
                           /* Choose of the integral method used :        */
#define CLOSED          1  /* either opened or closed filament            */
#define LOCAL_INDUCTION 1  /* either local induction                      */
#define CALL_AND_TING   0  /*     or Callegari and Ting                   */
#define DE_SINGU        0  /*     or a De-singularization  method         */
#define M1_KNIO_KLEIN   0  /*     M1 method of Knio/Klein                 */
#define STRUCTURE       1  /* 1 with the structure  or 0 without          */
#define SIMIL_PART      1  /* 1 similar part of vortex core               */
#define NON_SIMIL_PART  0  /* 1 non-similar part of the vortex core       */
#define UNIFORM_CORE    0  /* 1 uniform vortex core                       */
                           /*                                             */
                           /* Choose of the derivative method used :      */
#define SPECTRAL        0  /* either finite-diff. 0 or spectral 1         */
                           /* Choose of the time step method used :       */
#define EXPLICIT   0       /* either    explicit                          */
#define NEWTON     0       /*     or    implicit  with Newton solver      */
#define ADAMS_BASHFORTH 1  /*     or         Adams-Bashforth             */
#define GAUSS_SEIDEL   0   /* either Jacobi or Gauss-Seidel(explicit methods)*/
                           /*                                             */
                           /* Choose of the running  mode       :         */
#define GRAPHICS   1       /* if 1 then run with interactive graphics     */
#define MOVE       0       /* if 1 then you have an automatic motion of   */
                           /*  the graphic windows in the z */
                           /*  direction:see MOVE in the bottom of ezvortex.c*/
#define MOVIE      0       /* if 1 then snapshot the screen to do a movie */
#define COMPUTE    1       /* 1 if compute and store mode or 0 if history mode*/
#define CONV_ANALYSE  0     /* 1 to test the initial velocity convergence with np and n_b*/
/* ------------------------------------------------------------------------- */

/* 
 * I always define min, max, and make sure M_PI is defined 
 * ------------------------------------------------------- */
#define min(a,b)      ((a)>(b) ? (b) : (a))   
#define max(a,b)      ((a)>(b) ? (a) : (b))   
#ifndef M_PI
#define M_PI	      3.14159265358979323846   
#endif

#define M_GAMMA	      0.577215   

/* --------------------------------------------------- 
 * Global variables used throughout the EZ-Vortex Code 
 * (I use lots of them)
 * ---------------------------------------------------*/

#define  VECTOR(n)  (Real *)malloc((unsigned)(n)*sizeof(Real))


extern Real epsilon, nu_bar_param,  /*  common model parameters */
            *gamma_param, *m_0_param, /* non common model parameters */
            *delta_0_bar_param,
            *c_v, *c_w, *delta_bar,      /* useful parameter combinations */
            *S, *S_0, *m_flux, beta,            
            *fields, *u,            /* arrays for points of filaments and properties*/
            *int_S_o_S0_tmp, *int_S_o_S0, /* to store int S0/Sdt */
            *u_s, *u_ss,*sigma, *u_0, 
            *v, *v_0, *v_m,         /* Adams-Bashforth */
            *d, *c,  *a,                /* Non-similar core parameters */
            ds, dt, error_stop,     /* numerical parameters */
            error_stop,             /* error to stop the Newton iteration */
            *xm, *xm_2, *trigs;            /* Fourier int components */
extern long int  ifax[13];
extern int  np, np_fft, np_fft_o_2,      /* # of  points per filament */
            nf,                           /* # of  filaments */
            field_size,                  /* array size for each field */
            nsteps,                      /* # time steps to take */
            istep,                       /* current time step */
            n_laguerre_mode,
            n_b,                        /* number of periodic boxes */ 


            write_filament,              /* write filament flag */
            simulating_resolution,       /* graphics parameter */


            rotating_resolution,         /* graphics parameter */
            verbose;                     /* verbosity level */
extern Real xmin,xmax,ymin,ymax,zmin,zmax;  /* physical space window */




/* -------------------------------------------------------------------------
 * All index ranges throughout the code are expressed in terms of NP and NF.
 *  The code is generally more efficient if these are known numbers
 * at compile time, but then one must recompile for each change.  If the
 * values are known (eg specified on the compile line) then do nothing
 * here, otherwise define NP etc to be the *variables* np, etc.
 * ------------------------------------------------------------------------- */

/* ---- Dimensions for the filament -----*/

#ifndef NP
  #define NP  np
#endif

#ifndef NF
  #define NF  nf
#endif

#ifndef N_LAGUERRE_MODE
  #define N_LAGUERRE_MODE  n_laguerre_mode
#endif




/* ------------------------------------------------------------------------- 
 * Memory for the filament points (u and v) and the spatial sums (sigma_u    
 * and sigma_v) is allocated in Allocate_memory() in ezvortex.c.  These are      
 * allocated as long (single dimensional) arrays.  Here macros are defined   
 * so that one can easily reference a value corresponding to a particular    
 * grid point, i.e. macros are defined to treat all arrays as                
 * multi-dimensional. If you don't like the way I do this, you should be     
 * able to change it easily by making modifications here and in              
 * AllocateMem().  Let me know if you find a significant improvement.        
 *                                                                           
 * INDEX(i,j,k) converts grid point (i,j,k) to the array index.              
 *                                                                           
 * Ux(i,j),Uy(i,j),Uz(i,j) --  x,y,z coordinates of point i on curve U(j)       * Vx(i,j),Vy(i,j),Vz(i,j) --  x,y,z coordinates of point i on curve V(j)       * ------------------------------------------------------------------------- */
  
#define J_INC       ((NP+2)*3)
#define I_INC        3


#define INDEX(i,j)  ((i)*I_INC + (j)*J_INC)
#define INDEX_SCAL(i,j)  ((i)  + (j)*(NP+2))

#define Ux(i,j)           u[INDEX(i,j)]
#define Uy(i,j)           u[INDEX(i,j)+1]
#define Uz(i,j)           u[INDEX(i,j)+2]


#define Ux_s(i,j)        u_s[INDEX(i,j)]
#define Uy_s(i,j)        u_s[INDEX(i,j)+1]
#define Uz_s(i,j)        u_s[INDEX(i,j)+2]

#define Ux_ss(i,j)        u_ss[INDEX(i,j)]
#define Uy_ss(i,j)        u_ss[INDEX(i,j)+1]
#define Uz_ss(i,j)        u_ss[INDEX(i,j)+2]


#define SIGMA(i,j)        sigma[INDEX_SCAL(i,j)]

#define Ux_0(i,j)          u_0[INDEX(i,j)]
#define Uy_0(i,j)          u_0[INDEX(i,j)+1]
#define Uz_0(i,j)          u_0[INDEX(i,j)+2]

#define Vx(i,j)          v[INDEX(i,j)]
#define Vy(i,j)          v[INDEX(i,j)+1]
#define Vz(i,j)          v[INDEX(i,j)+2]

#define Vx_0(i,j)          v_0[INDEX(i,j)]
#define Vy_0(i,j)          v_0[INDEX(i,j)+1]
#define Vz_0(i,j)          v_0[INDEX(i,j)+2]

#define Vx_m(i,j)          v_m[INDEX(i,j)]
#define Vy_m(i,j)          v_m[INDEX(i,j)+1]
#define Vz_m(i,j)          v_m[INDEX(i,j)+2]


#define K_INC_SPECT      ((NP+1))
#define I_INC_SPECT      1
#define INDEX_SPECT(i,k)  ((i)*I_INC_SPECT +(k)*K_INC_SPECT)
#define U_spect(i,k)      u_spect[INDEX_SPECT(i,k)]
#define U_s_spect(i,k)    u_s_spect[INDEX_SPECT(i,k)]
#define U_ss_spect(i,k)   u_ss_spect[INDEX_SPECT(i,k)]



/* ------------------------------------------------------------------------- 
 *                          Core   coefficients   macros     
 * ------------------------------------------------------------------------- */
#define INDEX_CORE(in,j)  ((in)  + (j)*(N_LAGUERRE_MODE))
#define INDEX_A(in,im)  ((in)  + (im)*(N_LAGUERRE_MODE))
#define DD(in,j)         d[INDEX_CORE(in,j)] 
#define CC(in,j)         c[INDEX_CORE(in,j)] 
#define AA(in,im)        a[INDEX_A(in,im)]

/* ------------------------------------------------------------------------- 
 *                          Vector macros         
 * ------------------------------------------------------------------------- */


#define X 0                           /* Indices for vector components */
#define Y 1
#define Z 2
/* Vector operations used exclusively by the functions 
 *     */



#define SUB(ans,v1,v2) \
  (ans)[X] = (v1)[X] - (v2)[X]; \
  (ans)[Y] = (v1)[Y] - (v2)[Y]; \
  (ans)[Z] = (v1)[Z] - (v2)[Z]

#define ADD(ans,v1,v2) \
  (ans)[X] = (v1)[X] + (v2)[X]; \
  (ans)[Y] = (v1)[Y] + (v2)[Y]; \
  (ans)[Z] = (v1)[Z] + (v2)[Z]

#define DOT(v1, v2) ((v1)[X]*(v2)[X] + (v1)[Y]*(v2)[Y] + (v1)[Z]*(v2)[Z])

#define NORMV(v) sqrt((v)[X]*(v)[X] + (v)[Y]*(v)[Y] + (v)[Z]*(v)[Z])


#define CROSS_PROD(ans,v1,v2) \
  (ans)[X] = (v1)[Y] * (v2)[Z] - (v1)[Z] * (v2)[Y]; \
  (ans)[Y] = (v1)[Z] * (v2)[X] - (v1)[X] * (v2)[Z]; \
  (ans)[Z] = (v1)[X] * (v2)[Y] - (v1)[Y] * (v2)[X]





/* ------------------------------------------- 
 * Prototypes for public functions defined in: 
 * ------------------------------------------- */

/* ezvortex.c 
 * ---------- */
void Write_filament_data (float x0, float y0, float z0, 
			  float x1, float y1, float z1); 

 

/* ezstep3d.c 
 * ---------- */
void  Step      (FILE *fpp);
void  Step_ini  (void);




/* ezgraph3d.c 
 * ----------- */
void  Draw         (void);
void  Draw_ini     (int initial_field);
int   Event_check  (void);
void  QuitX        (void);
void  Save_image   (void);
void  Mover_continuous (Real m_x, Real m_y, Real m_z);

/* ezopengl.c 
 * ------------ */
void  EZplot3d  (unsigned int resolution, int field, Real *plot_length, Real *scale, Real *offset);
void  EZplot3d__ini    (void);

#endif /*  _EZVORTEX_  */

