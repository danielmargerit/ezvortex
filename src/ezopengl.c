/* ------------------------------------------------------------------------- */
/* ezopengl.c -- Ez routines to use opengl
 * 
 *
 * Copyright (C) 2000 Daniel Margerit
 *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 16:43:05 $
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glx.h>



#include "ezvortex.h"    
#include "ezgraph3d.h"
#include "ezopengl.h"

/* -------------------------------------------------------------------------
 *
 * This file contains all graphics functions that   
 *   render the filaments.
 *
 * The important things to know about this file are:
 * 
 * (1) All functions containing OpenGL calls start with Draw_.  These
 * functions and EZplot3d() are probably the only ones you would
 * need to change to affect the look of the graphics (contour colors etc.)
 *
 *
 * (2) There are three states the draw_mode can be in:
 *   OFF:       glEnd() has been called and nothing being passed to OpenGL.
 *   LINES:     glBegin(GL_LINES) has been called and points on line 
 *              segments are being passed.
 *   TRIANGLES: glBegin(TRIANGLES) has been called and triangle vertices 
 *              are being passed.
 *
 * 
 *
 * ------------------------------------------------------------------------- */


/* 
 * Global variables for this file only 
 * ----------------------------------- */



static int             draw_mode=OFF;
                            /* draw_mode can be in one of three states. start
			     * with it off */



/* ========================================================================= */
/* worm  parameters  */


/* ========================================================================= */




/* 
 * Private functions 
 * ----------------- */



/* OpenGL specific drawing functions */


static void       Set_draw_mode_lines     (GLReal lwidth);
static void       Set_draw_mode_off       (void);
static void       Draw_bounding_box       (Real *plot_length);
static void       Draw_filament           (GLReal point[2][3]);
static void       Draw_curve              (Real *plot_length, Real *scale, Real *offset, int nbr_max, int ind_filament);
static void       Draw_worm               (Real *plot_length, Real *scale, Real *offset, int nbr_max, int ntet, Real *core_radius, int flag, int ind_filament);





 /* This is the essential function. Here you can call the routine that you 
   *  want  to plot what you want; for instance, here it is Draw_bounding_box
   * that is called to draw a box. 
   * (This is the Marching_cube function of EZscroll) */
/* ========================================================================= */

void  EZplot3d  (unsigned int resolution, int field, Real *plot_length, Real *scale, Real *offset)
{

    Real core_radius =  CORE_RADIUS;
    int i;

   Draw_bounding_box(plot_length);

   if (field == NO_FIELD ) return;
   
 for(i=1;i<=NF;i++) {
    
   if (field == CURVE_DISPLAY ){
  
     Draw_curve(plot_length, scale, offset, AXIAL_RESOLUTION, i-1);
   }
   if (field == WORM_DISPLAY ){
        Draw_worm(plot_length, scale, offset, AXIAL_RESOLUTION, CIRCLE_RESOLUTION, &core_radius, 0, i-1);
     
   }

 }

    /* Finish up */
    Set_draw_mode_off();
  
}


/* ========================================================================= */
void EZplot3d_ini (void)
{
  
}


/* ========================================================================= */

static void Draw_filament (GLReal point[2][3])
{
  Set_draw_mode_lines(FILAMENT_WT);
  GLCOLOR3(FILAMENT_R, FILAMENT_G, FILAMENT_B);
  GLVERTEX3V(point[0]);
  GLVERTEX3V(point[1]);
}

/* ========================================================================= */

static void Set_draw_mode_lines (GLReal lwidth)
{
  if(draw_mode==LINES) return;

  if(draw_mode!=OFF) glEnd();
  glLineWidth(lwidth);
  glBegin(GL_LINES);
  draw_mode=LINES;
}
/* ========================================================================= */

static void Set_draw_mode_off (void)
{
  glEnd();
  draw_mode=OFF;
}
/* ========================================================================= */

void Draw_bounding_box (Real *plot_length)
{
  Set_draw_mode_off();
  Set_draw_mode_lines(BBOX_WT);

  GLCOLOR3(BBOX_R, BBOX_G, BBOX_B);

  GLVERTEX3(0.,0.,0.);  
  GLVERTEX3(plot_length[0],0.,0.);
  GLVERTEX3(plot_length[0],0.,0.);  
  GLVERTEX3(plot_length[0],plot_length[1],0.);
  GLVERTEX3(plot_length[0],plot_length[1],0.);  
  GLVERTEX3(0.,plot_length[1],0.);
  GLVERTEX3(0.,plot_length[1],0.);  
  GLVERTEX3(0.,0.,0.);

  GLVERTEX3(0.,0.,0.);  
  GLVERTEX3(0.,0.,plot_length[2]);
  GLVERTEX3(plot_length[0],0.,0.);  
  GLVERTEX3(plot_length[0],0.,plot_length[2]);
  GLVERTEX3(plot_length[0],plot_length[1],0.);  
  GLVERTEX3(plot_length[0],plot_length[1],plot_length[2]);
  GLVERTEX3(0.,plot_length[1],0.);  
  GLVERTEX3(0.,plot_length[1],plot_length[2]);

  GLVERTEX3(0.,0.,plot_length[2]);  
  GLVERTEX3(plot_length[0],0.,plot_length[2]);
  GLVERTEX3(plot_length[0],0.,plot_length[2]);  
  GLVERTEX3(plot_length[0],plot_length[1],plot_length[2]);
  GLVERTEX3(plot_length[0],plot_length[1],plot_length[2]);  
  GLVERTEX3(0.,plot_length[1],plot_length[2]);
  GLVERTEX3(0.,plot_length[1],plot_length[2]);  
  GLVERTEX3(0.,0.,plot_length[2]);

  Set_draw_mode_off();
}


/* ========================================================================= */
static void Draw_curve(Real *plot_length, Real *scale, Real *offset, int nbr_max, int ind_filament)
{
  static GLReal        point[2][3];       /* The line segment to be drawn. */  
  int i, i_tmp0, i_tmp1;
  Real rp;
       
  
    rp = (Real)(NP)/(Real)(nbr_max);
    for(i=1;i<=nbr_max-1;i++) {
          i_tmp0 = (int)(rp*(i));
          i_tmp1 = (int)(rp*(i+1));	      
	  point[0][0] = PX(Ux(i_tmp0,ind_filament));
	  point[0][1] = PY(Uy(i_tmp0,ind_filament));
	  point[0][2] = PZ(Uz(i_tmp0,ind_filament));
	  point[1][0] = PX(Ux(i_tmp1,ind_filament));
	  point[1][1] = PY(Uy(i_tmp1,ind_filament));
	  point[1][2] = PZ(Uz(i_tmp1,ind_filament));
	  Draw_filament(point);
          Set_draw_mode_off(); /* important to be set to have lighting in torus */
         
    } 

 if (NP!=nbr_max){
          point[0][0] = PX(Ux(1,ind_filament));
	  point[0][1] = PY(Uy(1,ind_filament));
	  point[0][2] = PZ(Uz(1,ind_filament));
	  point[1][0] = PX(Ux((int)(rp*(1)),ind_filament));
	  point[1][1] = PY(Uy((int)(rp*(1)),ind_filament));
	  point[1][2] = PZ(Uz((int)(rp*(1)),ind_filament));
	  Draw_filament(point);
          Set_draw_mode_off();

          point[1][0] = PX(Ux((int)(rp*(nbr_max-1)),ind_filament));
	  point[1][1] = PY(Uy((int)(rp*(nbr_max-1)),ind_filament));
	  point[1][2] = PZ(Uz((int)(rp*(nbr_max-1)),ind_filament));
          point[0][0] = PX(Ux(NP,ind_filament));
	  point[0][1] = PY(Uy(NP,ind_filament));
	  point[0][2] = PZ(Uz(NP,ind_filament));
	 
	  Draw_filament(point);
          Set_draw_mode_off();
 
   }
     
     
}  
/* ========================================================================= */



/*  Worm plot 3d curved cylinder (a worm) along the curve *u (list of NP  */
/*  stored points). The transversal radius  *core_radius may be uniform   */
/* if flag=0 or with axial variations if flag!=0. Here nbr_max is the number */
/* of axial points used to have the plot and ntet is the resolution, number */
/* of points to draw circles (12 is a good value) */

/* ========================================================================= */
static void Draw_worm(Real *plot_length, Real *scale, Real *offset, int nbr_max, int ntet, Real *core_radius, int flag, int ind_filament)
{

  static GLReal        point[2][3];       /* The line segment to be drawn. */  
  int i, i_tmp0, j;
  Real rp, *core_rad_tmp, *xc_tmp, *yc_tmp, *zc_tmp, *x_tmp, *y_tmp, *z_tmp;
  Real *nx_tmp, *ny_tmp, *nz_tmp;
  Real v_tmp[3], V_tmp, *x1, *x2, *y1, *y2, *z1, *z2, *dx_tmp, *dy_tmp, *dz_tmp;
  Real salpha, calpha, sbeta, cbeta, distance, valk;
  int *k_permut, kk;

  Real Tet_Step = 2.0 * M_PI / ntet;


  core_rad_tmp = (Real *) malloc((unsigned)(nbr_max)*sizeof(Real));
  xc_tmp       = (Real *) malloc((unsigned)(nbr_max)*sizeof(Real));
  yc_tmp       = (Real *) malloc((unsigned)(nbr_max)*sizeof(Real));
  zc_tmp       = (Real *) malloc((unsigned)(nbr_max)*sizeof(Real));

  x_tmp = (Real *) malloc((unsigned)(nbr_max)*(ntet+1)*sizeof(Real));
  y_tmp = (Real *) malloc((unsigned)(nbr_max)*(ntet+1)*sizeof(Real));
  z_tmp = (Real *) malloc((unsigned)(nbr_max)*(ntet+1)*sizeof(Real));


    dx_tmp = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
    dy_tmp = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
    dz_tmp = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
        x1 = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
        y1 = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
        z1 = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
        x2 = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
        y2 = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
        z2 = (Real *) malloc((unsigned)(ntet+1)*sizeof(Real));
  k_permut = (int *)  malloc((unsigned)(ntet+1)*sizeof(int));


#define I_INC_TMP        (ntet+1)     /* it is important to have the () here !!*/

#define INDEX_TMP(i,j)  (((i)-1)*I_INC_TMP+ (j))


#define XX(i,j)           x_tmp[INDEX_TMP((i),(j))]
#define YY(i,j)           y_tmp[INDEX_TMP((i),(j))]
#define ZZ(i,j)           z_tmp[INDEX_TMP((i),(j))]
#define NORM(a,b,c)       sqrt((a)*(a)+(b)*(b)+(c)*(c))
#define PRODSQR(a,b,c)    ((a)*(a)+(b)*(b)+(c)*(c))
  
#define NX(i,j)           nx_tmp[INDEX_TMP((i),(j))]
#define NY(i,j)           ny_tmp[INDEX_TMP((i),(j))]
#define NZ(i,j)           nz_tmp[INDEX_TMP((i),(j))]


  if (0==flag) { 
    for(i=1;i<=nbr_max;i++) {
      core_rad_tmp[i-1] = core_radius[0];}
  }
  else {
     core_rad_tmp = core_radius;
  }
    rp = (Real)(NP)/(Real)(nbr_max);


   for(i=1;i<=nbr_max;i++) {
           i_tmp0 = (int)(rp*(i));
	     
	   xc_tmp[i-1] = Ux(i_tmp0,ind_filament);
	   yc_tmp[i-1] = Uy(i_tmp0,ind_filament);
	   zc_tmp[i-1] = Uz(i_tmp0,ind_filament); 
      
    } 

   if (NP!=nbr_max){
           xc_tmp[0] = Ux(1,ind_filament);
	   yc_tmp[0] = Uy(1,ind_filament);
	   zc_tmp[0] = Uz(1,ind_filament); 
           xc_tmp[nbr_max-1] = Ux(NP,ind_filament);
	   yc_tmp[nbr_max-1] = Uy(NP,ind_filament);
	   zc_tmp[nbr_max-1] = Uz(NP,ind_filament); 
   }
    
  


for (i = 1; i <= nbr_max-1; i++) {
         i_tmp0 = i-1;
         v_tmp[0] = xc_tmp[i_tmp0+1]-xc_tmp[i_tmp0];
         v_tmp[1] = yc_tmp[i_tmp0+1]-yc_tmp[i_tmp0];
         v_tmp[2] = zc_tmp[i_tmp0+1]-zc_tmp[i_tmp0];
         V_tmp = NORM(v_tmp[0], v_tmp[1], v_tmp[2]);
         salpha= v_tmp[2]/V_tmp;
         calpha = sqrt(v_tmp[0]*v_tmp[0]+v_tmp[1]*v_tmp[1])/V_tmp;
	 sbeta = 0.;
         cbeta = 1.;
         if (0!=calpha) {
            sbeta = v_tmp[1]/V_tmp/calpha;
            cbeta = v_tmp[0]/V_tmp/calpha;}


       for (j = 0; j <= ntet-1; j++) {
	 Real tet = j * Tet_Step;
         Real ct = cos(tet);
         Real st = sin(tet);
         /* Here I could have used storage in a vector [ntet] for ct and st instead */
 
 

         XX(i,j) = ct*salpha*cbeta+st*sbeta;
         YY(i,j) = ct*salpha*sbeta-st*cbeta;
         ZZ(i,j) = -ct*calpha;  
     
       }
       XX(i,ntet) = XX(i,0);   /* I could have gone up to ntet in the previous loop */
       YY(i,ntet) = YY(i,0);   /* and do not applied this condition */
       ZZ(i,ntet) = ZZ(i,0);
       

}
 
 
   


/* move and scale circles properly */
   for (j = 0; j <= ntet; j++) {
         XX(nbr_max,j) = XX(nbr_max-1,j)*core_rad_tmp[nbr_max-1]+xc_tmp[nbr_max-1];
         YY(nbr_max,j) = YY(nbr_max-1,j)*core_rad_tmp[nbr_max-1]+yc_tmp[nbr_max-1];
         ZZ(nbr_max,j) = ZZ(nbr_max-1,j)*core_rad_tmp[nbr_max-1]+zc_tmp[nbr_max-1];
 
       }

 



     /* rotate to match with previous */
       for (j = 0; j <= ntet; j++) {
         x1[j] = XX(1,j); y1[j] = YY(1,j); z1[j] = ZZ(1,j);
       }

for (i = 2; i <= nbr_max-1; i++) {
         
       for (j = 0; j <= ntet; j++) {
        
         dx_tmp[j] = x1[0] - XX(i,j); 
	 dy_tmp[j] = y1[0] - YY(i,j);
	 dz_tmp[j] = z1[0] - ZZ(i,j);
	 dx_tmp[j] = PRODSQR(dx_tmp[j], dy_tmp[j], dz_tmp[j]);
       }
      
     

       kk =0; valk = dx_tmp[0] ;
       for (j = 1; j < ntet; j++) {
        if (dx_tmp[j]<valk){
           valk = dx_tmp[j];
           kk = j;
	}
       }
       
       for (j = 0; j <= ntet-1; j++) {
	 if ((kk + j)<ntet){
	 k_permut[j] = kk + j;
	 }
         else{
	 k_permut[j] = kk + j - (ntet);  /* this is not easy */
	 }

       }
         k_permut[ntet] = k_permut[0]; /* this is not easy */

     for (j = 0; j <= ntet; j++) {
        x2[j] = XX(i,k_permut[j]); 
        y2[j] = YY(i,k_permut[j]); 
        z2[j] = ZZ(i,k_permut[j]);       
     }
      

     for (j = 0; j <= ntet; j++) {	 

	 XX(i,j) = 0.5*(x1[j]+x2[j])*core_rad_tmp[i-1] + xc_tmp[i-1]; 
	 YY(i,j) = 0.5*(y1[j]+y2[j])*core_rad_tmp[i-1] + yc_tmp[i-1];
	 ZZ(i,j) = 0.5*(z1[j]+z2[j])*core_rad_tmp[i-1] + zc_tmp[i-1];
 
        x1[j] = x2[j];
        y1[j] = y2[j];
        z1[j] = z2[j];
        
   }
}

  



   for (j = 0; j <= ntet; j++) {
         dx_tmp[j] =  XX(nbr_max,k_permut[j]);
         dy_tmp[j] =  YY(nbr_max,k_permut[j]);
         dz_tmp[j] =  ZZ(nbr_max,k_permut[j]);
   }

    for (j = 0; j <= ntet; j++) {
         XX(nbr_max,j) = dx_tmp[j] ; 
	 YY(nbr_max,j) = dy_tmp[j];
	 ZZ(nbr_max,j) = dz_tmp[j];
   }

   for (j = 0; j <= ntet; j++) {
     
         XX(1,j) = XX(1,j)*core_rad_tmp[0] + xc_tmp[0]; 
	 YY(1,j) = YY(1,j)*core_rad_tmp[0] + yc_tmp[0];
	 ZZ(1,j) = ZZ(1,j)*core_rad_tmp[0] + zc_tmp[0];
   }
   /* check for periodicity */
     distance = PRODSQR(xc_tmp[0] - xc_tmp[nbr_max-1], yc_tmp[0] - yc_tmp[nbr_max-1], zc_tmp[0] - zc_tmp[nbr_max-1]);
       
     
    

     if (distance<1e-12){ 

    /* rotate to match with previous */
       for (j = 0; j <= ntet; j++) {
         x1[j] = XX(nbr_max,j); y1[j] = YY(nbr_max,j); z1[j] = ZZ(nbr_max,j);
       }
         
       for (j = 0; j <= ntet; j++) {
        
         dx_tmp[j] = x1[0] - XX(1,j); 
	 dy_tmp[j] = y1[0] - YY(1,j);
	 dz_tmp[j] = z1[0] - ZZ(1,j);
	 dx_tmp[j] = PRODSQR(dx_tmp[j], dy_tmp[j], dz_tmp[j]);
       }
      
     

       kk =0; valk = dx_tmp[0] ;
       for (j = 1; j < ntet; j++) {
        if (dx_tmp[j]<valk){
           valk = dx_tmp[j];
           kk = j;
	}
       }
       
       for (j = 0; j <= ntet-1; j++) {
	 if ((kk + j)<ntet){
	 k_permut[j] = kk + j;
	 }
         else{
	 k_permut[j] = kk + j - (ntet);  /* this is not easy */
	 }

       }
         k_permut[ntet] = k_permut[0]; /* this is not easy */



        for (j = 0; j <= ntet; j++) {
         XX(nbr_max,j) = 0.5*(XX(1,k_permut[j]) + XX(nbr_max,j)); 
	 YY(nbr_max,j) = 0.5*(YY(1,k_permut[j]) + YY(nbr_max,j));
	 ZZ(nbr_max,j) = 0.5*(ZZ(1,k_permut[j]) + ZZ(nbr_max,j));

	 XX(1,k_permut[j]) = XX(nbr_max,j); 
	 YY(1,k_permut[j]) = YY(nbr_max,j);
	 ZZ(1,k_permut[j]) = ZZ(nbr_max,j);
      
         }
     }


     /* now draw  the worm */






#if 1   /* plot the worm */

   glClearColor(0.9, 0.9, 0.9, 1.0);
     /*  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);*/




    glEnable(GL_DEPTH_TEST); 
   
 
    glColor3f(0.2, 0.2, 0.2);
  

  glPushMatrix();






for (i =1 ; i <=nbr_max-1  ; i++) {
          
#if 1   /*plot the color triangles */    
     glBegin(GL_TRIANGLE_STRIP);
    

      for (j = 0; j <= ntet; j++) {  
 
      
      GLReal x_0 = (GLReal)PX(XX(i,j));
      GLReal y_0 = (GLReal)PY(YY(i,j));
      GLReal z_0 = (GLReal)PZ(ZZ(i,j));
      GLReal x_1 = (GLReal)PX(XX(i+1,j));
      GLReal y_1 = (GLReal)PY(YY(i+1,j));
      GLReal z_1 = (GLReal)PZ(ZZ(i+1,j));
     
	 

 
         glTexCoord2f(i / (GLReal) nbr_max, j / (GLReal)ntet);
  
      glColor3f(x_0, y_0, z_0 ); 
  
  
    	 glVertex3f(x_0, y_0, z_0);
         glTexCoord2f((i + 1) / (GLReal)nbr_max , j / (GLReal)ntet);
 
          glColor3f(x_1, y_1, z_1); 
       
	 glVertex3f(x_1, y_1, z_1);
    }
    glEnd();



#endif    /* plot the color triangles */

#if 1   /* plot also the border  lines */

 glBegin(GL_LINES);
    
    for (j = 0; j <= ntet-1; j++) {  
 
      
      GLReal x_0 = (GLReal)PX(XX(i,j));
      GLReal y_0 = (GLReal)PY(YY(i,j));
      GLReal z_0 = (GLReal)PZ(ZZ(i,j));
      GLReal x_1 = (GLReal)PX(XX(i,j+1));
      GLReal y_1 = (GLReal)PY(YY(i,j+1));
      GLReal z_1 = (GLReal)PZ(ZZ(i,j+1));
     
	 

 
       
   
         glColor3f(0, 0, 0);  
    	 glVertex3f(x_0, y_0, z_0);      
	 glVertex3f(x_1, y_1, z_1);
    }






      for (j = 0; j <= ntet; j++) {  
 
      
      GLReal x_0 = (GLReal)PX(XX(i,j));
      GLReal y_0 = (GLReal)PY(YY(i,j));
      GLReal z_0 = (GLReal)PZ(ZZ(i,j));
      GLReal x_1 = (GLReal)PX(XX(i+1,j));
      GLReal y_1 = (GLReal)PY(YY(i+1,j));
      GLReal z_1 = (GLReal)PZ(ZZ(i+1,j));
     
	 

 
    
  
         glColor3f(0, 0, 0); 
    	 glVertex3f(x_0, y_0, z_0);       
	 glVertex3f(x_1, y_1, z_1);
    }
    glEnd();

#endif /* end of plot also the border lines */


}

  /* glDisable(GL_DEPTH_TEST);*/
   glDisable(GL_LIGHTING);
   glDisable(GL_TEXTURE_2D);
   glDisable(GL_BLEND);
   glDisable(GL_TEXTURE_GEN_S);
   glDisable(GL_TEXTURE_GEN_T);
  
         
   
     
#endif   


#if 0   /* plot also the central line */
        glDisable(GL_DEPTH_TEST);
    
    for(i=1;i<=nbr_max-1;i++) {
      GLReal x0 = PX(xc_tmp[i-1]);
      GLReal y0 = PY(yc_tmp[i-1]);
      GLReal z0 = PZ(zc_tmp[i-1]);
      GLReal x1 = PX(xc_tmp[i]);
      GLReal y1 = PY(yc_tmp[i]);
      GLReal z1 = PZ(zc_tmp[i]);

         	      
	  point[0][0] = x0;
	  point[0][1] = y0;
	  point[0][2] = z0;
	  point[1][0] = x1;
	  point[1][1] = y1;
	  point[1][2] = z1;
	  Draw_filament(point);
          Set_draw_mode_off(); /* important to be set to have lighting in torus */
         
    } 
       glEnable(GL_DEPTH_TEST);
     
     
#endif


    free(core_rad_tmp);
    free(xc_tmp); free(yc_tmp); free(zc_tmp);
    free(x_tmp); free(y_tmp); free(z_tmp);
    free(x1); free(y1); free(z1);
    free(x2); free(y2); free(z2);
    free(dx_tmp); free(dy_tmp); free(dz_tmp);
    free(k_permut);




}


/* ========================================================================= */









