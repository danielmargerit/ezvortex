/* ------------------------------------------------------------------------- */
/* ezopengl.h -- Macros for EZOPENGL
 *
 * Copyright (C) 2000 Daniel Margerit
 *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 16:43:23 $
 * ------------------------------------------------------------------------- */

#ifndef _EZOPENGL_
#define _EZOPENGL_




#define FILAMENT_R   0.0    /* R, */
#define FILAMENT_G   0.0    /* G, */
#define FILAMENT_B   1.0    /* B values for filament color: 0 <= R,G,B <= 1 */
#define FILAMENT_WT  3.0    /* filament line width */

#define BBOX_R       1.0    /* R, */
#define BBOX_G       0.0    /* G, */
#define BBOX_B       0.0    /* B values for bounding box: 0 <= R,G,B <= 1 */
#define BBOX_WT      2.0    /* bounding box line width */

/* ========================================================================= */
/* worm  lighting and texture parameters and other */

#define CORE_RADIUS 0.04
#define CIRCLE_RESOLUTION  12
#define AXIAL_RESOLUTION   NP /* this could be changed but <= NP*/




/* ------------------------------------------------------------------------- 
 * All floating point numbers in ezmarching are of type GLReal. Defined here
 * to be GLfloat. If the precision of the numerical data is double then
 * implicit conversion will take place when calculating the contour intersect
 * and finding the cube index. You may wish to experiment and try setting
 * this to double or to Real.  The functions glVertex* and glColor* functions
 * have many forms so I use macros to make changes easier. 
 *
 * Note that we often need to define points (and vectors) in 3 space
 * dimensions. These are used both in calculations and as data to pass to
 * OpenGL for rendering.  All such variables are defined as arrays of size
 * 3. In cases where we need arrays of coordinates, the index for each
 * dimension should be the last one. This means that the X, Y and Z
 * components are contiguous in memory and enables us to use glVertex3fv
 * instead of glVertex3f. 
 * ------------------------------------------------------------------------- */

typedef GLfloat       GLReal;

#define GLVERTEX3V(v)    glVertex3fv(v)
#define GLVERTEX3(x,y,z) glVertex3f((x),(y),(z))
#define GLCOLOR3(r,g,b)  glColor3f((r),(g),(b))

/* ------------------------------------------------------------------------- 
 *
 *            Probably you do not want to change anything below              
 *
 * ------------------------------------------------------------------------- */


typedef unsigned int  Axis;           /* Either X, Y or Z. */




#define EQUAL_ZERO(n)  ((n) == 0.0)   /* In this code, this should be
				       * ok because we always clamp
				       * values afterwards and we're
				       * not too * concerned about
				       * underflow and overflow. */

#define OFF       0                   /* Draw_modes */ 
#define LINES     1
#define TRIANGLES 2

#define X 0                           /* Indices for vector components */
#define Y 1
#define Z 2

#define A 0                           /* Params for implicit form of plane. */
#define B 1
#define C 2
#define D 3











#endif  /* _EZOPENGL_ */
