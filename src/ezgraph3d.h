/* ------------------------------------------------------------------------- */
/* ezgraph3d.h -- Macros for EZGRAPH3D
 *
 * Copyright (C) 2000 Daniel Margerit
 *
 * RCS Information
 * ---------------------------
 * $Revision:  $
 * $Date: 2000/06/29 16:43:23 $
 * ------------------------------------------------------------------------- */

#ifndef _EZGRAPH3D_
#define _EZGRAPH3D_

#define WINDOW_TITLE "EZ-Vortex"

#define WINX         0      /* Window location, in pixels, from screen left */
#define WINY         0      /* Window location, in pixels, from screen top. */
#define WM_CTRLS_POS 0      /* If WM_CTRL_POS is 0 then the window is placed
			     * at (WINX, WINY). If WM_CTRL_POS is 1 then WINX
			     * and WINY are ignored and the location is
			     * determined by the window manager. */
#define WINSIZE      500    /* Window is square of this size in pixels. */
#define PLOT_SIZE    1.0   /* This controls the size of the simulation
			     * volume in the view port: >1.0 for larger
			     * size, <1.0 for smaller size. */

#define PERSPECTIVE  1      /* 1 for perspective projection, 0 for orthographic.
			     * see myReshape(). */

#define DISTANCE     5.0    /* This is the distance that the screen is from
			     * from the eye (measured in units of the largest
			     * side of the simulations volume as it appears
			     * on the screen). This is only relevant for
			     * perspective projection. see myReshape() and
			     * polarView(). */

#define BACKGROUND   1.0    /* Background color (R=G=B=BACKGROUND, so 0.0 gives
			       BLACK, 1.0 gives WHITE) */

#define START_PAUSED 1      /* If 1 then window is opened in paused mode
			     * showing initial condition. */

#define INITIAL_THETA 130.0 /* Initial view (any value in [-180,180)). */
#define INITIAL_PHI   70.0  /* Initial view (any value in [-180,180)). */
#define ROTATE_SCALE   0.5  /* Number of degrees of rotation corresponding to
			     * one pixel of mouse motion. */

/* --------------------------------------------- 
 * You probably should not change anything below 
 * --------------------------------------------- */

#define TRUE             1   
#define FALSE            0

#define CURVE_DISPLAY    0 /* These are used to determine the display mode: */
#define WORM_DISPLAY     1 
#define NO_FIELD         -1 /* 0 for curve, 1 for worm , -1 for no display */

#define ISO_LIST         1

#define MODE_SIMULATING  1   
#define MODE_VIEWING     2   
#define MODE_ROTATING    3  

/*  Macros for converting from physical coordinates (x,y,z) to 
 *  pixel location (PX,PY) within one of the square plots.
 */
#define PX(x) ((Real)(plot_length[0]*scale[0]*((x)-offset[0])))      
#define PY(y) ((Real)(plot_length[1]*scale[1]*((y)-offset[1])))
#define PZ(z) ((Real)(plot_length[2]*scale[2]*((z)-offset[2])))


#endif /*  _EZGRAPH3D_  */
