/* fft.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <stdio.h>
#include <time.h>
#include <math.h>

/* Table of constant values */

static integer c__2 = 2;
static integer c__9 = 9;
static integer c__1 = 1;


/* PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE */
/*              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX */
/*              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE */
/*              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE */
/*              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER */
/*              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE */
/*              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS */
/*              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5. */
/*              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE */
/*              THAT IS MOSTLY WRITTEN IN CAL. */

/*              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES: */

/*            SUBROUTINE SET99 */
/*                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE */
/*                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES */
/*                (PROVIDED THAT N IS NOT CHANGED). */

/*            SUBROUTINES FFT99 AND FFT991 */
/*                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT */
/*                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE. */


/* ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH */

/*                   *FORTRAN,P=XLIB,SN=FFT99F */

/*              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY */
/*              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE */
/*              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE */
/*              ACCESSED USING */

/*                   FETCH P=CRAYLIB,SN=FFT99 */
/*                   FETCH P=CRAYLIB,SN=CAL99 */

/* USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1, */
/*              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF */
/*              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH */
/*              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS */
/*              OF LENGTH N IS */

/*                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)), */
/*                  +          WORK(M*(N+1)) */

/*                   CALL SET99 (TRIGS, IFAX, N) */
/*                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN) */

/*              SEE THE INDIVIDUAL WRITE-UPS FOR SET99, FFT99, AND */
/*              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE */
/*              ARGUMENTS. */

/* HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN */
/*              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED */
/*              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980. */

/* ----------------------------------------------------------------------- */

/* SUBROUTINE SET99 (TRIGS, IFAX, N) */

/* PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE */
/*              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT */
/*              ROUTINES (PROVIDED THAT N IS NOT CHANGED). */

/* ARGUMENT     IFAX(13),TRIGS(3*N/2+1) */
/* DIMENSIONS */

/* ARGUMENTS */

/* ON INPUT     TRIGS */
/*               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS */
/*               EVEN, OR 3*N/2+1 IF N/2 IS ODD. */

/*              IFAX */
/*               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED */
/*               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING */
/*               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION. */

/*              N */
/*               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR */
/*               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE */
/*               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE */
/*               DEFINITIONS OF THE TRANSFORMS). */

/* ON OUTPUT    IFAX */
/*               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE */
/*               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED */
/*               IN IFAX(2),IFAX(3),...  IF SET99 IS CALLED WITH N ODD, */
/*               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1) */
/*               IS SET TO -99. */

/*              TRIGS */
/*               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY */
/*               USED BY THE FFT ROUTINES. */

/* ----------------------------------------------------------------------- */

/* SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN) */
/*                       AND */
/* SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN) */

/* PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX */
/*              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE */
/*              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT */
/*              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE */
/*              GRIDPOINT VALUES (FFT99).  GIVEN A SET */
/*              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF */
/*              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE */
/*              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN */
/*              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS */
/*              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE */
/*              OPTIMIZED FOR USE ON THE CRAY-1. */

/* ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13) */
/* DIMENSIONS */

/* ARGUMENTS */

/* ON INPUT     A */
/*               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA */
/*               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY */
/*               THE RESULTS. */

/*              WORK */
/*               A WORK ARRAY OF DIMENSION M*(N+1) */

/*              TRIGS */
/*               AN ARRAY SET UP BY SET99, WHICH MUST BE CALLED FIRST. */

/*              IFAX */
/*               AN ARRAY SET UP BY SET99, WHICH MUST BE CALLED FIRST. */

/*              INC */
/*               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF */
/*               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR */
/*               CONSECUTIVELY STORED DATA). */

/*              JUMP */
/*               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF */
/*               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1, */
/*               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8 */
/*               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF */
/*               INC AND JUMP, SEE THE EXAMPLES BELOW. */

/*              N */
/*               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF */
/*               TRANSFORMS, BELOW). */

/*              M */
/*               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY. */

/*              ISIGN */
/*               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO */
/*                    GRIDPOINT VALUES. */
/*               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER */
/*                    COEFFICIENTS. */

/* ON OUTPUT    A */
/*               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED */
/*               EACH CONTAINING THE SEQUENCE: */

/*               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES) */

/*               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH */
/*               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES: */

/*               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0. */
/*               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0). */
/*                   (EXPLICIT CYCLIC CONTINUITY) */

/*               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY: */
/*                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N)) */
/*                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K) */
/*                 AND I=SQRT (-1) */

/*               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH */
/*               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS */
/*               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS */
/*               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS */
/*               A(K), B(K), 0 .LE. K .LE N/2. */

/*               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY: */
/*                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N)) */
/*                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1) */

/*               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1 */
/*               (OR VICE VERSA) RETURNS THE ORIGINAL DATA. */

/*               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL */
/*               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1, */
/*               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS. */

/* EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT */
/*               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF */
/*               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE */
/*               ARRANGED LIKE THIS: */

/* FIRST DATA   A(1)=    . . .                A(66)=             A(70) */
/* VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS) */

/* SECOND DATA  A(71)=   . . .                                  A(140) */
/* VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS) */

/*               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1, */
/*               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC */
/*               CONTINUITY). */

/*               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS: */

/*                FIRST         SECOND                          LAST */
/*                DATA          DATA                            DATA */
/*                VECTOR        VECTOR                          VECTOR */

/*                 A(1)=         A(2)=                           A(19)= */

/*                 X(63)         X(63)       . . .               X(63) */
/*        A(20)=   X(0)          X(0)        . . .               X(0) */
/*        A(39)=   X(1)          X(1)        . . .               X(1) */
/*                  .             .                               . */
/*                  .             .                               . */
/*                  .             .                               . */

/*               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING */
/*               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH */
/*               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT */
/*               DATA VECTOR. */

/* ----------------------------------------------------------------------- */

/*     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM */
/*     CORRESPONDING TO OLD SCALAR ROUTINE FFT9 */
/*     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM */
/*     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12 */
/*     (1970), 315-337) */

/*     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA */
/*     WORK IS AN AREA OF SIZE (N+1)*LOT */
/*     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES */
/*     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2 */
/*     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR' */
/*         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA) */
/*     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR */
/*     N IS THE LENGTH OF THE DATA VECTORS */
/*     LOT IS THE NUMBER OF DATA VECTORS */
/*     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT */
/*           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL */

/*     ORDERING OF COEFFICIENTS: */
/*         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2) */
/*         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED */

/*     ORDERING OF DATA: */
/*         X(N-1),X(0),X(1),X(2),...,X(N),X(0) */
/*         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED */

/*     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN */
/*     PARALLEL */

/*     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER */

/*     DEFINITION OF TRANSFORMS: */
/*     ------------------------- */

/*     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N)) */
/*         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K) */

/*     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N)) */
/*               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N)) */

integer nabs(integer x);



/* Subroutine */ int fft99a_(a, work, trigs, inc, jump, n, lot)
real *a, *work, *trigs;
integer *inc, *jump, *n, *lot;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static real c__;
    static integer k, l;
    static real s;
    static integer ia, ib, ja, jb, iabase, ibbase, nh, jabase, jbbase, nx, 
	    ink;


/*     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1 */
/*     (SPECTRAL TO GRIDPOINT TRANSFORM) */

    /* Parameter adjustments */
    --trigs;
    --work;
    --a;

    /* Function Body */
    nh = *n / 2;
    nx = *n + 1;
    ink = *inc + *inc;

/*     A(0) AND A(N/2) */
    ia = 1;
    ib = *n * *inc + 1;
    ja = 1;
    jb = 2;
/* DIR$ IVDEP */
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	work[ja] = a[ia] + a[ib];
	work[jb] = a[ia] - a[ib];
	ia += *jump;
	ib += *jump;
	ja += nx;
	jb += nx;
/* L10: */
    }

/*     REMAINING WAVENUMBERS */
    iabase = (*inc << 1) + 1;
    ibbase = (*n - 2) * *inc + 1;
    jabase = 3;
    jbbase = *n - 1;

    i__1 = nh;
    for (k = 3; k <= i__1; k += 2) {
	ia = iabase;
	ib = ibbase;
	ja = jabase;
	jb = jbbase;
	c__ = trigs[*n + k];
	s = trigs[*n + k + 1];
/* DIR$ IVDEP */
	i__2 = *lot;
	for (l = 1; l <= i__2; ++l) {
	    work[ja] = a[ia] + a[ib] - (s * (a[ia] - a[ib]) + c__ * (a[ia + *
		    inc] + a[ib + *inc]));
	    work[jb] = a[ia] + a[ib] + (s * (a[ia] - a[ib]) + c__ * (a[ia + *
		    inc] + a[ib + *inc]));
	    work[ja + 1] = c__ * (a[ia] - a[ib]) - s * (a[ia + *inc] + a[ib + 
		    *inc]) + (a[ia + *inc] - a[ib + *inc]);
	    work[jb + 1] = c__ * (a[ia] - a[ib]) - s * (a[ia + *inc] + a[ib + 
		    *inc]) - (a[ia + *inc] - a[ib + *inc]);
	    ia += *jump;
	    ib += *jump;
	    ja += nx;
	    jb += nx;
/* L20: */
	}
	iabase += ink;
	ibbase -= ink;
	jabase += 2;
	jbbase += -2;
/* L30: */
    }

    if (iabase != ibbase) {
	goto L50;
    }
/*     WAVENUMBER N/4 (IF IT EXISTS) */
    ia = iabase;
    ja = jabase;
/* DIR$ IVDEP */
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	work[ja] = a[ia] * (float)2.;
	work[ja + 1] = a[ia + *inc] * (float)-2.;
	ia += *jump;
	ja += nx;
/* L40: */
    }

L50:
    return 0;
} /* fft99a_ */

/* Subroutine */ int fft99b_(work, a, trigs, inc, jump, n, lot)
real *work, *a, *trigs;
integer *inc, *jump, *n, *lot;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static real c__;
    static integer k, l;
    static real s, scale;
    static integer ia, ib, ja, jb, iabase, ibbase, nh, jabase, jbbase, nx, 
	    ink;


/*     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1 */
/*     (GRIDPOINT TO SPECTRAL TRANSFORM) */

    /* Parameter adjustments */
    --trigs;
    --a;
    --work;

    /* Function Body */
    nh = *n / 2;
    nx = *n + 1;
    ink = *inc + *inc;

/*     A(0) AND A(N/2) */
    scale = (float)1. / (real) (*n);
    ia = 1;
    ib = 2;
    ja = 1;
    jb = *n * *inc + 1;
/* DIR$ IVDEP */
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	a[ja] = scale * (work[ia] + work[ib]);
	a[jb] = scale * (work[ia] - work[ib]);
	a[ja + *inc] = (float)0.;
	a[jb + *inc] = (float)0.;
	ia += nx;
	ib += nx;
	ja += *jump;
	jb += *jump;
/* L10: */
    }

/*     REMAINING WAVENUMBERS */
    scale *= (float).5;
    iabase = 3;
    ibbase = *n - 1;
    jabase = (*inc << 1) + 1;
    jbbase = (*n - 2) * *inc + 1;

    i__1 = nh;
    for (k = 3; k <= i__1; k += 2) {
	ia = iabase;
	ib = ibbase;
	ja = jabase;
	jb = jbbase;
	c__ = trigs[*n + k];
	s = trigs[*n + k + 1];
/* DIR$ IVDEP */
	i__2 = *lot;
	for (l = 1; l <= i__2; ++l) {
	    a[ja] = scale * (work[ia] + work[ib] + (c__ * (work[ia + 1] + 
		    work[ib + 1]) + s * (work[ia] - work[ib])));
	    a[jb] = scale * (work[ia] + work[ib] - (c__ * (work[ia + 1] + 
		    work[ib + 1]) + s * (work[ia] - work[ib])));
	    a[ja + *inc] = scale * (c__ * (work[ia] - work[ib]) - s * (work[
		    ia + 1] + work[ib + 1]) + (work[ib + 1] - work[ia + 1]));
	    a[jb + *inc] = scale * (c__ * (work[ia] - work[ib]) - s * (work[
		    ia + 1] + work[ib + 1]) - (work[ib + 1] - work[ia + 1]));
	    ia += nx;
	    ib += nx;
	    ja += *jump;
	    jb += *jump;
/* L20: */
	}
	iabase += 2;
	ibbase += -2;
	jabase += ink;
	jbbase -= ink;
/* L30: */
    }

    if (iabase != ibbase) {
	goto L50;
    }
/*     WAVENUMBER N/4 (IF IT EXISTS) */
    ia = iabase;
    ja = jabase;
    scale *= (float)2.;
/* DIR$ IVDEP */
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	a[ja] = scale * work[ia];
	a[ja + *inc] = -scale * work[ia + 1];
	ia += nx;
	ja += *jump;
/* L40: */
    }

L50:
    return 0;
} /* fft99b_ */

/* Subroutine */ int fft991_(a, work, trigs, ifax, inc, jump, n, lot, isign)
real *a, *work, *trigs;
integer *ifax, *inc, *jump, *n, *lot, *isign;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer nfax, i__, j, k, l, m, ibase, jbase;
    extern /* Subroutine */ int fft99a_(), fft99b_();
    static integer ia, ib, la, nh, nx;
    extern /* Subroutine */ int vpassm_();
    static integer igo, ink;

   static integer i;


/*     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC */
/*     FAST FOURIER TRANSFORM */

/*     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO */
/*     THAT IN MRFFT2 */

/*     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM */
/*     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12 */
/*     (1970), 315-337) */

/*     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA */
/*     WORK IS AN AREA OF SIZE (N+1)*LOT */
/*     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES */
/*     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2 */
/*     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR' */
/*         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA) */
/*     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR */
/*     N IS THE LENGTH OF THE DATA VECTORS */
/*     LOT IS THE NUMBER OF DATA VECTORS */
/*     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT */
/*           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL */

/*     ORDERING OF COEFFICIENTS: */
/*         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2) */
/*         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED */

/*     ORDERING OF DATA: */
/*         X(0),X(1),X(2),...,X(N-1) */

/*     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN */
/*     PARALLEL */

/*     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER */

/*     DEFINITION OF TRANSFORMS: */
/*     ------------------------- */

/*     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N)) */
/*         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K) */

/*     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N)) */
/*               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N)) */



    /* Parameter adjustments */
    --ifax;
    --trigs;
    --work;
    --a;
   
   
   


    /* Function Body */
    nfax = ifax[1];
    nx = *n + 1;
    nh = *n / 2;
    ink = *inc + *inc;

#if 0 
      printf("%d\n",*jump);
      printf("%d\n",*inc);
      exit(0);
#endif

    if (*isign == 1) {
	goto L30;
    }

/*     IF NECESSARY, TRANSFER DATA TO WORK AREA */
    igo = 50;
    if (nfax % 2 == 1) {
	goto L40;
    }
    ibase = 1;
    jbase = 1;
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	i__ = ibase;
	j = jbase;
/* DIR$ IVDEP */
	i__2 = *n;
	for (m = 1; m <= i__2; ++m) {
	    work[j] = a[i__];
	    i__ += *inc;
	    ++j;
/* L10: */
	}
	ibase += *jump;
	jbase += nx;
/* L20: */
    }
    
    igo = 60;
    goto L40;

/*     PREPROCESSING (ISIGN=+1) */
/*     ------------------------ */

L30:
    fft99a_(&a[1], &work[1], &trigs[1], inc, jump, n, lot);
    igo = 60;

/*     COMPLEX TRANSFORM */
/*     ----------------- */

L40:
    ia = 1;
    la = 1;
    i__1 = nfax;
    for (k = 1; k <= i__1; ++k) {
	if (igo == 60) {
	    goto L60;
	}
/* L50: */
	vpassm_(&a[ia], &a[ia + *inc], &work[1], &work[2], &trigs[1], &ink, &
		c__2, jump, &nx, lot, &nh, &ifax[k + 1], &la);
	igo = 60;
	goto L70;
L60:
	vpassm_(&work[1], &work[2], &a[ia], &a[ia + *inc], &trigs[1], &c__2, &
		ink, &nx, jump, lot, &nh, &ifax[k + 1], &la);
	igo = 50;
L70:
	la *= ifax[k + 1];
/* L80: */
    }

    if (*isign == -1) {
	goto L130;
    }

/*     IF NECESSARY, TRANSFER DATA FROM WORK AREA */
    if (nfax % 2 == 1) {
	goto L110;
    }
    ibase = 1;
    jbase = 1;
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	i__ = ibase;
	j = jbase;
/* DIR$ IVDEP */
	i__2 = *n;
	for (m = 1; m <= i__2; ++m) {
	    a[j] = work[i__];
	    ++i__;
	    j += *inc;
/* L90: */
	}
	ibase += nx;
	jbase += *jump;
/* L100: */
    }

/*     FILL IN ZEROS AT END */
L110:
    ib = *n * *inc + 1;
/* DIR$ IVDEP */
    i__1 = *lot;
    for (l = 1; l <= i__1; ++l) {
	a[ib] = (float)0.;
	a[ib + *inc] = (float)0.;
	ib += *jump;
/* L120: */
    }
    goto L140;

/*     POSTPROCESSING (ISIGN=-1): */
/*     -------------------------- */

L130:
    fft99b_(&work[1], &a[1], &trigs[1], inc, jump, n, lot);

L140:
    return 0;
} /* fft991_ */

/* Subroutine */ int set99_(trigs, ifax, n)
real *trigs;
integer *ifax, *n;
{
    /* Initialized data */

    static integer mode = 3;

#if 0
    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();
#endif

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int fftrig_(), fax_();

    /* Fortran I/O blocks */
    static cilist io___49 = { 0, 6, 0, 0, 0 };



/* MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE */
/* TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT */
/* DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE */
/* WAS WRITTEN. */

    /* Parameter adjustments */
    --ifax;
    --trigs;

    /* Function Body */
   
    fax_(&ifax[1], n, &mode);
    i__ = ifax[1];
    if (ifax[i__ + 1] > 5 || *n <= 4) {
	ifax[1] = -99;
    }
    if (ifax[1] <= 0) {
        printf("SET99 -- INVALID N");
        exit(0);
#if 0
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, " SET99 -- INVALID N", (ftnlen)19);
	e_wsle();
	s_stop("SET99", (ftnlen)5);
#endif
    }
    fftrig_(&trigs[1], n, &mode);
    return 0;
} /* set99_ */

/* Subroutine */ int fax_(ifax, n, mode)
integer *ifax, *n, *mode;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer nfax, item, i__, k, l, istop, ii, nn, inc;

    /* Parameter adjustments */
    --ifax;

    /* Function Body */

    
    nn = *n;
  
    if (nabs(*mode) == 1) {
	goto L10;
    }
    if (nabs(*mode) == 8) {
	goto L10;
    }
    nn = *n / 2;
    if (nn + nn == *n) {
	goto L10;
    }
    ifax[1] = -99;
    return 0;
L10:
    k = 1;
/*     TEST FOR FACTORS OF 4 */
L20:
    if (nn % 4 != 0) {
	goto L30;
    }
    ++k;
    ifax[k] = 4;
    nn /= 4;
    if (nn == 1) {
	goto L80;
    }
    goto L20;
/*     TEST FOR EXTRA FACTOR OF 2 */
L30:
    if (nn % 2 != 0) {
	goto L40;
    }
    ++k;
    ifax[k] = 2;
    nn /= 2;
    if (nn == 1) {
	goto L80;
    }
/*     TEST FOR FACTORS OF 3 */
L40:
    if (nn % 3 != 0) {
	goto L50;
    }
    ++k;
    ifax[k] = 3;
    nn /= 3;
    if (nn == 1) {
	goto L80;
    }
    goto L40;
/*     NOW FIND REMAINING FACTORS */
L50:
    l = 5;
    inc = 2;
/*     INC ALTERNATELY TAKES ON VALUES 2 AND 4 */
L60:
    if (nn % l != 0) {
	goto L70;
    }
    ++k;
    ifax[k] = l;
    nn /= l;
    if (nn == 1) {
	goto L80;
    }
    goto L60;
L70:
    l += inc;
    inc = 6 - inc;
    goto L60;
L80:
    ifax[1] = k - 1;
/*     IFAX(1) CONTAINS NUMBER OF FACTORS */
    nfax = ifax[1];
/*     SORT FACTORS INTO ASCENDING ORDER */
    if (nfax == 1) {
	goto L110;
    }
    i__1 = nfax;
    for (ii = 2; ii <= i__1; ++ii) {
	istop = nfax + 2 - ii;
	i__2 = istop;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    if (ifax[i__ + 1] >= ifax[i__]) {
		goto L90;
	    }
	    item = ifax[i__];
	    ifax[i__] = ifax[i__ + 1];
	    ifax[i__ + 1] = item;
L90:
	    ;
	}
/* L100: */
    }
L110:
    return 0;
} /* fax_ */

/* Subroutine */ int fftrig_(trigs, n, mode)
real *trigs;
integer *n, *mode;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double asin(), cos(), sin();

    /* Local variables */
    static integer i__, l;
    static real angle;
    static integer imode, la, nh;
    static real pi;
    static integer nn;
    static real del;

    /* Parameter adjustments */
    --trigs;

    /* Function Body */
    pi = asin((float)1.) * (float)2.;
    imode = nabs(*mode);
    nn = *n;
    if (imode > 1 && imode < 6) {
	nn = *n / 2;
    }
    del = (pi + pi) / (real) nn;
    l = nn + nn;
    i__1 = l;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	angle = (real) (i__ - 1) * (float).5 * del;
	trigs[i__] = cos(angle);
	trigs[i__ + 1] = sin(angle);
/* L10: */
    }
    if (imode == 1) {
	return 0;
    }
    if (imode == 8) {
	return 0;
    }
    del *= (float).5;
    nh = (nn + 1) / 2;
    l = nh + nh;
    la = nn + nn;
    i__1 = l;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	angle = (real) (i__ - 1) * (float).5 * del;
	trigs[la + i__] = cos(angle);
	trigs[la + i__ + 1] = sin(angle);
/* L20: */
    }
    if (imode <= 3) {
	return 0;
    }
    del *= (float).5;
    la += nn;
    if (*mode == 5) {
	goto L40;
    }
    i__1 = nn;
    for (i__ = 2; i__ <= i__1; ++i__) {
	angle = (real) (i__ - 1) * del;
	trigs[la + i__] = sin(angle) * (float)2.;
/* L30: */
    }
    return 0;
L40:
    del *= (float).5;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	angle = (real) (i__ - 1) * del;
	trigs[la + i__] = sin(angle);
/* L50: */
    }
    return 0;
} /* fftrig_ */

/* Subroutine */ int vpassm_(a, b, c__, d__, trigs, inc1, inc2, inc3, inc4, 
	lot, n, ifac, la)
real *a, *b, *c__, *d__, *trigs;
integer *inc1, *inc2, *inc3, *inc4, *lot, *n, *ifac, *la;
{
    /* Initialized data */

    static real sin36 = (float).587785252292473;
    static real cos36 = (float).809016994374947;
    static real sin72 = (float).951056516295154;
    static real cos72 = (float).309016994374947;
    static real sin60 = (float).866025403784437;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer iink, jink, jump, i__, j, k, l, m, ibase, jbase;
    static real c1, c2, c3, c4, s1, s2, s3, s4;
    static integer ia, ja, ib, jb, kb, ic, jc, kc, id, jd, kd, ie, je, ke, 
	    la1, ijk, igo;


/*     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA" */
/*     PERFORMS ONE PASS THROUGH DATA */
/*     AS PART OF MULTIPLE COMPLEX FFT ROUTINE */
/*     A IS FIRST REAL INPUT VECTOR */
/*     B IS FIRST IMAGINARY INPUT VECTOR */
/*     C IS FIRST REAL OUTPUT VECTOR */
/*     D IS FIRST IMAGINARY OUTPUT VECTOR */
/*     TRIGS IS PRECALCULATED TABLE OF SINES " COSINES */
/*     INC1 IS ADDRESSING INCREMENT FOR A AND B */
/*     INC2 IS ADDRESSING INCREMENT FOR C AND D */
/*     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S */
/*     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S */
/*     LOT IS THE NUMBER OF VECTORS */
/*     N IS LENGTH OF VECTORS */
/*     IFAC IS CURRENT FACTOR OF N */
/*     LA IS PRODUCT OF PREVIOUS FACTORS */

    /* Parameter adjustments */
    --trigs;
    --d__;
    --c__;
    --b;
    --a;

    /* Function Body */

    m = *n / *ifac;
    iink = m * *inc1;
    jink = *la * *inc2;
    jump = (*ifac - 1) * jink;
    ibase = 0;
    jbase = 0;
    igo = *ifac - 1;
    if (igo > 4) {
	return 0;
    }
    switch ((int)igo) {
	case 1:  goto L10;
	case 2:  goto L50;
	case 3:  goto L90;
	case 4:  goto L130;
    }

/*     CODING FOR FACTOR 2 */

L10:
    ia = 1;
    ja = 1;
    ib = ia + iink;
    jb = ja + jink;
    i__1 = *la;
    for (l = 1; l <= i__1; ++l) {
	i__ = ibase;
	j = jbase;
/* DIR$ IVDEP */
	i__2 = *lot;
	for (ijk = 1; ijk <= i__2; ++ijk) {
	    c__[ja + j] = a[ia + i__] + a[ib + i__];
	    d__[ja + j] = b[ia + i__] + b[ib + i__];
	    c__[jb + j] = a[ia + i__] - a[ib + i__];
	    d__[jb + j] = b[ia + i__] - b[ib + i__];
	    i__ += *inc3;
	    j += *inc4;
/* L15: */
	}
	ibase += *inc1;
	jbase += *inc2;
/* L20: */
    }
    if (*la == m) {
	return 0;
    }
    la1 = *la + 1;
    jbase += jump;
    i__1 = m;
    i__2 = *la;
    for (k = la1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
	kb = k + k - 2;
	c1 = trigs[kb + 1];
	s1 = trigs[kb + 2];
	i__3 = *la;
	for (l = 1; l <= i__3; ++l) {
	    i__ = ibase;
	    j = jbase;
/* DIR$ IVDEP */
	    i__4 = *lot;
	    for (ijk = 1; ijk <= i__4; ++ijk) {
		c__[ja + j] = a[ia + i__] + a[ib + i__];
		d__[ja + j] = b[ia + i__] + b[ib + i__];
		c__[jb + j] = c1 * (a[ia + i__] - a[ib + i__]) - s1 * (b[ia + 
			i__] - b[ib + i__]);
		d__[jb + j] = s1 * (a[ia + i__] - a[ib + i__]) + c1 * (b[ia + 
			i__] - b[ib + i__]);
		i__ += *inc3;
		j += *inc4;
/* L25: */
	    }
	    ibase += *inc1;
	    jbase += *inc2;
/* L30: */
	}
	jbase += jump;
/* L40: */
    }
    return 0;

/*     CODING FOR FACTOR 3 */

L50:
    ia = 1;
    ja = 1;
    ib = ia + iink;
    jb = ja + jink;
    ic = ib + iink;
    jc = jb + jink;
    i__2 = *la;
    for (l = 1; l <= i__2; ++l) {
	i__ = ibase;
	j = jbase;
/* DIR$ IVDEP */
	i__1 = *lot;
	for (ijk = 1; ijk <= i__1; ++ijk) {
	    c__[ja + j] = a[ia + i__] + (a[ib + i__] + a[ic + i__]);
	    d__[ja + j] = b[ia + i__] + (b[ib + i__] + b[ic + i__]);
	    c__[jb + j] = a[ia + i__] - (a[ib + i__] + a[ic + i__]) * (float)
		    .5 - sin60 * (b[ib + i__] - b[ic + i__]);
	    c__[jc + j] = a[ia + i__] - (a[ib + i__] + a[ic + i__]) * (float)
		    .5 + sin60 * (b[ib + i__] - b[ic + i__]);
	    d__[jb + j] = b[ia + i__] - (b[ib + i__] + b[ic + i__]) * (float)
		    .5 + sin60 * (a[ib + i__] - a[ic + i__]);
	    d__[jc + j] = b[ia + i__] - (b[ib + i__] + b[ic + i__]) * (float)
		    .5 - sin60 * (a[ib + i__] - a[ic + i__]);
	    i__ += *inc3;
	    j += *inc4;
/* L55: */
	}
	ibase += *inc1;
	jbase += *inc2;
/* L60: */
    }
    if (*la == m) {
	return 0;
    }
    la1 = *la + 1;
    jbase += jump;
    i__2 = m;
    i__1 = *la;
    for (k = la1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
	kb = k + k - 2;
	kc = kb + kb;
	c1 = trigs[kb + 1];
	s1 = trigs[kb + 2];
	c2 = trigs[kc + 1];
	s2 = trigs[kc + 2];
	i__3 = *la;
	for (l = 1; l <= i__3; ++l) {
	    i__ = ibase;
	    j = jbase;
/* DIR$ IVDEP */
	    i__4 = *lot;
	    for (ijk = 1; ijk <= i__4; ++ijk) {
		c__[ja + j] = a[ia + i__] + (a[ib + i__] + a[ic + i__]);
		d__[ja + j] = b[ia + i__] + (b[ib + i__] + b[ic + i__]);
		c__[jb + j] = c1 * (a[ia + i__] - (a[ib + i__] + a[ic + i__]) 
			* (float).5 - sin60 * (b[ib + i__] - b[ic + i__])) - 
			s1 * (b[ia + i__] - (b[ib + i__] + b[ic + i__]) * (
			float).5 + sin60 * (a[ib + i__] - a[ic + i__]));
		d__[jb + j] = s1 * (a[ia + i__] - (a[ib + i__] + a[ic + i__]) 
			* (float).5 - sin60 * (b[ib + i__] - b[ic + i__])) + 
			c1 * (b[ia + i__] - (b[ib + i__] + b[ic + i__]) * (
			float).5 + sin60 * (a[ib + i__] - a[ic + i__]));
		c__[jc + j] = c2 * (a[ia + i__] - (a[ib + i__] + a[ic + i__]) 
			* (float).5 + sin60 * (b[ib + i__] - b[ic + i__])) - 
			s2 * (b[ia + i__] - (b[ib + i__] + b[ic + i__]) * (
			float).5 - sin60 * (a[ib + i__] - a[ic + i__]));
		d__[jc + j] = s2 * (a[ia + i__] - (a[ib + i__] + a[ic + i__]) 
			* (float).5 + sin60 * (b[ib + i__] - b[ic + i__])) + 
			c2 * (b[ia + i__] - (b[ib + i__] + b[ic + i__]) * (
			float).5 - sin60 * (a[ib + i__] - a[ic + i__]));
		i__ += *inc3;
		j += *inc4;
/* L65: */
	    }
	    ibase += *inc1;
	    jbase += *inc2;
/* L70: */
	}
	jbase += jump;
/* L80: */
    }
    return 0;

/*     CODING FOR FACTOR 4 */

L90:
    ia = 1;
    ja = 1;
    ib = ia + iink;
    jb = ja + jink;
    ic = ib + iink;
    jc = jb + jink;
    id = ic + iink;
    jd = jc + jink;
    i__1 = *la;
    for (l = 1; l <= i__1; ++l) {
	i__ = ibase;
	j = jbase;
/* DIR$ IVDEP */
	i__2 = *lot;
	for (ijk = 1; ijk <= i__2; ++ijk) {
	    c__[ja + j] = a[ia + i__] + a[ic + i__] + (a[ib + i__] + a[id + 
		    i__]);
	    c__[jc + j] = a[ia + i__] + a[ic + i__] - (a[ib + i__] + a[id + 
		    i__]);
	    d__[ja + j] = b[ia + i__] + b[ic + i__] + (b[ib + i__] + b[id + 
		    i__]);
	    d__[jc + j] = b[ia + i__] + b[ic + i__] - (b[ib + i__] + b[id + 
		    i__]);
	    c__[jb + j] = a[ia + i__] - a[ic + i__] - (b[ib + i__] - b[id + 
		    i__]);
	    c__[jd + j] = a[ia + i__] - a[ic + i__] + (b[ib + i__] - b[id + 
		    i__]);
	    d__[jb + j] = b[ia + i__] - b[ic + i__] + (a[ib + i__] - a[id + 
		    i__]);
	    d__[jd + j] = b[ia + i__] - b[ic + i__] - (a[ib + i__] - a[id + 
		    i__]);
	    i__ += *inc3;
	    j += *inc4;
/* L95: */
	}
	ibase += *inc1;
	jbase += *inc2;
/* L100: */
    }
    if (*la == m) {
	return 0;
    }
    la1 = *la + 1;
    jbase += jump;
    i__1 = m;
    i__2 = *la;
    for (k = la1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
	kb = k + k - 2;
	kc = kb + kb;
	kd = kc + kb;
	c1 = trigs[kb + 1];
	s1 = trigs[kb + 2];
	c2 = trigs[kc + 1];
	s2 = trigs[kc + 2];
	c3 = trigs[kd + 1];
	s3 = trigs[kd + 2];
	i__3 = *la;
	for (l = 1; l <= i__3; ++l) {
	    i__ = ibase;
	    j = jbase;
/* DIR$ IVDEP */
	    i__4 = *lot;
	    for (ijk = 1; ijk <= i__4; ++ijk) {
		c__[ja + j] = a[ia + i__] + a[ic + i__] + (a[ib + i__] + a[id 
			+ i__]);
		d__[ja + j] = b[ia + i__] + b[ic + i__] + (b[ib + i__] + b[id 
			+ i__]);
		c__[jc + j] = c2 * (a[ia + i__] + a[ic + i__] - (a[ib + i__] 
			+ a[id + i__])) - s2 * (b[ia + i__] + b[ic + i__] - (
			b[ib + i__] + b[id + i__]));
		d__[jc + j] = s2 * (a[ia + i__] + a[ic + i__] - (a[ib + i__] 
			+ a[id + i__])) + c2 * (b[ia + i__] + b[ic + i__] - (
			b[ib + i__] + b[id + i__]));
		c__[jb + j] = c1 * (a[ia + i__] - a[ic + i__] - (b[ib + i__] 
			- b[id + i__])) - s1 * (b[ia + i__] - b[ic + i__] + (
			a[ib + i__] - a[id + i__]));
		d__[jb + j] = s1 * (a[ia + i__] - a[ic + i__] - (b[ib + i__] 
			- b[id + i__])) + c1 * (b[ia + i__] - b[ic + i__] + (
			a[ib + i__] - a[id + i__]));
		c__[jd + j] = c3 * (a[ia + i__] - a[ic + i__] + (b[ib + i__] 
			- b[id + i__])) - s3 * (b[ia + i__] - b[ic + i__] - (
			a[ib + i__] - a[id + i__]));
		d__[jd + j] = s3 * (a[ia + i__] - a[ic + i__] + (b[ib + i__] 
			- b[id + i__])) + c3 * (b[ia + i__] - b[ic + i__] - (
			a[ib + i__] - a[id + i__]));
		i__ += *inc3;
		j += *inc4;
/* L105: */
	    }
	    ibase += *inc1;
	    jbase += *inc2;
/* L110: */
	}
	jbase += jump;
/* L120: */
    }
    return 0;

/*     CODING FOR FACTOR 5 */

L130:
    ia = 1;
    ja = 1;
    ib = ia + iink;
    jb = ja + jink;
    ic = ib + iink;
    jc = jb + jink;
    id = ic + iink;
    jd = jc + jink;
    ie = id + iink;
    je = jd + jink;
    i__2 = *la;
    for (l = 1; l <= i__2; ++l) {
	i__ = ibase;
	j = jbase;
/* DIR$ IVDEP */
	i__1 = *lot;
	for (ijk = 1; ijk <= i__1; ++ijk) {
	    c__[ja + j] = a[ia + i__] + (a[ib + i__] + a[ie + i__]) + (a[ic + 
		    i__] + a[id + i__]);
	    d__[ja + j] = b[ia + i__] + (b[ib + i__] + b[ie + i__]) + (b[ic + 
		    i__] + b[id + i__]);
	    c__[jb + j] = a[ia + i__] + cos72 * (a[ib + i__] + a[ie + i__]) - 
		    cos36 * (a[ic + i__] + a[id + i__]) - (sin72 * (b[ib + 
		    i__] - b[ie + i__]) + sin36 * (b[ic + i__] - b[id + i__]))
		    ;
	    c__[je + j] = a[ia + i__] + cos72 * (a[ib + i__] + a[ie + i__]) - 
		    cos36 * (a[ic + i__] + a[id + i__]) + (sin72 * (b[ib + 
		    i__] - b[ie + i__]) + sin36 * (b[ic + i__] - b[id + i__]))
		    ;
	    d__[jb + j] = b[ia + i__] + cos72 * (b[ib + i__] + b[ie + i__]) - 
		    cos36 * (b[ic + i__] + b[id + i__]) + (sin72 * (a[ib + 
		    i__] - a[ie + i__]) + sin36 * (a[ic + i__] - a[id + i__]))
		    ;
	    d__[je + j] = b[ia + i__] + cos72 * (b[ib + i__] + b[ie + i__]) - 
		    cos36 * (b[ic + i__] + b[id + i__]) - (sin72 * (a[ib + 
		    i__] - a[ie + i__]) + sin36 * (a[ic + i__] - a[id + i__]))
		    ;
	    c__[jc + j] = a[ia + i__] - cos36 * (a[ib + i__] + a[ie + i__]) + 
		    cos72 * (a[ic + i__] + a[id + i__]) - (sin36 * (b[ib + 
		    i__] - b[ie + i__]) - sin72 * (b[ic + i__] - b[id + i__]))
		    ;
	    c__[jd + j] = a[ia + i__] - cos36 * (a[ib + i__] + a[ie + i__]) + 
		    cos72 * (a[ic + i__] + a[id + i__]) + (sin36 * (b[ib + 
		    i__] - b[ie + i__]) - sin72 * (b[ic + i__] - b[id + i__]))
		    ;
	    d__[jc + j] = b[ia + i__] - cos36 * (b[ib + i__] + b[ie + i__]) + 
		    cos72 * (b[ic + i__] + b[id + i__]) + (sin36 * (a[ib + 
		    i__] - a[ie + i__]) - sin72 * (a[ic + i__] - a[id + i__]))
		    ;
	    d__[jd + j] = b[ia + i__] - cos36 * (b[ib + i__] + b[ie + i__]) + 
		    cos72 * (b[ic + i__] + b[id + i__]) - (sin36 * (a[ib + 
		    i__] - a[ie + i__]) - sin72 * (a[ic + i__] - a[id + i__]))
		    ;
	    i__ += *inc3;
	    j += *inc4;
/* L135: */
	}
	ibase += *inc1;
	jbase += *inc2;
/* L140: */
    }
    if (*la == m) {
	return 0;
    }
    la1 = *la + 1;
    jbase += jump;
    i__2 = m;
    i__1 = *la;
    for (k = la1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
	kb = k + k - 2;
	kc = kb + kb;
	kd = kc + kb;
	ke = kd + kb;
	c1 = trigs[kb + 1];
	s1 = trigs[kb + 2];
	c2 = trigs[kc + 1];
	s2 = trigs[kc + 2];
	c3 = trigs[kd + 1];
	s3 = trigs[kd + 2];
	c4 = trigs[ke + 1];
	s4 = trigs[ke + 2];
	i__3 = *la;
	for (l = 1; l <= i__3; ++l) {
	    i__ = ibase;
	    j = jbase;
/* DIR$ IVDEP */
	    i__4 = *lot;
	    for (ijk = 1; ijk <= i__4; ++ijk) {
		c__[ja + j] = a[ia + i__] + (a[ib + i__] + a[ie + i__]) + (a[
			ic + i__] + a[id + i__]);
		d__[ja + j] = b[ia + i__] + (b[ib + i__] + b[ie + i__]) + (b[
			ic + i__] + b[id + i__]);
		c__[jb + j] = c1 * (a[ia + i__] + cos72 * (a[ib + i__] + a[ie 
			+ i__]) - cos36 * (a[ic + i__] + a[id + i__]) - (
			sin72 * (b[ib + i__] - b[ie + i__]) + sin36 * (b[ic + 
			i__] - b[id + i__]))) - s1 * (b[ia + i__] + cos72 * (
			b[ib + i__] + b[ie + i__]) - cos36 * (b[ic + i__] + b[
			id + i__]) + (sin72 * (a[ib + i__] - a[ie + i__]) + 
			sin36 * (a[ic + i__] - a[id + i__])));
		d__[jb + j] = s1 * (a[ia + i__] + cos72 * (a[ib + i__] + a[ie 
			+ i__]) - cos36 * (a[ic + i__] + a[id + i__]) - (
			sin72 * (b[ib + i__] - b[ie + i__]) + sin36 * (b[ic + 
			i__] - b[id + i__]))) + c1 * (b[ia + i__] + cos72 * (
			b[ib + i__] + b[ie + i__]) - cos36 * (b[ic + i__] + b[
			id + i__]) + (sin72 * (a[ib + i__] - a[ie + i__]) + 
			sin36 * (a[ic + i__] - a[id + i__])));
		c__[je + j] = c4 * (a[ia + i__] + cos72 * (a[ib + i__] + a[ie 
			+ i__]) - cos36 * (a[ic + i__] + a[id + i__]) + (
			sin72 * (b[ib + i__] - b[ie + i__]) + sin36 * (b[ic + 
			i__] - b[id + i__]))) - s4 * (b[ia + i__] + cos72 * (
			b[ib + i__] + b[ie + i__]) - cos36 * (b[ic + i__] + b[
			id + i__]) - (sin72 * (a[ib + i__] - a[ie + i__]) + 
			sin36 * (a[ic + i__] - a[id + i__])));
		d__[je + j] = s4 * (a[ia + i__] + cos72 * (a[ib + i__] + a[ie 
			+ i__]) - cos36 * (a[ic + i__] + a[id + i__]) + (
			sin72 * (b[ib + i__] - b[ie + i__]) + sin36 * (b[ic + 
			i__] - b[id + i__]))) + c4 * (b[ia + i__] + cos72 * (
			b[ib + i__] + b[ie + i__]) - cos36 * (b[ic + i__] + b[
			id + i__]) - (sin72 * (a[ib + i__] - a[ie + i__]) + 
			sin36 * (a[ic + i__] - a[id + i__])));
		c__[jc + j] = c2 * (a[ia + i__] - cos36 * (a[ib + i__] + a[ie 
			+ i__]) + cos72 * (a[ic + i__] + a[id + i__]) - (
			sin36 * (b[ib + i__] - b[ie + i__]) - sin72 * (b[ic + 
			i__] - b[id + i__]))) - s2 * (b[ia + i__] - cos36 * (
			b[ib + i__] + b[ie + i__]) + cos72 * (b[ic + i__] + b[
			id + i__]) + (sin36 * (a[ib + i__] - a[ie + i__]) - 
			sin72 * (a[ic + i__] - a[id + i__])));
		d__[jc + j] = s2 * (a[ia + i__] - cos36 * (a[ib + i__] + a[ie 
			+ i__]) + cos72 * (a[ic + i__] + a[id + i__]) - (
			sin36 * (b[ib + i__] - b[ie + i__]) - sin72 * (b[ic + 
			i__] - b[id + i__]))) + c2 * (b[ia + i__] - cos36 * (
			b[ib + i__] + b[ie + i__]) + cos72 * (b[ic + i__] + b[
			id + i__]) + (sin36 * (a[ib + i__] - a[ie + i__]) - 
			sin72 * (a[ic + i__] - a[id + i__])));
		c__[jd + j] = c3 * (a[ia + i__] - cos36 * (a[ib + i__] + a[ie 
			+ i__]) + cos72 * (a[ic + i__] + a[id + i__]) + (
			sin36 * (b[ib + i__] - b[ie + i__]) - sin72 * (b[ic + 
			i__] - b[id + i__]))) - s3 * (b[ia + i__] - cos36 * (
			b[ib + i__] + b[ie + i__]) + cos72 * (b[ic + i__] + b[
			id + i__]) - (sin36 * (a[ib + i__] - a[ie + i__]) - 
			sin72 * (a[ic + i__] - a[id + i__])));
		d__[jd + j] = s3 * (a[ia + i__] - cos36 * (a[ib + i__] + a[ie 
			+ i__]) + cos72 * (a[ic + i__] + a[id + i__]) + (
			sin36 * (b[ib + i__] - b[ie + i__]) - sin72 * (b[ic + 
			i__] - b[id + i__]))) + c3 * (b[ia + i__] - cos36 * (
			b[ib + i__] + b[ie + i__]) + cos72 * (b[ic + i__] + b[
			id + i__]) - (sin36 * (a[ib + i__] - a[ie + i__]) - 
			sin72 * (a[ic + i__] - a[id + i__])));
		i__ += *inc3;
		j += *inc4;
/* L145: */
	    }
	    ibase += *inc1;
	    jbase += *inc2;
/* L150: */
	}
	jbase += jump;
/* L160: */
    }
    return 0;
} /* vpassm_ */



integer nabs(integer x) {return x > 0 ? x : -x;}
