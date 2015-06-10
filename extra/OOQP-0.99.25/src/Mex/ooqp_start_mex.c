#include "mexUtility.h"
#include "cQpGenSparse.h"
#include "string.h"
#include "OoqpMonitorData.h"

#define INT4

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

enum TerminationCode 
{
  SUCCESSFUL_TERMINATION = 0,
  NOT_FINISHED,
  MAX_ITS_EXCEEDED,
  INFEASIBLE,
  UNKNOWN
};

int mexDoMonitor( OoqpMonitorData * data );
void QpGenFinishWithSp ( double   x[],  
			 double   w[], double gamma[],
			 double   v[], double phi[],
			 double   y[], 
			 double   z[],  
			 double   s[],
			 double   t[], double lambda[],
			 double   u[], double pi[],
			 double mutol, double artol,
			 QpGenContext * ctx, int * status_code );

/* ****************************************************************
 * mex stub function: mexFunction
 * **************************************************************** */
void mexFunction( int nlhs, mxArray * plhs[], int nrhs,
    const mxArray * prhs[] ) {

    double      *dummyDoublePtr;

    double      *c;                       /* Input    argument 0   */
    double      *Q;                       /* Input    argument 1   */
    int         *Q_ir;                   
    int         *Q_jc;                   
    int          Q_nz;                   
    double      *xlow;                    /* Input    argument 2   */
    double      *ixlow;                   /* Input    argument 3   */
    double      *xupp;                    /* Input    argument 4   */
    double      *ixupp;                   /* Input    argument 5   */
    double      *A;                       /* Input    argument 6   */
    int         *A_ir;                   
    int         *A_jc;                   
    int          A_nz;                   
    double      *dA;                      /* Input    argument 7   */
    double      *C;                       /* Input    argument 8   */
    int         *C_ir;                   
    int         *C_jc;                   
    int          C_nz;                   
    double      *clow;                    /* Input    argument 9   */
    double      *iclow;                   /* Input    argument 10  */
    double      *cupp;                    /* Input    argument 11  */
    double      *icupp;                   /* Input    argument 12  */
    double      *xi;                      /* Input    argument 13  */
    double      *si;                      /* Input    argument 14  */
    double      *yi;                      /* Input    argument 15  */
    double      *zi;                      /* Input    argument 16  */
    double      *vi;                      /* Input    argument 17  */
    double      *gammai;                  /* Input    argument 18  */
    double      *wi;                      /* Input    argument 19  */
    double      *phii;                    /* Input    argument 20  */
    double      *ti;                      /* Input    argument 21  */
    double      *lambdai;                 /* Input    argument 22  */
    double      *ui;                      /* Input    argument 23  */
    double      *pii;                     /* Input    argument 24  */

    int          status_code;             /* Output   argument 0   */
    mxArray     *mexstatus_code;         
    double      *x;                       /* Output   argument 1   */
    mxArray     *mexx;                   
    double      *gamma;                   /* Output   argument 2   */
    mxArray     *mexgamma;               
    double      *phi;                     /* Output   argument 3   */
    mxArray     *mexphi;                 
    double      *y;                       /* Output   argument 4   */
    mxArray     *mexy;                   
    double      *z;                       /* Output   argument 5   */
    mxArray     *mexz;                   
    double      *lambda;                  /* Output   argument 6   */
    mxArray     *mexlambda;              
    double      *pi;                      /* Output   argument 7   */
    mxArray     *mexpi;                  

    char         doPrint[32];             /* Optional argument 0   */
    double       mutol = 1e-8;            /* Optional argument 1   */
    double       rtol  = 1e-8;            /* Optional argument 2   */

    int         nx;                       /* Array dimension     */
    int         my;                       /* Array dimension     */
    int         mz;                       /* Array dimension     */

    int         i;
    char        *cxlow = NULL, *cxupp = NULL, *cclow = NULL, *ccupp = NULL;
    int         iPrint;

    QpGenContext ctx;
    OoqpMonitorData monitorCtx;

#ifdef HAVE_GETRUSAGE
	double solve_time;
	struct rusage before_solve, after_solve;
#endif

    if ( nrhs < 25 || nrhs > 28  ) {
        mexErrMsgTxt( "Wrong number of input arguments" );
    }
    if ( nlhs > 8 ) {
        mexErrMsgTxt( "Too many output arguments" );
    }

    /* c : Input argument 0 */
    nx = mxGetM( prhs[0] );
    assertColDim( prhs[0], 1, "c" );
    assertDense( prhs[0], "c" );
    c = mxGetPr( prhs[0] );

    /* Q : Input argument 1 */
    assertRowDim( prhs[1], nx, "Q" );
    assertColDim( prhs[1], nx, "Q" );
    assertSparse( prhs[1], "Q" );
    Q = mxGetPr( prhs[1] );
    Q_ir = mxGetIr( prhs[1] );
    Q_jc = mxGetJc( prhs[1] );
    Q_nz  = Q_jc[nx];
    /* sparse data will be copied */

    /* xlow : Input argument 2 */
    assertRowDim( prhs[2], nx, "xlow" );
    assertColDim( prhs[2], 1, "xlow" );
    assertDense( prhs[2], "xlow" );
    xlow = mxGetPr( prhs[2] );

    /* ixlow : Input argument 3 */
    assertRowDim( prhs[3], nx, "ixlow" );
    assertColDim( prhs[3], 1, "ixlow" );
    assertDense( prhs[3], "ixlow" );
    ixlow = mxGetPr( prhs[3] );

    /* xupp : Input argument 4 */
    assertRowDim( prhs[4], nx, "xupp" );
    assertColDim( prhs[4], 1, "xupp" );
    assertDense( prhs[4], "xupp" );
    xupp = mxGetPr( prhs[4] );

    /* ixupp : Input argument 5 */
    assertRowDim( prhs[5], nx, "ixupp" );
    assertColDim( prhs[5], 1, "ixupp" );
    assertDense( prhs[5], "ixupp" );
    ixupp = mxGetPr( prhs[5] );

    /* A : Input argument 6 */
    assertRowDim( prhs[6], nx, "A" );
    my = mxGetN( prhs[6] );
    assertSparse( prhs[6], "A" );
    A = mxGetPr( prhs[6] );
    A_ir = mxGetIr( prhs[6] );
    A_jc = mxGetJc( prhs[6] );
    A_nz  = A_jc[my];
    /* sparse data will be copied */

    /* dA : Input argument 7 */
    assertRowDim( prhs[7], my, "dA" );
    assertColDim( prhs[7], 1, "dA" );
    assertDense( prhs[7], "dA" );
    dA = mxGetPr( prhs[7] );

    /* C : Input argument 8 */
    assertRowDim( prhs[8], nx, "C" );
    mz = mxGetN( prhs[8] );
    assertSparse( prhs[8], "C" );
    C = mxGetPr( prhs[8] );
    C_ir = mxGetIr( prhs[8] );
    C_jc = mxGetJc( prhs[8] );
    C_nz  = C_jc[mz];
    /* sparse data will be copied */

    /* clow : Input argument 9 */
    assertRowDim( prhs[9], mz, "clow" );
    assertColDim( prhs[9], 1, "clow" );
    assertDense( prhs[9], "clow" );
    clow = mxGetPr( prhs[9] );

    /* iclow : Input argument 10 */
    assertRowDim( prhs[10], mz, "iclow" );
    assertColDim( prhs[10], 1, "iclow" );
    assertDense( prhs[10], "iclow" );
    iclow = mxGetPr( prhs[10] );

    /* cupp : Input argument 11 */
    assertRowDim( prhs[11], mz, "cupp" );
    assertColDim( prhs[11], 1, "cupp" );
    assertDense( prhs[11], "cupp" );
    cupp = mxGetPr( prhs[11] );

    /* icupp : Input argument 12 */
    assertRowDim( prhs[12], mz, "icupp" );
    assertColDim( prhs[12], 1, "icupp" );
    assertDense( prhs[12], "icupp" );
    icupp = mxGetPr( prhs[12] );

    /* xi : Input argument 13 */
    assertRowDim( prhs[13], nx, "xi" );
    assertColDim( prhs[13], 1, "xi" );
    assertDense( prhs[13], "xi" );
    xi = mxGetPr( prhs[13] );

    /* si : Input argument 14 */
    assertRowDim( prhs[14], mz, "si" );
    assertColDim( prhs[14], 1, "si" );
    assertDense( prhs[14], "si" );
    si = mxGetPr( prhs[14] );

    /* yi : Input argument 15 */
    assertRowDim( prhs[15], my, "yi" );
    assertColDim( prhs[15], 1, "yi" );
    assertDense( prhs[15], "yi" );
    yi = mxGetPr( prhs[15] );

    /* zi : Input argument 16 */
    assertRowDim( prhs[16], mz, "zi" );
    assertColDim( prhs[16], 1, "zi" );
    assertDense( prhs[16], "zi" );
    zi = mxGetPr( prhs[16] );

    /* vi : Input argument 17 */
    assertRowDim( prhs[17], nx, "vi" );
    assertColDim( prhs[17], 1, "vi" );
    assertDense( prhs[17], "vi" );
    vi = mxGetPr( prhs[17] );

    /* gammai : Input argument 18 */
    assertRowDim( prhs[18], nx, "gammai" );
    assertColDim( prhs[18], 1, "gammai" );
    assertDense( prhs[18], "gammai" );
    gammai = mxGetPr( prhs[18] );

    /* wi : Input argument 19 */
    assertRowDim( prhs[19], nx, "wi" );
    assertColDim( prhs[19], 1, "wi" );
    assertDense( prhs[19], "wi" );
    wi = mxGetPr( prhs[19] );

    /* phii : Input argument 20 */
    assertRowDim( prhs[20], nx, "phii" );
    assertColDim( prhs[20], 1, "phii" );
    assertDense( prhs[20], "phii" );
    phii = mxGetPr( prhs[20] );

    /* ti : Input argument 21 */
    assertRowDim( prhs[21], mz, "ti" );
    assertColDim( prhs[21], 1, "ti" );
    assertDense( prhs[21], "ti" );
    ti = mxGetPr( prhs[21] );

    /* lambdai : Input argument 22 */
    assertRowDim( prhs[22], mz, "lambdai" );
    assertColDim( prhs[22], 1, "lambdai" );
    assertDense( prhs[22], "lambdai" );
    lambdai = mxGetPr( prhs[22] );

    /* ui : Input argument 23 */
    assertRowDim( prhs[23], mz, "ui" );
    assertColDim( prhs[23], 1, "ui" );
    assertDense( prhs[23], "ui" );
    ui = mxGetPr( prhs[23] );

    /* pii : Input argument 24 */
    assertRowDim( prhs[24], mz, "pii" );
    assertColDim( prhs[24], 1, "pii" );
    assertDense( prhs[24], "pii" );
    pii = mxGetPr( prhs[24] );

    /* Output argument 0, "status_code", is scalar */

    /* x : output argument 1 */
    mexx = mxCreateDoubleMatrix( nx, 1, mxREAL );
    x =  mxGetPr( mexx );

    /* gamma : output argument 2 */
    mexgamma = mxCreateDoubleMatrix( nx, 1, mxREAL );
    gamma =  mxGetPr( mexgamma );

    /* phi : output argument 3 */
    mexphi = mxCreateDoubleMatrix( nx, 1, mxREAL );
    phi =  mxGetPr( mexphi );

    /* y : output argument 4 */
    mexy = mxCreateDoubleMatrix( my, 1, mxREAL );
    y =  mxGetPr( mexy );

    /* z : output argument 5 */
    mexz = mxCreateDoubleMatrix( mz, 1, mxREAL );
    z =  mxGetPr( mexz );

    /* lambda : output argument 6 */
    mexlambda = mxCreateDoubleMatrix( mz, 1, mxREAL );
    lambda =  mxGetPr( mexlambda );

    /* pi : output argument 7 */
    mexpi = mxCreateDoubleMatrix( mz, 1, mxREAL );
    pi =  mxGetPr( mexpi );

    if ( 25 < nrhs ) {
        /* doPrint : Optional argument 25 */
        assertString( prhs[25], 32, "doPrint" );
        mx2char( 1, 32, mxGetN(prhs[25]), 
                 (mxChar *) mxGetPr(prhs[25]), doPrint );
    }
    if ( 26 < nrhs ) {
        /* mutol : Optional argument 26 */
        assertScalar( prhs[26], "mutol" );
        mutol =  *mxGetPr ( prhs[26] );
    }
    if ( 27 < nrhs ) {
        /* rtol : Optional argument 27 */
        assertScalar( prhs[27], "rtol" );
        rtol =  *mxGetPr ( prhs[27] );
    }


    /* Check the optional argument */
    if( 25 >= nrhs ) {
      /* there is no optional argument */
      iPrint = 0;
    } else {
      iPrint = -1;
      if( 0 == strncmp( "yes                             ",
			doPrint, 32 ) ) iPrint = 1;
      if( 0 == strncmp( "on                              ",
			doPrint, 32 ) ) iPrint = 1;
      if( 0 == strncmp( "no                              ",
			doPrint, 32 ) ) iPrint = 0;
      if( 0 == strncmp( "off                             ",
			doPrint, 32 ) ) iPrint = 0;
      if( iPrint == -1 ) {
	mexErrMsgTxt( "Valid values of argument 13 \"doPrint\" are"
		      " \"yes\", \"no\", \"on\" or \"off\"\n" );
      }
    }

    cxlow = (char *) mxCalloc( nx, sizeof( char ) );
    cxupp = (char *) mxCalloc( nx, sizeof( char ) );
    for( i = 0; i < nx; i++ ) {
      cxlow[i] = (ixlow[i] == 0 ) ? 0 : 1;
      cxupp[i] = (ixupp[i] == 0 ) ? 0 : 1;
    }
    if( mz > 0 ) {
      cclow = (char *) mxCalloc( mz, sizeof( char ) );
      ccupp = (char *) mxCalloc( mz, sizeof( char ) );
      for( i = 0; i < mz; i++ ) {
	cclow[i] = (iclow[i] == 0) ? 0 : 1;
	ccupp[i] = (icupp[i] == 0) ? 0 : 1;
      }
    }

    QpGenHbGondzioSetup( c,      nx,
			 Q_jc,   Q_ir,  Q,
			 xlow,   cxlow, 
			 xupp,   cxupp,
			 A_jc,   my,    A_ir, A, dA,
			 C_jc,   mz,    C_ir, C,
			 clow,   cclow,
			 cupp,  ccupp,
			 &ctx,
			 &status_code );

    if( iPrint ) QpGenAddMonitor( &ctx, mexDoMonitor, &monitorCtx );
 
    if( status_code != 0 )
      mexErrMsgTxt( "Could not begin solving the problem\n." );
    
    for( i = 0; i < nx; i++ ) {
      x[i] = xi[i]; gamma[i] = gammai[i], phi[i] = phii[i];
    }
    for( i = 0; i < my; i++ ) {
      y[i] = yi[i];
    }
    for( i = 0; i < mz; i++ ) {
      z[i] = zi[i]; lambda[i] = lambdai[i]; pi[i]     = pii[i];
    }
#ifdef HAVE_GETRUSAGE
    getrusage( RUSAGE_SELF, &before_solve );
#endif
    QpGenFinishWithSp( x,  vi, gamma, wi, phi, y, z, 
		       si, ti, lambda, ui, pi, mutol, rtol,
		       &ctx, &status_code );
#ifdef HAVE_GETRUSAGE
	getrusage( RUSAGE_SELF, &after_solve );

	solve_time =
	  (after_solve.ru_utime.tv_sec - before_solve.ru_utime.tv_sec)
	  + (after_solve.ru_utime.tv_usec - before_solve.ru_utime.tv_usec)
	  / 1000000.0;

	if( iPrint &&  status_code == 0 ) 
	  mexPrintf( "QP solved in %10.3g seconds.\n", solve_time );
#endif

    if( status_code == 1024 ) 
      mexErrMsgTxt( "The given starting point was not interior" );

    if( status_code != 0 ) mexWarnMsgTxt( "Could not solve the problem\n." );
    QpGenCleanup( &ctx );


    /* Free any scratch arrays */
   if( ccupp != NULL ) mxFree( ccupp ); ccupp = NULL;
    if( cclow != NULL ) mxFree( cclow ); cclow = NULL;
    mxFree( cxupp ); cxupp = NULL;
    mxFree( cxlow ); cxlow = NULL;

    /* status_code : output argument 0 */
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    dummyDoublePtr = mxGetPr( plhs[0] );
    *dummyDoublePtr = (double) status_code;

    /* x : output argument 1 */
    if ( 1 < nlhs ) {
        plhs[1] = mexx;
    } else {
        mxDestroyArray( mexx );
    }

    /* gamma : output argument 2 */
    if ( 2 < nlhs ) {
        plhs[2] = mexgamma;
    } else {
        mxDestroyArray( mexgamma );
    }

    /* phi : output argument 3 */
    if ( 3 < nlhs ) {
        plhs[3] = mexphi;
    } else {
        mxDestroyArray( mexphi );
    }

    /* y : output argument 4 */
    if ( 4 < nlhs ) {
        plhs[4] = mexy;
    } else {
        mxDestroyArray( mexy );
    }

    /* z : output argument 5 */
    if ( 5 < nlhs ) {
        plhs[5] = mexz;
    } else {
        mxDestroyArray( mexz );
    }

    /* lambda : output argument 6 */
    if ( 6 < nlhs ) {
        plhs[6] = mexlambda;
    } else {
        mxDestroyArray( mexlambda );
    }

    /* pi : output argument 7 */
    if ( 7 < nlhs ) {
        plhs[7] = mexpi;
    } else {
        mxDestroyArray( mexpi );
    }

}
/* ****************************************************************
 * end of mex stub: mexFunction
 * **************************************************************** */

int mexDoMonitor( OoqpMonitorData * data )
{
  switch( data->level ) {
  case 0 : case 1: { 
    mexPrintf( "\nDuality Gap: %g\n", data->gap );
	if( data->i > 1 ) {
	  mexPrintf( " alpha = %g\n", data->alpha );
	}
    mexPrintf( " *** Iteration %d ***\n", data->i );
	mexPrintf( " mu = %9.4e relative residual norm = %9.4e\n",
			   data->mu,
			   data->rnorm / data->dataNorm );

    if( data->level == 1) { 
      /* Termination has been detected by the status check; print
       * appropriate message */
      switch( data->status_code ) {
      case SUCCESSFUL_TERMINATION:
		mexPrintf( "\n *** SUCCESSFUL TERMINATION ***\n" );
		break;
      case MAX_ITS_EXCEEDED:
		mexPrintf( "\n *** MAXIMUM ITERATIONS REACHED *** \n" );
		break;
      case INFEASIBLE:
		mexPrintf( "\n *** TERMINATION: PROBABLY INFEASIBLE *** \n" );
		break;
	  default:
		mexPrintf( "\n *** TERMINATION: STATUS UNKNOWN *** \n" );
		break;
      } /* end switch(statusCode) */
    }
  } break; /* end case 0: case 1: */
  } /* end switch(level) */
  
  return 0;
}

