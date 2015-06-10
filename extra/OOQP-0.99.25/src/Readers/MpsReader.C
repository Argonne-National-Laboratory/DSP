/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "MpsReader.h"
#include "OoqpVector.h"
#include "DoubleMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include <cstring>
#include <cerrno>
#include <cassert>
#include <cstdlib>
#include <cctype>
#include <cmath>

extern int gOoqpPrintLevel;

enum{ DATALINE = 1, HEADERLINE };
enum{ kBadRowType = -1, kFreeRow, kLessRow, kGreaterRow, kEqualRow,
        kLessRowWithRange, kGreaterRowWithRange };
enum { kLowerBound, kUpperBound, kFixedBound, kFreeBound, kMInftyBound, kPInftyBound };
 
const int READERROR = mpsioerr;

struct MpsRowInfo {
  char name[17];
  int  kind;
  int  nnz;
};

struct MpsColInfo {
  char name[17];
  int  nnz;
};

int MpsRowTypeFromCode( char code[4] );
int MpsRowTypeFromCode2( char code );

void doubleLexSort( int first[], int n, int second[], double data[] );

MpsReader::MpsReader( FILE * file_ )
{
  file  = file_;
  iline = 0;
  rowInfo     = 0;
  rowRemap    = 0;
  colInfo     = 0;
  rowTable    = 0;
  colTable    = 0;
  infilename  =  0;
  // -1 indicates that we have not yet determined these values
  nnzA = -1; nnzC = -1; nnzQ = -1;
  my   = -1; mz   = -1;
  memcpy(this->objectiveSense, "MIN", 3);
}

int isOnlySpaces( char str[], int start, int finish )
{
  int i;
  for( i = start; i <= finish; i++ ) {
    if( str[i] != ' ' ) return 0;
  }
  return 1;
}

// isLJustOf - is str the left justification of prefix is a field of len?
int isLJustOf( const char str[], const char prefix[], int len )
{
  int i;
  int result = 1;
  for( i = 0; i < len; i++ ) {
    if( prefix[i] == '\0' ) break;
    if( prefix[i] != str[i] ) {
      result = 0;
      break;
    }
  }
  // Did we use all characters of the prefix?
  if ( prefix[i] != '\0' ) {
    // No we didn't, so str is not the ljust of prefix
    result = 0;
  }
  if( result == 1 ) { 
    // Ok so far. Make sure all remaining characters in str are ' '
    for( ; i < len; i++ ) { // for all remaining characters in str
      if( str[i] != ' ' ) {
	result = 0;
	break;
      }
    } // end for all remaining characters in str
  } // end if Ok so far
  return result;
}

double asDouble( char str[], int len, int& ierr )
{
  char * pstr = str;
  int lpstr   = len;
  char * endptr;
  ierr = 0;

  while( ' ' == *pstr && lpstr > 0 ) {
    pstr++; lpstr--;
  }
  if( lpstr == 0 ) {
    ierr = 1;
    return 0.0;
  }

  double value = strtod( pstr, &endptr );
  while( *endptr == ' ' ) endptr++;
  if( endptr != pstr + lpstr ) {
    ierr = 1;
    return 0.0;
  }
  
  return value;
}

void MpsReader::readColsSection( OoqpVector& c_,
				 GenMatrix& A, GenMatrix& C,
				 char line[],
                                 int& ierr, int& kindOfLine )
{
  // Create a few temporaries
  if( nnzA < 0 || nnzC < 0 ) {
    int nnzA_, nnzC_, nnzQ_; // Force the computation of these cached values
    this->numberOfNonZeros( nnzQ_, nnzA_, nnzC_ );
  }

  SimpleVectorHandle c( new SimpleVector( totalCols ) );

  int    *irowC = 0, *jcolC = 0;
  double *dC = 0;
  if( nnzC > 0 ) {
    irowC = new int[nnzC];
    jcolC = new int[nnzC];
    dC    = new double[nnzC];
  }
  int    *irowA = 0, *jcolA = 0;
  double *dA    = 0;
  if( nnzA > 0 ) {
    irowA = new int[nnzA];
    jcolA = new int[nnzA];
    dA    = new double[nnzA];
  }
  this->readColsSection( c->elements(), irowA, jcolA, dA, irowC, jcolC, dC,
			 line, ierr, kindOfLine );

  this->stuffMatrix( A, irowA, nnzA, jcolA, dA );
  this->stuffMatrix( C, irowC, nnzC, jcolC, dC );

  c_.copyFrom( *c );

  delete [] dA; delete [] jcolA; delete [] irowA;
  delete [] dC; delete [] jcolC; delete [] irowC;
}

void MpsReader::readColsSection( double c[],
				 int irowA[], int jcolA[], double dA[],
				 int irowC[], int jcolC[], double dC[],
				 char line[],
                                 int& ierr, int& kindOfLine )
{
  char   blank[4], colname[16], row[2][16];
  double val[2];
  int    hasSecondValue;

  // reset the column position
  ierr = fseek( file, columnFilePosition, SEEK_SET );
  if( ierr != 0 ) {
    ierr = mpsioerr;
    return;
  }
  iline = firstColumnLine;

  assert( rowRemap ); // We have already mapped the rownums to their position
                      // in A and C
  int colnum = -1;
  char oldColumnName[16] = "";

  int nea = 0; // we have not yet read any elements of A
  int nec = 0; // we have not yet read any elements of C

  for( int i = 0; i < totalCols; i++ ) {
    c[i] = 0.0;
  }

  while( DATALINE == (kindOfLine = GetLine( line ) ) ) {
    // we are still in the columns section.
    // We have already checked the syntax in this->scanColsSection.
    // Any (non-io) errors here are program errors, not syntax errors.
    ierr = this->ParseDataLine2( line, blank, colname, row[0], &val[0],
				hasSecondValue, row[1], &val[1] );
    assert( ierr != mpssyntaxerr );
    if( ierr != mpsok ) break;

    if (strcmp(oldColumnName, colname) != 0) {
      // we are not already working on this column
      strncpy( oldColumnName, colname, 16 );
      colnum = GetIndex( colTable, colname );
      assert( colnum >= 0 );
    }

    int nvals = (hasSecondValue) ? 2 : 1;
    int i;
    for( i = 0; i < nvals; i++ ) {
      // all rows specified
      int rownum = GetIndex( rowTable, row[i] );
      assert( rownum >= 0 );
      switch( rowInfo[rownum].kind ) {
	// on the kind of row
      case kFreeRow:
	if( 0 == strcmp( objectiveName, row[i] ) ) {
	  // T`his free row is the objective
	  c[colnum] = val[i];
	} // otherwise, skip it.
	break;
      case kEqualRow:
	this->insertElt( irowA, nnzA, jcolA, dA, nea,
			 rowRemap[rownum], colnum, val[i], ierr );
	assert( 0 == ierr );
	break;
      default:
	this->insertElt( irowC, nnzC, jcolC, dC, nec,
			 rowRemap[rownum], colnum, val[i], ierr );
	assert( 0 == ierr );
      } // end switch on the kind of row
    } // end for all rows specified
  } // end while we are still in the cols section.

  assert( nnzC == nec );
  assert( nnzA == nea );

  if( nnzA > 0 ) doubleLexSort( irowA, nnzA, jcolA, dA );
  if( nnzC > 0 ) doubleLexSort( irowC, nnzC, jcolC, dC );
}

void MpsReader::remapRows()
{
  // At this point,  we actually know which rows are equality 
  // constraints and which aren't. Remap them.
  rowRemap = new int[totalRows];
  int equalityRow = 0, inequalityRow = 0;
  int i;
  for( i = 0; i < totalRows; i++ ) {
    switch ( rowInfo[i].kind ) {
    case kFreeRow: break;
    case kEqualRow: 
      rowRemap[i] = equalityRow;
      equalityRow++;
      break;
    default:
      rowRemap[i] = inequalityRow;
      inequalityRow++;
      break;
    }
  }
}

void MpsReader::stuffMatrix( GenMatrix& A,
			     int irow[], int nnz, int jcol[], double dA[] )
{
  int info = 0;
  if( nnz > 0 ) A.putSparseTriple( irow, nnz, jcol, dA, info );
  assert( info == 0 ); // If not, there is a program error.
}

void MpsReader::stuffMatrix( SymMatrix& Q,
			     int irow[], int nnz, int jcol[], double dA[] )
{
  int info = 0;
  if( nnz > 0 ) Q.putSparseTriple( irow, nnz, jcol, dA, info );
  assert( info == 0 ); // If not, there is a program error.
}

void MpsReader::insertElt( int irow[], int len, int jcol[],
			   double dval[], int& ne,
			   int row, int col, double val,
			   int& ierr )
{
  if( ne >= len ) {
    ierr = mpssyntaxerr;
    return;
  } 
  ierr = mpsok;
  irow[ne] = row;
  jcol[ne] = col;
  dval[ne] = val;
  
  ne++;
}

void MpsReader::readRHSSection( double b[],
				double clow[], char iclow[],
				double cupp[], char icupp[],
				char line[], int& ierr, int& kindOfLine )
{
  char * seenRow = new char[totalRows];
  int i;
  for( i = 0; i < totalRows; i++ ) seenRow[i] = 0; // we haven't seen any yet

  if( my < 0 || mz < 0 ) {
    int nx_, my_, mz_; // Force the computation of the cached values
    this->getSizes( nx_, my_, mz_ );
  }
  for( i = 0; i < my; i++ ) {
    b[i] = 0;
  }
  for( i = 0; i < mz; i++ ) {
    clow[i] = 0; iclow[i] = 0;
    cupp[i] = 0; icupp[i] = 0;
  }
  objminus = 0.0;

  assert( rowRemap );

  for( i = 0; i < totalRows; i++ ) {
    switch( rowInfo[i].kind ) {
    case kFreeRow: case kEqualRow:
      break;
    case kLessRow:
      icupp[rowRemap[i]] = 1;
      break;
    case kGreaterRow:
      iclow[rowRemap[i]] = 1;
      break;
    default:
      iclow[rowRemap[i]] = 1;
      icupp[rowRemap[i]] = 1;
      break;
    }
  }

  char currentRHS[16] = "";
  char blank[4];
  char rhsName[16]="";
  char row[2][16];
  double val[2];
  int hasSecondValue;
  while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) {
    ierr = this->ParseDataLine2( line, blank, rhsName, row[0], &val[0],
				hasSecondValue, row[1], &val[1] ); 
    if( ierr != mpsok ) return;
    /* if( !isOnlySpaces( blank, 0, 1 ) ) {
      fprintf( stderr, "A code field is unexpected on line %d.\n", iline );
      ierr = mpssyntaxerr;
      return;
    } */
    if( 0 != strcmp( rhsName, currentRHS ) ) {
      if( 0 == strcmp( currentRHS, "" ) ) {
	strncpy( currentRHS, rhsName, 16 );
      } else { 
	fprintf( stderr, "Multiple rhs were specified.\n"
		 "The first rhs, \"%s\", will be used.\n",
		 currentRHS );
	// Skip the rest.
	while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) ;
	return;
      }
    }
      
    int nvals = (hasSecondValue) ? 2 : 1;
    for( i = 0; i < nvals; i++ ) {
      // all values specified
      int rownum = GetIndex( rowTable, row[i] );

      if( rownum < 0 ) {
	fprintf( stderr, "Unrecognized row name, \"%s\", on line %d.\n",
		 row[i], iline );
	ierr = mpssyntaxerr;
	return;
      }
      if( seenRow[rownum] ) {
	fprintf( stderr, "Multiple rhs were specified for row %s, "
		 "most recently at line %d.\n", rowInfo[i].name, iline );
	ierr = mpssyntaxerr;
	return;
      }
      seenRow[rownum] = 1;
      
      switch( rowInfo[rownum].kind ) {
      case kFreeRow:
	if( 0 == strcmp( objectiveName, rowInfo[rownum].name ) ) {
	  objminus = val[i];

        if( !strncmp( this->objectiveSense, "MAX", 3))
		objminus *= -1.0;
	}
	break;
      case kLessRow: 
      case kLessRowWithRange:
	cupp[ rowRemap[rownum] ] = val[i];
	break;
      case kGreaterRow:
      case kGreaterRowWithRange:
	clow[ rowRemap[rownum] ] = val[i];
	break;
      case kEqualRow:
	b[ rowRemap[rownum] ] = val[i];
	break;
      default: // Can't get here
	assert( 0 );
	break;
      }
    } // end for all values specified.
  }
  delete [] seenRow;
}


void MpsReader::readRHSSection( OoqpVector& b_,
				SimpleVector& clow, OoqpVector& iclow_,
				SimpleVector& cupp, OoqpVector& icupp_,
				char line[], int& ierr, int& kindOfLine )
{
  char   *iclow = 0, *icupp = 0;
  double *db    = 0, *dclow = 0, *dcupp = 0;
  if( my < 0 || mz < 0 ) {
    int nx_, my_, mz_; // Force the computation of the cached values
    this->getSizes( nx_, my_, mz_ );
  }
  SimpleVectorHandle b( new SimpleVector( my ) );
  if( my > 0 ) db = b->elements();
  if( mz > 0 ) {
    dclow = &clow[0];     dcupp = &cupp[0];
    iclow = new char[mz]; icupp = new char[mz];
  } 

  this->readRHSSection( db, dclow, iclow, dcupp, icupp,
			line, ierr, kindOfLine );

  b_.copyFrom( *b );

  iclow_.copyFromArray( iclow );
  icupp_.copyFromArray( icupp );

  delete [] iclow; delete [] icupp;
}

void MpsReader::readRangesSection( SimpleVector& clow, SimpleVector& cupp,
				   char line[], int& iErr, int& kindOfLine )
{
  double *dclow = 0, *dcupp = 0;
  
  if( clow.n > 0 ) dclow = &clow[0];
  if( cupp.n > 0 ) dcupp = &cupp[0];

  this->readRangesSection( dclow, dcupp, line, iErr, kindOfLine );
}

void MpsReader::readRangesSection( double clow[], double cupp[],
				   char line[], int& iErr, int& kindOfLine )
{
  char currentRange[10] = "";
  // The ranges section has already been scanned. Any syntax errors
  // left are programming errors
  char blank[4], rangeName[16]="", row[2][16];
  double val[2];
  int hasSecondValue;
  while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) {
    iErr = this->ParseDataLine2( line, blank, rangeName,
				row[0], &val[0],
				hasSecondValue,
				row[1], &val[1]);
    assert( iErr != mpssyntaxerr );
    if( iErr != mpsok ) return;

    if( 0 != strcmp( rangeName, currentRange) ) {
      // This is a new section of range values
      if( 0 == strcmp( currentRange, "" ) ) { 
	// This is the first range
	strncpy( currentRange, rangeName, 16 );
      } else {
	// This is the second range we have seen. We support only
	// one range, so skip the rest of the section
	while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) ;
	break;
      } // end else this is the second range we have seen
    } // end if this is a new section of range values
    int nvals = (hasSecondValue) ? 2 : 1;
    int i;
    for( i = 0; i < nvals; i++ ) {
      int rownum = GetIndex( rowTable, row[i] );
      assert( rownum >= 0 );
      int icrow  = rowRemap[rownum];
      switch( rowInfo[rownum].kind ) {
      case kGreaterRowWithRange: 
	cupp[icrow] = clow[icrow] + fabs( val[i] );
	break;

      case kLessRowWithRange:    
	clow[icrow] = cupp[icrow] - fabs( val[i] );
	break;
      default: // Can't get here.
	assert( 0 );
	break;
      }
    }
  }
}

void MpsReader::readBoundsSection( OoqpVector& xlow_, OoqpVector& ixlow_,
				   OoqpVector& xupp_, OoqpVector& ixupp_,
				   char line[], int& ierr, int& kindOfLine )
{
  // Force the computation of cached values
  if( my < 0 ) {
    int nx_, my_, mz_;
    this->getSizes( nx_, my_, mz_ );
  }
  // Create some temporary simple vectors
  SimpleVectorHandle xlow( new SimpleVector( totalCols ) );
  char *ixlow          = new char[ totalCols ];
  SimpleVectorHandle xupp( new SimpleVector( totalCols ) );
  char *ixupp          = new char[ totalCols ];

  this->readBoundsSection( xlow->elements(), ixlow, xupp->elements(), ixupp,
			   line, ierr, kindOfLine );

  xlow_ .copyFrom( *xlow );
  ixlow_.copyFromArray( ixlow );
  xupp_ .copyFrom( *xupp );
  ixupp_.copyFromArray( ixupp );

  delete [] ixlow; 
  delete [] ixupp;
}

void MpsReader::defaultBounds( OoqpVector& xlow_, OoqpVector& ixlow_,
			       OoqpVector& xupp_, OoqpVector& ixupp_ )
{
  // Force the computation of cached values
  if( my < 0 ) {
    int nx_, my_, mz_;
    this->getSizes( nx_, my_, mz_ );
  }
  // Create some temporary simple vectors
  SimpleVectorHandle xlow( new SimpleVector( totalCols ) );
  char *ixlow          = new char[ totalCols ];
  SimpleVectorHandle xupp( new SimpleVector( totalCols ) );
  char *ixupp          = new char[ totalCols ];

  this->defaultBounds( xlow->elements(), ixlow, xupp->elements(), ixupp );

  xlow_ .copyFrom( *xlow );
  ixlow_.copyFromArray( ixlow );
  xupp_ .copyFrom( *xupp );
  ixupp_.copyFromArray( ixupp );

  delete [] ixlow; 
  delete [] ixupp;
}

void MpsReader::defaultBounds( double xlow[], char ixlow[],
			      double xupp[], char ixupp[] )
{
  int i;
  for( i = 0; i < totalCols; i++ ) {
    xlow[i] = 0.0; ixlow[i] = 1; // Initially there is a lower bound
    xupp[i] = 0.0; ixupp[i] = 0; // but no upper bound
  }
}
void MpsReader::readBoundsSection( double xlow[], char ixlow[],
				   double xupp[], char ixupp[],
				   char line[], int& ierr, int& kindOfLine )
{
  int code;
  char bound[16], col[16];
  double val;

  char * lboundSpecified = new char[totalCols];
  char * uboundSpecified = new char[totalCols];

  int i;
  this->defaultBounds( xlow, ixlow, xupp, ixupp );
  for( i = 0; i < totalCols; i++ ) {
    lboundSpecified[i] = 0; // Nothing specified yet
    uboundSpecified[i] = 0;
  }

  while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) {
    ierr = this->ParseBoundsLine2( line, code, bound, col, &val );

    // we are reading datalines
    if( ierr != mpsok ) return;

    int colnum = GetIndex( colTable, col );
    if( colnum < 0 ) {
      fprintf( stderr, "Unrecognized column name on line %d.\n", iline );
      ierr = mpssyntaxerr;
      return;
    }
    int conflictingBound = 0;
    switch( code ) {
    case kLowerBound:
      if( lboundSpecified[colnum] ) {
	conflictingBound = 1;
      } else {
	ixlow[colnum] = 1; xlow[colnum] = val; lboundSpecified[colnum] = 1;
      }
      break;
    case kUpperBound:
      if( uboundSpecified[colnum] ) {
	conflictingBound = 1;
      } else {
	ixupp[colnum] = 1; xupp[colnum] = val; uboundSpecified[colnum] = 1;
      }
      break;
    case kFixedBound:
      if( lboundSpecified[colnum] || uboundSpecified[colnum] ) {
	conflictingBound = 1;
      } else {
	ixlow[colnum] = 1; xlow[colnum] = val; lboundSpecified[colnum] = 1;
	ixupp[colnum] = 1; xupp[colnum] = val; uboundSpecified[colnum] = 1;
      }
      break;
    case kFreeBound:
      if( uboundSpecified[colnum] || lboundSpecified[colnum] ) {
	conflictingBound = 1;
      } else {
	ixlow[colnum] = 0; xlow[colnum] = 0.0; lboundSpecified[colnum] = 1;
	ixupp[colnum] = 0; xupp[colnum] = 0.0; uboundSpecified[colnum] = 1;
      }
      break;
    case kMInftyBound:
      if( lboundSpecified[colnum] ) {
	conflictingBound = 1;
      } else {
	ixlow[colnum] = 0; xlow[colnum] = 0.0; lboundSpecified[colnum] = 1;
	
     // Commented out to ensure that is is consistent with how
     // kPInftyBound is treated.
    /* if( !uboundSpecified[colnum] ) {
	      xupp[colnum] = 0;
          ixupp[colnum] = 1; 
	    } */
      }
      break;
    case kPInftyBound:
      if( uboundSpecified[colnum] ) {
	conflictingBound = 1;
      } else {
	ixupp[colnum] = 0; xupp[colnum] = 0.0; uboundSpecified[colnum] = 1;
      }
      break;
    }
    if( conflictingBound ) {
      fprintf( stderr, "The bound specified on line %d for variable %s\n"
	       "may conflict with some eariler bound on the same variable.\n",
	       iline, col );
      ierr = mpssyntaxerr;
      return;
    }
  } // end while we are reading datalines

  for( i = 0; i < totalCols; i++ ) {
    if( ixlow[i] && ixupp[i] && xlow[i] > xupp[i] ) {
      fprintf( stderr, "The lower bound for variable \"%s\" is greater than\n"
	       "its upper bound.\n", colInfo[i].name );
      ierr = mpssyntaxerr;
      return;
    }
  }

  delete [] lboundSpecified;
  delete [] uboundSpecified;
}

void MpsReader::readHessSection( SymMatrix& Q,
				 char line[], int& ierr, int& kindOfLine )
{
  if( nnzA < 0 || nnzC < 0 ) {
    int nnzA_, nnzC_, nnzQ_; // Force the computation of these cached values
    this->numberOfNonZeros( nnzQ_, nnzA_, nnzC_ );
  }
  if( nnzQ == 0 ) {
    // The Hessian section must be empty. Skip out early
    kindOfLine = this->GetLine( line );
    return;
  }
  int * irowQ = new int[nnzQ];
  int * jcolQ = new int[nnzQ];
  double * dQ = new double[nnzQ];

  this->readHessSection( irowQ, jcolQ, dQ, line, ierr, kindOfLine );

  this->stuffMatrix( Q, irowQ, nnzQ, jcolQ, dQ );

  delete [] dQ;
  delete [] jcolQ;
  delete [] irowQ;
}

void MpsReader::readHessSection( int irowQ[], int jcolQ[], double dQ[],
				 char line[], int& ierr, int& kindOfLine )
{
  char colname[16];
  char name[2][16];
  char code[4];
  double val[2];
  int hasSecondValue;
  
  if( nnzA < 0 || nnzC < 0 ) {
    int nnzA_, nnzC_, nnzQ_; // Force the computation of these cached values
    this->numberOfNonZeros( nnzQ_, nnzA_, nnzC_ );
  }
  if( nnzQ == 0 ) {
    // The Hessian section must be empty. Skip out early
    kindOfLine = this->GetLine( line );
    return;
  }

  int neq     = 0; // nothing in there yet.
  char oldColName[16] = "";

  int colnum = -1;
  while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) {
    // we are still in the Hessian section
    // This section has already been scanned, so any syntax errors at
    // this point are program errors.
    ierr = this->ParseDataLine2(line, code, colname, name[0], &val[0],
			       hasSecondValue, name[1], &val[1] );
    assert( ierr != mpssyntaxerr );
    if( ierr != mpsok ) break;
    // are we already working on this column?
    if( 0 != strcmp( oldColName, colname) ) {
      // it is a new column
      colnum = GetIndex( colTable, colname );
      assert( colnum >= 0 );
    }
    int nvals = (hasSecondValue) ? 2 : 1;
    int i;
    for( i = 0; i < nvals; i++ ) {
      int rownum = GetIndex( colTable, name[i] );
      assert( rownum >= 0 );
      this->insertElt( irowQ, nnzQ, jcolQ, dQ, neq,
		       rownum, colnum, val[i], ierr );
      assert( ierr == 0 );
    }
  } // while we are still in the Hessian section
  assert( neq == nnzQ );
  doubleLexSort( irowQ, nnzQ, jcolQ, dQ );
}

void MpsReader::scanRangesSection( char line[200], 
                                   int& iErr, int& kindOfLine )
{
  // currentRange holds the name of the current range. There is no
  // current range at this point, so set it to empty.
  char currentRange[16] = "";
  char blank[4];
  char rangeName[16]="", row[2][16];
  int hasSecondValue;
  double val[2];
  int nvals;

  
  char * seenRow = new char[totalRows];
  if( !seenRow ) { iErr = mpsmemoryerr; return; };

  int i;
  for( i = 0; i < totalRows; i++ ) seenRow[i] = 0;
  iErr = mpsok;
  while( iErr == mpsok &&
	 DATALINE == (kindOfLine = this->GetLine( line ) ) ) {
    // we are reading data lines
    iErr = this->ParseDataLine2( line, blank, rangeName,
				row[0], &val[0],
				hasSecondValue,
				row[1], &val[1]);
    if( iErr != mpsok ) break;
    if( 0 != strcmp( rangeName, currentRange) ) {
      // This is a new section of range values
      if( 0 == strcmp( currentRange, "" ) ) { 
	// This is the first range
	strncpy( currentRange, rangeName, 16 );
      } else {
	// This is the second range we have seen. We support only
	// one range, so skip the rest of the section
	fprintf( stderr, "Warning, we only support one set of "
		 "values in the RANGES section.\n" );
	while( DATALINE == (kindOfLine = this->GetLine( line ) ) ) ;
	break;
      } // end else this is the second range we have seen
    } // end if this is a new section of range values
    nvals = (hasSecondValue) ? 2 : 1;
    for( i = 0; i < nvals; i++ ) {
      // all rows specified
      int rownum = GetIndex( rowTable, row[i] );
      if( rownum < 0 ) {
	fprintf( stderr, 
		 "Unrecognized row name %s at line %d.\n",
		 row[i], iline );
	iErr = mpssyntaxerr; break;
      }
      if( seenRow[rownum] ) {
	fprintf( stderr, "The range for row %d has been specified twice, "
		 "most recently at line %d.\n", rownum, iline );
	iErr = mpssyntaxerr;
	break;
      }
      seenRow[rownum]++;
      this->rowHasRange( rownum, val[i], iErr );
      if( iErr != mpsok ) break;
    } // end for all rows specified.
  } // end while we are reading data lines
    
  delete [] seenRow;
}

void MpsReader::rowHasRange( int rownum, double val, int& iErr )
{
  iErr = mpsok;

  if( val == 0 ) {
    fprintf( stderr, "A zero value has been specified "
	     "for a range on line %d.\n", iline );
    iErr = mpssyntaxerr;
    return;
  }

  int kind = rowInfo[rownum].kind;
  switch( kind ) {
  case kLessRow:
    kind = kLessRowWithRange;
    break;
  case kGreaterRow:
    kind = kGreaterRowWithRange;
    break;
  case kEqualRow:
    if( val > 0 ) {
      kind = kGreaterRowWithRange;
    } else {
      kind = kLessRowWithRange;
    }
    break;
  case kFreeRow:
    fprintf( stderr, "A range has been specified for a free "
	     "row on line %d.\n", iline );
    iErr = mpssyntaxerr;
    return;
    break;
  default:
    assert( 0 && "Can't get here." );
    break;
  }
  rowInfo[rownum].kind = kind;
}

void MpsReader::expectHeader( int lineType, const char expectName[],
			     char line[], int& ierr )
{
  ierr = 0;
  if( lineType == HEADERLINE ) {
    char name[16];
    ierr = this->ParseHeaderLine( line, name ); 
    if( ierr == mpsok ) {
      if( !isLJustOf( name, expectName, 16 ) ) {
	ierr = mpssyntaxerr;
	fprintf( stderr, "Expected %s at line %d, got %s.\n",
		 expectName, iline, name );
      }
    }
  } else {
    ierr = mpssyntaxerr;
    fprintf( stderr, "Expected %s at line %d.\n", expectName, iline );
  }
}

void MpsReader::expectHeader2( int lineType, const char expectName[],
			     char line[], int& ierr )
{
  ierr = 0;
  if( lineType == HEADERLINE ) {
    char name[16];
    ierr = this->ParseHeaderLine2( line, name ); 
    if( ierr != mpsok ) {
        ierr = mpssyntaxerr;
        fprintf( stderr, "Expected %s at line %d.\n", expectName, iline );
        }
  }
}

int MpsReader::acceptHeader( int lineType, const char acceptName[],
			     char line[], int& ierr )
{
  if( lineType == HEADERLINE ) {
    char name[16];
    ierr = this->ParseHeaderLine( line, name );
    if( ierr != mpsok ) return 0;
    if( isLJustOf( name, acceptName, 16 ) ) {
      return 1;
    } else {
      return 0;
    }
  } else {
    fprintf( stderr, "Expected a new section to start at line %d.\n", iline );
    ierr = mpssyntaxerr;
    return 0;
  }
}

int MpsReader::acceptHeader2( int lineType, const char acceptName[],
			     char line[], int& ierr )
    {
    if( lineType == HEADERLINE ){
        char name[16];
        ierr = this->ParseHeaderLine2( line, name );

        if( ierr == mpsok && !strcmp( name, acceptName))
            return 1;
        else
            return 0;
        }
    else{
        fprintf( stderr, "Expected a new section to start at line %d.\n", iline );
        ierr = mpssyntaxerr;
        return 0;
        }
    }

void MpsReader::scanHessSection( char line[200], 
                                 int& iErr, int& linetype )
{
  char code[4], name[2][16], colname[16];
  double val[2];
  int hasSecondValue;
  char oldColName[16] = "";
  int colnum = -1;

  int  * lastSeenRow = 0;
  linetype = mpssyntaxerr; // If this doesn't get set to something else,
  // it is an error
  
  int i, nvals;
  iErr = mpsok;

  lastSeenRow  = new int[totalCols];
  if( !lastSeenRow ) {
    iErr = mpsmemoryerr;
  } else {
    for( i = 0; i < totalCols; i++ ) {
      lastSeenRow[i] = -1;
    }
    while((linetype = this->GetLine(line)) == DATALINE ) {
      // we are reading data lines
      iErr = this->ParseDataLine2(line, code, colname, name[0], &val[0],
                                 hasSecondValue, name[1], &val[1] );
      if( iErr != mpsok ) break;
      // Are we already working on this column
      if( 0 != strcmp( oldColName, colname) ) {
	// it is a new column
	int lastcolnum = colnum;
	colnum = GetIndex( colTable, colname );
	if( colnum < 0 ) { 
	  fprintf( stderr, "Unrecognized column name %s in line %d.\n",
		   colname, iline ); 
	  iErr = mpssyntaxerr; break;
	}
	if( colnum <= lastcolnum ) {
          fprintf( stderr, "Column out of order at line %d.\n", iline );
	  iErr = mpssyntaxerr; break;
        }
        // Mark this column as seen
        strncpy( oldColName, colname, 16 );
      }
      nvals = (hasSecondValue) ? 2 : 1;
      for( i = 0; i < nvals; i++ ) {
        int rownum = GetIndex( colTable, name[i] );
        if( rownum < 0 ) {
          fprintf( stderr, 
                   "Unrecognized variable name %s at line %d.\n",
                   name[i], iline );
          iErr = mpssyntaxerr; break;
        }
        if( lastSeenRow[rownum] == colnum ) {
          fprintf( stderr,
                   "Element (%s, %s) specified twice.\n", 
                   name[i], colname );
          iErr = mpssyntaxerr; break;
        }
        if( colnum > rownum ) {
          fprintf( stderr, 
                   "Error line %d: element (%s, %s)\nis not in "
                   "the lower triangle of QUADOBJ.\n",
		   iline, name[i], colname );
          iErr = mpssyntaxerr; break;
        }
        // All is well, record seeing this element in *rownum*
        colInfo[rownum].nnz++;
        lastSeenRow[rownum] = colnum;
      }
      if( iErr != mpsok ) break;
    } // end while we are reading data lines
  }
  delete [] lastSeenRow;

  if( iErr == mpsok ) {
    // There isn't already an error.
    switch( linetype ) {
    case HEADERLINE:   iErr = mpsok;         break;
    case mpssyntaxerr: iErr = mpssyntaxerr;  break;
    case mpsioerr:     iErr = mpsioerr;      break;
    default:           iErr = mpsunknownerr; break;
    }
  }
}

MpsReader * MpsReader::newReadingFile( char filename[], int& iErr )
{ 
  MpsReader * reader;

  iErr                  = mpsunknownerr;
  FILE * file;          char * resolvedName;

  MpsReader::findFile( file, resolvedName, filename );
  if( !file ) { 
    iErr = mpsfileopenerr; 
    return 0; 
  } else {
    reader = new MpsReader( file );
    reader->infilename = resolvedName;
  }

  // now get on with the reading
  reader->scanFile( iErr );

  if ( 0 == iErr ) {
    // looks OK, return pointer to the MpsReader object.
    // Force the computation of the non-zeros
    int dummy1, dummy2, dummy3;

    reader->numberOfNonZeros( dummy1, dummy2, dummy3 );
    return reader;
  } else {
    // there's been a problem. release the space and return a null.
    int closeErr;
    reader->releaseFile( closeErr );
    delete reader;
    return 0;
  }
}
  
void MpsReader::findFile( FILE*& file, char*& resolvedName,
			  char filename[] )
{
  file = NULL; resolvedName = 0;
  int lfilename = strlen( filename );

  file = fopen( filename, "r" );

  if( file ) {
    resolvedName = new char[lfilename + 1];
    strcpy( resolvedName, filename );
    return;
  }
  // Otherwise try with an .mps suffix

  if( lfilename < 4 || 0 != strcmp( &filename[lfilename-4], ".mps") ) {
    // the file doesn't already have an .mps suffix
    // append one
    resolvedName = new char[lfilename + 4];
    strcat(resolvedName, ".mps");
    file = fopen( resolvedName, "r" );
    if( !file ) {
      delete [] resolvedName; resolvedName = 0;
    }
  } // end if we didn't find the file.
}

void MpsReader::readProblemName( char line[], int& iErr, int kindOfLine )
{
  int extra_crud;
  char tag[6];

  if( HEADERLINE == kindOfLine ) {
    this->string_copy( tag, &line[0], 4 );  /* characters  1 - 4 to name */
    if( !isLJustOf( "NAME", tag, 4 ) ) {
      fprintf( stderr, "Expected NAME on line %d, got %s.\n",
	       iline, tag );
      iErr = mpssyntaxerr;
      return;
    }
    extra_crud = !isOnlySpaces( line, 4, 13 );
    if( !extra_crud ) {
      this->string_copy( problemName, &line[14], 8 );
    }
    extra_crud = extra_crud || !isOnlySpaces( line, 22, 60 );
    if( extra_crud ) {
      fprintf( stderr, "Extra characters in NAME field on line %d.\n",
	       iline );
      fprintf( stderr, "These will be ignored. Only the first 8 characters are significant: '%s'.\n", problemName);
      // iErr = mpssyntaxerr;
      // return;
    }
  } else {
    fprintf( stderr, "Expected NAME on line %d.\n", iline );
    iErr = mpssyntaxerr;
    return;
  }
  iErr = mpsok;
  return;
}

void MpsReader::readProblemName2( char line[], int& iErr, int kindOfLine )
{
    int i = 0;
    char *token;
    char *arrayOfTokens[2];
    int ierr =0;
    char tempLine[200];
    char tag[6];

    if( HEADERLINE == kindOfLine ) {
        strncpy( tempLine, line, 200);
        token = strtok( tempLine, " ");

        // Split the extracted line into tokens delimited by space...
        arrayOfTokens[0] = token;
        arrayOfTokens[1] = strtok( NULL, " ");

        // Field 1: Tag "Name"
        if( arrayOfTokens[0] != NULL){
            strcpy(tag, arrayOfTokens[0]);
            
            if( strncmp( "NAME", tag, 4 ) ) {
                fprintf( stderr, "Expected NAME on line %d, got %s.\n",
                   iline, tag );
                iErr = mpssyntaxerr;
                return;
                }
            }
        else{
            fprintf( stderr, "Empty tag NAME on line %d.\n", iline );
            iErr = mpssyntaxerr;
            return;
            }

        // Field 2: Problem Name
        if( arrayOfTokens[1] != NULL){
            this->string_copy(problemName, arrayOfTokens[1], 16);
            
            if( strlen( arrayOfTokens[1]) > 16){
                fprintf( stderr, "Extra characters in NAME field on line %d.\n",
                   iline );
                fprintf( stderr, "These will be ignored. Only the first 16 characters are significant: '%s'.\n", problemName);
                }
            }
        else{
            fprintf( stderr, "Empty problem name on line %d.\n", iline );
            iErr = mpssyntaxerr;
            return;
            }
        }
    else {
        fprintf( stderr, "Expected tag NAME on line %d.\n", iline );
        iErr = mpssyntaxerr;
        return;
        }

  iErr = mpsok;
  return;
}

void MpsReader::readObjectiveSense( char line[], int& iErr, int kindOfLine )
{
    int i = 0;
    char *token;
    char *arrayOfTokens[1];
    int ierr =0;
    char tempLine[200];

    if( DATALINE == kindOfLine ) {
        strncpy( tempLine, line, 200);
        token = strtok( tempLine, " ");

        // Split the extracted line into tokens delimited by space...
        arrayOfTokens[0] = token;

        // Field 1: "MAX" or "MIN"
        if( arrayOfTokens[0] != NULL){
            strncpy( objectiveSense, arrayOfTokens[0], 3);

            if( strncmp( "MAX", objectiveSense, 3 ) && strncmp( "MIN", objectiveSense, 3) ) {
                fprintf( stderr, "Expected objective sense MAX or MIN on line %d, got %s.\n",
                   iline, objectiveSense );
                iErr = mpssyntaxerr;
                return;
                }
            }
        else{
            fprintf( stderr, "Empty objective sense on line %d.\n", iline );
            iErr = mpssyntaxerr;
            return;
            }
        }
    else {
        fprintf( stderr, "Expected objective sense MAX or MIN on line %d.\n", iline );
        iErr = mpssyntaxerr;
        return;
        }

  iErr = mpsok;
  return;
}

void MpsReader::scanFile( int& iErr )
{
  char line[200];
  int kindOfLine;

  // Problem name
  kindOfLine = this->GetLine(line);
  this->readProblemName2( line, iErr, kindOfLine );
  if( iErr != mpsok ) return;
  kindOfLine = this->GetLine(line);

  // Objective sense - MAX or MIN (optional)
  if( this->acceptHeader2( kindOfLine, "OBJSENSE", line, iErr ) ) {
        kindOfLine = this->GetLine(line);
        this->readObjectiveSense( line, iErr, kindOfLine );

        if( iErr != mpsok ) return;
        kindOfLine = this->GetLine(line);
    }

  // ROWS section
  this->expectHeader2( kindOfLine, "ROWS", line, iErr );
  if( iErr != mpsok ) return;
  this->readRowsSection( line, iErr,  kindOfLine );
  if( iErr != mpsok ) return;

  // COLUMNS section
  this->expectHeader2( kindOfLine, "COLUMNS", line, iErr );
  if( iErr != mpsok ) return;
  this->scanColsSection( line, iErr,  kindOfLine );
  if( iErr != mpsok ) return;
  
  // RHS section - required, but we skip it during the scan phase
  this->expectHeader2( kindOfLine, "RHS", line, iErr );
  if( iErr != mpsok ) return;
  while((kindOfLine = this->GetLine(line)) == DATALINE) ;
  
  // RANGES (optional)
  if( this->acceptHeader2( kindOfLine, "RANGES", line, iErr ) ) {
    this->scanRangesSection( line, iErr, kindOfLine );
  }
  if( iErr != mpsok ) return;

  this->remapRows();

  // BOUNDS (optional)
  if( this->acceptHeader2( kindOfLine, "BOUNDS", line, iErr ) ) {
    // skip it
    while((kindOfLine = this->GetLine(line)) == DATALINE) ;
  }
  if( iErr != mpsok ) return;

  // QUADOBJ (optional)
  if( this->acceptHeader2( kindOfLine, "QUADOBJ", line, iErr ) ) {
    this->scanHessSection( line, iErr, kindOfLine );
  }
  //szhu - extend to QMATRIX
  else if( this->acceptHeader2( kindOfLine, "QMATRIX", line, iErr ) ) {
    this->scanHessSection( line, iErr, kindOfLine );
  }
  if( iErr != mpsok ) return;

  this->expectHeader2( kindOfLine, "ENDATA", line, iErr );
}


void MpsReader::readRowsSection( char line[200], 
				 int& iErr, int& linetype )
{
  char rname[16], code[4];
  const int rowsGuess          = 1000;
  const double rowsBlockFactor = 1.5;

  int lrowInfo   = rowsGuess;
  totalRows      = 0;
  rowInfo        = 0;
  rowTable       = 0;

  // Remember code for rhs and bound
  rowInfo = new MpsRowInfo[rowsGuess];
  if( rowInfo == 0 ) {
    fprintf( stderr, "Out of memory\n" );
    iErr = mpsmemoryerr; return;
  }

  int foundObjective = 0;

  while ((linetype = this->GetLine(line)) == DATALINE) {
    // we are reading datalines
    iErr = this->ParseRowsLine2(line, code, rname );
    if( iErr != mpsok ) break;

    int rowType = MpsRowTypeFromCode2( code[0] );
    if( rowType == kBadRowType ) {
      fprintf( stderr, "Unrecognized row type\n"); 
      iErr = mpssyntaxerr;
      break;
    }
    // Insert the row
    if( rowType == kFreeRow ) { // This may be the objective
      if( foundObjective ) {
	if( foundObjective < 2 ) {
	  fprintf( stderr, 
		   "Warning: More than one objective function was specified.\n"
		   "The first one found, \"%s\", will be used.\n",
		   objectiveName );
	}
      } else {
	this->string_copy( objectiveName, rname, 16 );
      }
      foundObjective++;
    }

    // Reallocate if necessary
    if( totalRows >= lrowInfo ) {
      int lNewRowInfo = (int) (lrowInfo * rowsBlockFactor);
      MpsRowInfo * newRowInfo;

      try {
	newRowInfo = new MpsRowInfo[lNewRowInfo];
      } catch(...) {
	cerr << "Out of memory in MpsReader.\n";
	iErr = mpsmemoryerr;
        break;
      }
      for( int k = 0; k < lrowInfo; k++ ) {
	newRowInfo[k] = rowInfo[k];
      }
      delete [] rowInfo;
      rowInfo = newRowInfo;
      lrowInfo = lNewRowInfo;
    }
    this->string_copy( rowInfo[totalRows].name, rname, 16 );
    rowInfo[totalRows].kind    = rowType;
    rowInfo[totalRows].nnz     = 0;
    
    
    totalRows++;
  } // end while we are reading data lines

  if( iErr == 0 ) {
//      rowInfo = 
//        (MpsRowInfo *) realloc( rowInfo, totalRows * sizeof( MpsRowInfo ) );
    // This is smaller, so it shouldn't fail.
    rowTable = NewHashTable( 2 * totalRows );
    if ( rowTable ) {
      int i;
      for( i = 0; i < totalRows; i++ ) {
        if (Insert( rowTable, rowInfo[i].name, i) == 1) {
          fprintf( stderr,
                   "The row name %s was used twice.\n", rowInfo[i].name );
          iErr = mpssyntaxerr;
          break;
        }
      }
    } else {
      iErr = mpsmemoryerr;
    }
  }
  if( iErr == 0 ) {
    // There isn't already an error.
    switch( linetype ) {
    case HEADERLINE:   iErr = mpsok;         break;
    case mpssyntaxerr: iErr = mpssyntaxerr;  break;
    case mpsioerr:     iErr = mpsioerr;      break;
    default:           iErr = mpsunknownerr; break;
    }
  }
  if( iErr != 0 ) {
    delete [] rowInfo;
    if( rowTable ) DeleteHashTable( rowTable );
  }
}

void MpsReader::scanColsSection( char line[200], 
				 int& iErr, int& linetype )
{
  char colname[16], code[4], name[2][16];
  int hasSecondValue;
  double val[2];
  int * lastSeen;

  int colsGuess = 1000;
  const double colsBlockFactor = 1.5;
  char oldColumnName[17] = "";
  int nvals;
  int colnum = -1;
  iErr       = mpsok;

  int i;

  // record the current position
  firstColumnLine = iline;
  columnFilePosition = ftell( file );
  if( -1 == columnFilePosition ) {
    fprintf( stderr, "Could not record the current file position.\n" );
    iErr = mpsioerr;
    return;
  }

  // allocate some space to the colnames array (make an initial guess
  // of the size)

  colTable  = 0; totalCols = 0; colnum = 0;
  lastSeen  = new int[totalRows];

  if(2*totalRows >  colsGuess) colsGuess = 2*totalRows;
  int lcolnames = colsGuess;
  colInfo   = new MpsColInfo[colsGuess];

  if( !colInfo || !lastSeen ) {
    iErr = mpsmemoryerr;
  } else { // memory was allocated sucessfully
    // Nothing has yet been seen
    for( i = 0; i < totalRows; i++ ) {
      lastSeen[i] = -1;
    }

    while ((linetype = this->GetLine(line)) == DATALINE) {
      // we are reading datalines
      iErr = this->ParseDataLine2(line, code, colname, name[0], &val[0],
                                 hasSecondValue, name[1], &val[1]);
      if( iErr != mpsok ) break;
      // are we already working on this column?
      if (strcmp(oldColumnName, colname) != 0) {
        /* it's a new column */
        colnum = totalCols++;
        
        this->string_copy( oldColumnName, colname, 16);
        /* register its name */
        this->string_copy( colInfo[colnum].name, colname, 16 );
        colInfo[colnum].nnz  = 0;
        
        if( totalCols>= lcolnames ) {
          // We must reallocate
	  int lNewColInfo = (int) (lcolnames * colsBlockFactor) ;
	  MpsColInfo * newColInfo;
	  try {
	    newColInfo = new MpsColInfo[lNewColInfo];
	  } catch( ... ) {
	    iErr = mpsmemoryerr;
	    break;
	  }
	  for( int k = 0; k < lcolnames; k++ ) {
	    newColInfo[k] = colInfo[k];
	  }
	  delete [] colInfo;

          colInfo = newColInfo;
	  lcolnames = lNewColInfo;
        } // end if we must reallocate
      }  // end if it's a new column
      
      nvals = hasSecondValue ? 2 : 1;
      for( i = 0; i < nvals; i++ ) {
        // What row is the value in?
        int rownum = GetIndex( rowTable, name[i] );
        if( rownum < 0 ) {
          iErr = mpssyntaxerr;
          fprintf( stderr, "Unrecognized row name" );
          break;
        }
        if( lastSeen[rownum] == colnum ) {
          iErr = mpssyntaxerr;
          fprintf( stderr,
                   "Error on line %d: "
                   "row %s was already specified for column %s.\n",
                   iline, name[i], colname );
          break;
        }
        rowInfo[rownum].nnz++;
        lastSeen[rownum] = colnum;
      }
      if( iErr != mpsok ) break;
    } // end while we are reading datalines.
  }

  delete [] lastSeen;

  if( iErr == mpsok ) {
    // No error yet.
    // Shrink the colInfo to actual size
//      colInfo = 
//        (MpsColInfo *) realloc( colInfo, totalCols * sizeof( MpsColInfo ) );
    
    colTable = NewHashTable( 2 * totalCols );
    if( !colTable ) { 
      iErr = mpsmemoryerr;
    } else { // The column table was successfully allocated.
      for( i = 0; i < totalCols; i++ ) { // loop over all columns
        if (Insert( colTable, colInfo[i].name, i) == 1) {
          // This column name was already found in the hash table.
          fprintf( stderr, "Column %s was specified twice.\n",
                   colInfo[i].name );
          iErr = mpssyntaxerr; break;
        }
      } // end loop over all columns
    } // end else the column table was allocated
  } // end if no error yet
  if( iErr == mpsok ) {
    // There isn't already an error.
    switch( linetype ) {
    case HEADERLINE:   iErr = mpsok;         break;
    case mpssyntaxerr: iErr = mpssyntaxerr;  break;
    case mpsioerr:     iErr = mpsioerr;      break;
    default:           iErr = mpsunknownerr; break;
    }
  }
  if( iErr != 0 ) {
    delete [] rowInfo;
    if( rowTable) { DeleteHashTable( rowTable ); rowTable = 0; }
  }
}


int MpsReader::string_copy( char dest[], char str[], int max)
{
  /* Copy string to dest converting non-printable characters to spaces.
   * Null-terminate dest at max+1 */
  
  strncpy( dest, str, max );
  dest[max] = '\0';

  return 0;
}

void MpsReader::getSizes( int& nx_, int& my_, int& mz_ )
{
  if( my < 0 || mz < 0 ) {
    my = 0; mz = 0;
    int i;
    for( i = 0; i < totalRows; i++ ) {
    if( rowInfo[i].kind == kEqualRow ) {
      my++;
    } else if ( rowInfo[i].kind == kFreeRow ) {
      // Do nothing
    } else {
      mz++;
    }
    }
  }
  nx_ = totalCols;
  my_ = my;
  mz_ = mz;
}

void MpsReader::numbersOfNonZeros( int lnnzQ[], int lnnzA[], int lnnzC[] )
{
    int i;
    for( i = 0; i < totalRows; i++ ) {
      if( rowInfo[i].kind == kEqualRow ) {
	lnnzA[ rowRemap[i] ] = rowInfo[i].nnz;
      } else if ( rowInfo[i].kind == kFreeRow ) {
	// Do nothing
      } else {
	lnnzC[ rowRemap[i] ] = rowInfo[i].nnz;
      }
    }
    for( i = 0; i < totalCols; i++ ) {
      lnnzQ[i] = colInfo[i].nnz;
    }
}

void MpsReader::numberOfNonZeros( int& nnzQ_, int& nnzA_, int& nnzC_ )
{
  if( nnzQ < 0 || nnzA < 0 || nnzC < 0 ) {
    nnzQ = 0; nnzA = 0; nnzC = 0;
    
    int i;
    for( i = 0; i < totalRows; i++ ) {
      if( rowInfo[i].kind == kEqualRow ) {
	nnzA += rowInfo[i].nnz;
      } else if ( rowInfo[i].kind == kFreeRow ) {
	// Do nothing
      } else {
	nnzC += rowInfo[i].nnz;
      }
    }
    for( i = 0; i < totalCols; i++ ) {
      nnzQ += colInfo[i].nnz;
    }
  }

  nnzA_ = nnzA;
  nnzC_ = nnzC;
  nnzQ_ = nnzQ;
}

void MpsReader::readQpGen( double     c[],  
			   int    irowQ[],  int  jcolQ[],  double dQ[],
			   double  xlow[],  char ixlow[],
			   double  xupp[],  char ixupp[],
			   int    irowA[],  int  jcolA[],  double dA[],
			   double     b[],
			   int    irowC[],  int  jcolC[],  double dC[],
			   double  clow[],  char iclow[],
			   double  cupp[],  char icupp[],
			   int&    ierr )
{
  char line[200];
  int  kindOfLine;
  this->readColsSection( c,
			 irowA, jcolA, dA,
			 irowC, jcolC, dC,
			 line, ierr, kindOfLine );
  if( ierr != mpsok ) return;

  // RHS section - required
  this->expectHeader2( kindOfLine, "RHS", line, ierr );
  if( ierr != mpsok ) return;
  this->readRHSSection( b, clow, iclow, cupp, icupp,
			line, ierr, kindOfLine );

  // RANGES (optional)
  if( this->acceptHeader2( kindOfLine, "RANGES", line, ierr ) ) {
    this->readRangesSection( clow, cupp,
			     line, ierr, kindOfLine );
  }
  if( ierr != mpsok ) return;

  // BOUNDS (optional)
  if( this->acceptHeader2( kindOfLine, "BOUNDS", line, ierr ) ) {
    this->readBoundsSection( xlow, ixlow, xupp, ixupp,
			     line, ierr, kindOfLine );
  } else {
    this->defaultBounds( xlow, ixlow, xupp, ixupp );
  }
  if( ierr != mpsok ) return;
  // QUADOBJ (optional)
  if( this->acceptHeader2( kindOfLine, "QUADOBJ", line, ierr ) ) {
    this->readHessSection( irowQ, jcolQ, dQ, line, ierr, kindOfLine );
  }
  //szhu - extend to QMATRIX
  else if( this->acceptHeader2( kindOfLine, "QMATRIX", line, ierr ) ) {
    this->readHessSection( irowQ, jcolQ, dQ, line, ierr, kindOfLine );
  }
  if( ierr != mpsok ) return;

  this->expectHeader2( kindOfLine, "ENDATA", line, ierr );

}

void MpsReader::readQpBound( OoqpVector& c, SymMatrix& Q,
			     OoqpVector& xlow, OoqpVector& ixlow,
			     OoqpVector& xupp, OoqpVector& ixupp,
			     int& iErr )
{
  if( my > 0 || mz > 0 ) { iErr = 1024; return; }
    char line[200];
  int  kindOfLine;

  {
    SimpleVectorHandle sc( new SimpleVector(totalCols) );

    this->readColsSection( sc->elements(),
			   0, 0, 0, // elements of A
			   0, 0, 0, // elements of C
			   line, iErr,  kindOfLine );
    c.copyFrom( *sc );
  }
  if( iErr != mpsok ) return;
  
  // RHS section - required
  this->expectHeader2( kindOfLine, "RHS", line, iErr );
  if( iErr != mpsok ) return;
  kindOfLine = this->GetLine( line );
  assert(HEADERLINE == kindOfLine );
  // BOUNDS (optional)
  if( this->acceptHeader2( kindOfLine, "BOUNDS", line, iErr ) ) {
    this->readBoundsSection( xlow, ixlow, xupp, ixupp,
			     line, iErr, kindOfLine );
  } else {
    this->defaultBounds( xlow, ixlow, xupp, ixupp );
  }
  if( iErr != mpsok ) return;
  
  // QUADOBJ (optional)
  if( this->acceptHeader2( kindOfLine, "QUADOBJ", line, iErr ) ) {
    this->readHessSection( Q, line, iErr, kindOfLine );
  }
  //szhu - extend to QMATRIX
  else if( this->acceptHeader2( kindOfLine, "QMATRIX", line, iErr ) ) {
    this->readHessSection( Q, line, iErr, kindOfLine );
  }
  if( iErr != mpsok ) return;
  
  this->expectHeader2( kindOfLine, "ENDATA", line, iErr );
}

void MpsReader::readQpGen( OoqpVector& c, SymMatrix& Q,
                           OoqpVector& xlow, OoqpVector& ixlow,
			   OoqpVector& xupp, OoqpVector& ixupp,
                           GenMatrix& A, OoqpVector& b,
                           GenMatrix& C,
			   OoqpVector& clow_, OoqpVector& iclow,
			   OoqpVector& cupp_, OoqpVector& icupp,
                           int& iErr )
{
  char line[200];
  int  kindOfLine;

  this->readColsSection( c, A, C, line, iErr,  kindOfLine );
  if( iErr != mpsok ) return;
  
  { 
    SimpleVectorHandle clow( new SimpleVector( mz ) );
    SimpleVectorHandle cupp( new SimpleVector( mz ) );

    // RHS section - required
    this->expectHeader2( kindOfLine, "RHS", line, iErr );
    if( iErr != mpsok ) return;
    this->readRHSSection( b, *clow, iclow, *cupp, icupp,
			  line, iErr, kindOfLine );
    if( iErr != mpsok ) return;
    
    // RANGES (optional)
    if( this->acceptHeader2( kindOfLine, "RANGES", line, iErr ) ) {
      this->readRangesSection( *clow, *cupp,
			       line, iErr, kindOfLine );
    }
    clow_.copyFrom( *clow );
    cupp_.copyFrom( *cupp );
  }
  if( iErr != mpsok ) return;

  // BOUNDS (optional)
  if( this->acceptHeader2( kindOfLine, "BOUNDS", line, iErr ) ) {
    this->readBoundsSection( xlow, ixlow, xupp, ixupp,
			     line, iErr, kindOfLine );
  } else {
    this->defaultBounds( xlow, ixlow, xupp, ixupp );
  }
  if( iErr != mpsok ) return;

  // QUADOBJ (optional)
  if( this->acceptHeader2( kindOfLine, "QUADOBJ", line, iErr ) ) {
    this->readHessSection( Q, line, iErr, kindOfLine );
  }
  //szhu - extend to QMATRIX
  else if( this->acceptHeader2( kindOfLine, "QMATRIX", line, iErr ) ) {
    this->readHessSection( Q, line, iErr, kindOfLine );
  }
  if( iErr != mpsok ) return;

  this->expectHeader2( kindOfLine, "ENDATA", line, iErr );

}

void MpsReader::releaseFile( int& ierr )
{
  if ( file == stdin || // Don't close stdin.
       0 == fclose( file ) ) {
    file = 0;
    ierr = 0;
  } else {
    ierr = errno;
  }
}


MpsReader::~MpsReader()
{
  // The user should have already closed the file 
  // (we don't want to close it in the destructor, because we don't want
  // to throw and error from the destructor if the file doesn't close.)
  assert( file == 0 && "You forgot to call MpsReader::releaseFile");

  if( colTable ) DeleteHashTable( colTable );
  if( rowTable ) DeleteHashTable( rowTable );
  
  delete [] rowInfo;
  delete [] colInfo;

  delete [] infilename;
  delete [] rowRemap;
}

// Universal end-of-line
inline int isEndOfLine( char c, FILE * file )
{
  if( '\n' == c ) return 1;
  if( '\r' == c ) {
    int newChar = getc( file );
    if( newChar != '\n' ) {
      ungetc( newChar, file );
    }
    return 1;
  }
  return 0;
}

int MpsReader::GetLine_old(char * line )
{
  int             i, j, terminated;
  const int length = 61;

  do {
    int c;

    iline++;
 
    i = 0; terminated = 0;
    while( i < length && EOF != ( c = getc( file ) ) ) {
      if( isEndOfLine( c, file ) ) {
	terminated = 1;
	break;
      }
      if( !isprint( c ) ) c = ' '; // Eliminate non-printing characters
      line[i] = c;
      
      i++;
    }
    if( !terminated ) {
      // The line wasn't terminated, or was possibly just terminated
      // after the 61st character.
      if( i == 0 ) {
	// Nothing was read. This is an error.
        fprintf( stderr, "Unexpected end-of-file at line %d.\n", iline );
        return READERROR;
      } else if( i == length ) { 
	// We read 61 characters without an error or a '\n'
	if( line[0] == '*' ) { 
	  // This is a comment line. There is no restriction on length.
	  // Just throw away characters til the end of line
          while( EOF != ( c = getc( file ) ) && c != '\n' ) ;
	} else {
	  // This is not a comment line. The only characters
	  // beyond the 61st column should be ' '
	  int extracrud = 0;
	  while( EOF != ( c = getc( file ) ) ) {
	    if( isEndOfLine( c, file ) ) break;
	    if( c != ' ' ) extracrud = 1;
	  }
	  if(  extracrud ) {
	    fprintf( stderr,
		     "Extra characters beyond column 61 in line %d.\n",
		     iline );
	    return mpssyntaxerr;
	  }
	} // end else this is not a comment line
      } // end else we read 61 characters
      // Otherwise, the file just doesn't have a terminating newline,
      // which is acceptable.
    }
    
    // Fill in the rest of line with spaces
    for( j = i; j < length; j++ ) line[j] = ' ';
    // Terminate the line
    line[length] = '\0';
  } while ( line[0] == '*' ); // Disgard comment lines
  
  if(line[0] == ' ') return DATALINE;
  else return HEADERLINE;
}

int MpsReader::GetLine(char * line )
{
  int i, j, terminated;

  /* 
   * This has been arbitrarily increased by 90 characters to allow 
   * for the increase in length of names from 8 to 16 characters,
   * the increase in length of numbers from 12 to 25 characters,
   * and possibility of having extra spaces between fields. 
   */
  const int length = 150;

  do {
    int c;

    iline++;
 
    i = 0; terminated = 0;
    while( i < length && EOF != ( c = getc( file ) ) ) {
      if( isEndOfLine( c, file ) ) {
	terminated = 1;
	break;
      }
      if( !isprint( c ) ) c = ' '; // Eliminate non-printing characters
      line[i] = c;
      
      i++;
    }
    if( !terminated ) {
      // The line wasn't terminated, or was possibly just terminated
      // after the 150th character.
      if( i == 0 ) {
	// Nothing was read. This is an error.
        fprintf( stderr, "Unexpected end-of-file at line %d.\n", iline );
        return READERROR;
      } else if( i == length ) { 
	// We read 150 characters without an error or a '\n'
	if( line[0] == '*' ) { 
	  // This is a comment line. There is no restriction on length.
	  // Just throw away characters til the end of line
          while( EOF != ( c = getc( file ) ) && c != '\n' ) ;
	} else {
	  // This is not a comment line. The only characters
	  // beyond the 150th column should be ' '
	  int extracrud = 0;
	  while( EOF != ( c = getc( file ) ) ) {
	    if( isEndOfLine( c, file ) ) break;
	    if( c != ' ' ) extracrud = 1;
	  }
	  if(  extracrud ) {
	    fprintf( stderr,
		     "Line %d has exceeded the maximum permissible characters.\n",
		     iline );
	    return mpssyntaxerr;
	  }
	} // end else this is not a comment line
      } // end else we read 150 characters
      // Otherwise, the file just doesn't have a terminating newline,
      // which is acceptable.
    }
    
    // Fill in the rest of line with spaces
    for( j = i; j < length; j++ ) line[j] = ' ';
    // Terminate the line
    line[length] = '\0';
  } while ( line[0] == '*' ); // Disgard comment lines
  
  if(line[0] == ' ') return DATALINE;
  else return HEADERLINE;
}


int MpsReader::ParseHeaderLine(char line[], char entry[] )
{
  // the comments in this function are 1-indexed. I.e. line[0] is
  // character one.
  int extra_crud = 0;

  this->string_copy(entry, &line[0],  16);  // characters  1 - 16 to entry1 
  extra_crud = !isOnlySpaces( line, 8, 60 );

  if( extra_crud ) {
    fprintf( stderr,
             "Extra characters outside prescribed fields at line %d.\n",
             iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}


int MpsReader::ParseHeaderLine2(char line[], char entry[] )
{
    int i = 0;
    char *token;
    int ierr =0;
    char tempLine[200];

    strncpy( tempLine, line, 200);

    /* Split the line delimited by space */
    token = strtok( tempLine, " ");

    /* Copy the token into the input argument */
    strcpy( entry, token);

    if( strlen( token) > 16 ) {
        fprintf( stderr,
                "Extra characters outside prescribed fields at line %d.\n",
                iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}


int MpsReader::ParseRowsLine2( char line[],  char code[], char name1[] )
{
    int i = 0;
    char *token;
    char *arrayOfTokens[2];
    int ierr =0;
    char tempLine[200];

    strncpy( tempLine, line, 200);

    token = strtok( tempLine, " ");

    // Split the extracted line into tokens delimited by space...
    arrayOfTokens[0] = token;
    arrayOfTokens[1] = strtok( NULL, " ");

    // Field 1: Row type
    if( arrayOfTokens[0] != NULL){
        strcpy(code, arrayOfTokens[0]);         
        }
    else{
        fprintf( stderr, "Empty row type field on line %d.\n", iline );
        return mpssyntaxerr;
        }

    // Field 2: Row name
    if( arrayOfTokens[1] != NULL){
        strcpy(name1, arrayOfTokens[1]);         
        }
    else{
        fprintf( stderr, "Empty row name field on line %d.\n", iline );
        return mpssyntaxerr;
        }
    
  return mpsok;
}


int MpsReader::ParseRowsLine( char line[],  char code[], char name1[] )
{
  int extra_crud = 0;
  assert( line[0] == ' ' );

  this->string_copy(code, &line[1], 2); // characters  2 -  3 to code  
  extra_crud = extra_crud || ' ' != line[3];

  this->string_copy(name1, &line[4], 8);    // characters  5 - 12 to name1 
  extra_crud = extra_crud || ' ' != line[12] || ' ' != line[13];
  if( isOnlySpaces( name1, 0, 7 ) ) { // name1 cannot be empty
    fprintf( stderr, "Empty row name field on line %d.\n", iline );
    return mpssyntaxerr;
  }
  
  extra_crud = extra_crud || !isOnlySpaces( line, 12, 60 );
  if( extra_crud ) {
    fprintf( stderr,
             "Extra characters outside prescribed fields at line %d.\n",
             iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}



int MpsReader::ParseBoundsLine( char line[], int& code, char name1[],
				char name2[], double * val )
{
  int extra_crud = 0;
  char codeStr[4];
  char valstr[16];
  int ierr;

  assert( line[0] == ' ' );
  this->string_copy(codeStr, &line[1], 2); // characters  2 -  3 to codeStr
  codeStr[0] = toupper( codeStr[0] ); codeStr[1] = toupper( codeStr[1] );
  extra_crud = extra_crud || ' ' != line[3];
    if( 0 == strcmp( codeStr, "LO" ) ) {
      code = kLowerBound;
    } else if ( 0 == strcmp( codeStr, "UP" ) ) {
      code = kUpperBound;
    } else if ( 0 == strcmp( codeStr, "FX" ) ) {
      code = kFixedBound;
      // fprintf( stderr, "This code does not support fixed variables.\n" );
      // return mpssyntaxerr;
    } else if ( 0 == strcmp( codeStr, "FR" ) ) {
      code = kFreeBound;
    } else if ( 0 == strcmp( codeStr, "MI" ) ) {
      code = kMInftyBound;
    } else if ( 0 == strcmp( codeStr, "PL" ) ) {
      code = kPInftyBound;
    } else {
      fprintf( stderr, "Bad type of bound specified on line %d.\n", iline );
      return mpssyntaxerr;
    }

  this->string_copy(name1, &line[4], 8);    // characters  5 - 12 to name1 
  // name1 is allowed to be empty.
  extra_crud = extra_crud || ' ' != line[12] || ' ' != line[13];

  this->string_copy(name2, &line[14], 8);   // characters 15 - 22 to name2 
  extra_crud = extra_crud || ' ' != line[22] || ' ' != line[23];
  if( isOnlySpaces( name2, 0, 7 ) ) { // name2 cannot be empty
    fprintf( stderr, "Empty second name field on line %d.\n", iline );
    return mpssyntaxerr;
  }
  switch( code ) {
  case kFreeBound:
  case kMInftyBound:
  case kPInftyBound:
    // this is it, nothing else may be specified
    extra_crud = extra_crud || !isOnlySpaces( line, 22, 60 );
    break;
  default:
    // A value must be specified.
    this->string_copy(valstr, &line[24], 12); 
    // characters 25 - 36 to valstr
    *val = asDouble( valstr, 12, ierr );

    if( 0 != ierr ) {
      fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
      return mpssyntaxerr;
    }
	
    extra_crud = extra_crud || !isOnlySpaces( line, 36, 60 );

    break;
  }
  if( extra_crud ) {
    fprintf( stderr,
             "Extra characters outside prescribed fields at line %d.\n",
             iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}


int MpsReader::ParseBoundsLine2( char line[], int& code, char name1[],
				char name2[], double * val )
    {
    int i = 0;
    char *token;
    char *arrayOfTokens[4];
    int ierr =0;
    char tempLine[200];
    char *endptr;
    char codeStr[4];
    bool boundsLabelPresent = false;

    *val = 0.0;

    strncpy( tempLine, line, 200);

    token = strtok( tempLine, " ");
    arrayOfTokens[0] = token;

    // Split the extracted line into tokens delimited by space...
    for( i = 1; i < 4; i++){       
        arrayOfTokens[i] = strtok( NULL, " ");
        }

    // Field 1: Specifies the types of bound
    if( arrayOfTokens[0] != NULL){
        strcpy(codeStr, arrayOfTokens[0]);         
        }
    else{
        fprintf( stderr, "Empty bound type on line %d.\n", iline );
        return mpssyntaxerr;
        }

    /* Convert the code to upper case if needed */
    codeStr[0] = toupper( codeStr[0] );
    codeStr[1] = toupper( codeStr[1] );

    if( 0 == strcmp( codeStr, "LO" ) ) {
      code = kLowerBound;
    } else if ( 0 == strcmp( codeStr, "UP" ) ) {
      code = kUpperBound;
    } else if ( 0 == strcmp( codeStr, "FX" ) ) {
      code = kFixedBound;
    } else if ( 0 == strcmp( codeStr, "FR" ) ) {
      code = kFreeBound;
    } else if ( 0 == strcmp( codeStr, "MI" ) ) {
      code = kMInftyBound;
    } else if ( 0 == strcmp( codeStr, "PL" ) ) {
      code = kPInftyBound;
    } else {
      fprintf( stderr, "Bad type of bound specified on line %d.\n", iline );
      return mpssyntaxerr;
    }

	// The presence of kFreeBound, kMInftyBound or kPInftyBound and 3 tokens indicates
	// that a bounds label is present.
	if( (code == kFreeBound) || (code == kMInftyBound) || (code == kPInftyBound)){
		if( arrayOfTokens[2] != NULL)
			boundsLabelPresent = true;
		}
	// The presence of kLowerBound, kUpperBound or kFixedBound and 4 tokens indicates 
	// that a bounds label is present
	else if( (code == kLowerBound) || (code == kUpperBound) || (code == kFixedBound)){
		if( arrayOfTokens[3] != NULL)
			boundsLabelPresent = true;
		}

    // Field 2: Bounds Label
    if( boundsLabelPresent){
	if( arrayOfTokens[1] != NULL){
		strcpy(name1, arrayOfTokens[1]);        
		}
	else{
		fprintf( stderr, "Empty first name field on line %d.\n", iline );
		return mpssyntaxerr;
		}
	}

	// Set the token index dynamically depending on whether or not the Bounds label is present
	int tokIndex = 1;
	if( boundsLabelPresent){
		tokIndex = 2;
		}

    // Field 3: Column Label
    if( arrayOfTokens[tokIndex] != NULL){
        strcpy(name2, arrayOfTokens[tokIndex]); 
		tokIndex++;
        }
    else{
        fprintf( stderr, "Empty second name field on line %d.\n", iline );
        return mpssyntaxerr;
        }

  switch( code ) {
  case kFreeBound:
  case kMInftyBound:
  case kPInftyBound:
    // this is it, nothing else may be specified
    break;

  default:

    // Field 4 (Optional): Bound value
    if( arrayOfTokens[tokIndex] != NULL){
        *val = strtod( arrayOfTokens[tokIndex], &endptr );

        if( endptr[0] != ' ' && endptr[0] != '\0')
            ierr = 1; // This works because we have already tokenized based on space delimiters

		if( 0 != ierr ){
            fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
            return mpssyntaxerr;
            }
        }
    break;
  }

return mpsok;
}

/*
//szhu - extend to character 60
int MpsReader::ParseBoundsLine( char line[], int& code, char name1[],
				char name2[], double * val )
{
  int extra_crud = 0;
  char codeStr[4];
  char valstr[36];  //extend to character 60
  int ierr;

  assert( line[0] == ' ' );
  this->string_copy(codeStr, &line[1], 2); // characters  2 -  3 to codeStr
  codeStr[0] = toupper( codeStr[0] ); codeStr[1] = toupper( codeStr[1] );
  extra_crud = extra_crud || ' ' != line[3];
    if( 0 == strcmp( codeStr, "LO" ) ) {
      code = kLowerBound;
    } else if ( 0 == strcmp( codeStr, "UP" ) ) {
      code = kUpperBound;
    } else if ( 0 == strcmp( codeStr, "FX" ) ) {
      code = kFixedBound;
      // fprintf( stderr, "This code does not support fixed variables.\n" );
      // return mpssyntaxerr;
    } else if ( 0 == strcmp( codeStr, "FR" ) ) {
      code = kFreeBound;
    } else if ( 0 == strcmp( codeStr, "MI" ) ) {
      code = kMInftyBound;
    } else if ( 0 == strcmp( codeStr, "PL" ) ) {
      code = kPInftyBound;
    } else {
      fprintf( stderr, "Bad type of bound specified on line %d.\n", iline );
      return mpssyntaxerr;
    }

  this->string_copy(name1, &line[4], 8);    // characters  5 - 12 to name1 
  // name1 is allowed to be empty.
  extra_crud = extra_crud || ' ' != line[12] || ' ' != line[13];

  this->string_copy(name2, &line[14], 8);   // characters 15 - 22 to name2 
  extra_crud = extra_crud || ' ' != line[22] || ' ' != line[23];
  if( isOnlySpaces( name2, 0, 7 ) ) { // name2 cannot be empty
    fprintf( stderr, "Empty second name field on line %d.\n", iline );
    return mpssyntaxerr;
  }
  switch( code ) {
  case kFreeBound:
  case kMInftyBound:
  case kPInftyBound:
    // this is it, nothing else may be specified
    extra_crud = extra_crud || !isOnlySpaces( line, 22, 60 );
    break;
  default:
    // A value must be specified.
	this->string_copy(valstr, &line[24], 36); 
    // characters 25 - 60 to valstr
    *val = asDouble( valstr, 36, ierr );

    if( 0 != ierr ) {
      fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
      return mpssyntaxerr;
    }

    break;
  }
  if( extra_crud ) {
    fprintf( stderr,
             "Extra characters outside prescribed fields at line %d.\n",
             iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}

*/


int MpsReader::ParseDataLine( char line[],  char code[],
                              char name1[], char name2[], double * val1,
                              int& hasSecondValue,
                              char name3[], double * val2)
{
  char            valstr1[16], valstr2[16];
  
  int extra_crud = 0, ierr;
  *val1 = 0.0;
  *val2 = 0.0;
  
  assert( line[0] == ' ' );
  this->string_copy(code, &line[1], 2); // characters  2 -  3 to code 
  extra_crud = extra_crud || ' ' != line[3];

  this->string_copy(name1, &line[4], 8);    // characters  5 - 12 to name1
  extra_crud = extra_crud || ' ' != line[12] || ' ' != line[13];
  // name1 is allowed to be empty.

  this->string_copy(name2, &line[14], 8);   // characters 15 - 22 to name2
  extra_crud = extra_crud || ' ' != line[22] || ' ' != line[23];
  if( isOnlySpaces( name2, 0, 7 ) ) { // name2 cannot be empty
    fprintf( stderr, "Empty second name field on line %d.\n", iline );
    return mpssyntaxerr;
  }

  this->string_copy(valstr1, &line[24], 12); // characters 25 - 36 to valstr 
  extra_crud = extra_crud ||
    ' ' != line[36] || ' ' != line[37] || ' ' != line[38];
  *val1 = asDouble( valstr1, 12, ierr );

  if( 0 != ierr ) {
    fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
    return mpssyntaxerr;
  }

  if( !isOnlySpaces( line, 39, 46 ) ) {
    hasSecondValue = 1;
    this->string_copy(name3, &line[39], 8); 
    // characters 40 - 47 to name3 
    extra_crud = extra_crud || ' ' != line[47] || ' ' != line[48];
    
    this->string_copy(valstr2, &line[49], 12); 
    // characters 50 - 61 to valstr 
    
    *val2 = asDouble( valstr2, 12, ierr );
    if( 0 != ierr ) {
      fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
      return mpssyntaxerr;
    }
  } else {
    hasSecondValue = 0;
    extra_crud = extra_crud || !isOnlySpaces( line, 47, 60 );
  }
  if( extra_crud ) {
    fprintf( stderr,
             "Extra characters outside prescribed fields at line %d.\n",
             iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}

int MpsReader::ParseDataLine2( char line[],  char code[],
                              char name1[], char name2[], double * val1,
                              int& hasSecondValue,
                              char name3[], double * val2)
    {
    int i = 0;
    char *token;
    char *arrayOfTokens[5];
    int ierr =0;
    char tempLine[200];
    char *endptr;

    *val1 = 0.0;
    *val2 = 0.0;
    hasSecondValue = 0;

    //Name1 is optional
    bool hasName1 = true;
    int numberOfTokens = 0;
    int tokIndex = 0;

    strncpy( tempLine, line, 200);

    token = strtok( tempLine, " ");
    arrayOfTokens[0] = token;
    if( arrayOfTokens[0] != NULL) 
    	numberOfTokens++;

    // Split the extracted line into tokens delimited by space...
    for( i = 1; i < 5; i++){       
        arrayOfTokens[i] = strtok( NULL, " ");

	if( arrayOfTokens[i] != NULL)
		numberOfTokens++;
        }

    // An even number of tokens indicates that name1 is missing
    if( (numberOfTokens % 2) == 0){
    	hasName1 = false;
	}

    // Field 2: Column/RHS/Right-hand side range vector Identifier
    if( hasName1 && arrayOfTokens[tokIndex] != NULL){
        strcpy(name1, arrayOfTokens[tokIndex]);         
        tokIndex++;
	}

    // Field 3: Row identifier
    if( arrayOfTokens[tokIndex] != NULL){
        strcpy(name2, arrayOfTokens[tokIndex]);
	tokIndex++;
        }
    else{
        fprintf( stderr, "Empty second name field on line %d.\n", iline );
        return mpssyntaxerr;
        }

    // Field 4: Value of matrix coefficient specified by fields 2 and 3
    if( arrayOfTokens[tokIndex] != NULL){

        *val1 = strtod( arrayOfTokens[tokIndex], &endptr );
	tokIndex++;

        if( endptr[0] != ' ' && endptr[0] != '\0')
            ierr = 1; // This works because we have already tokenized based on space delimiters

		if( 0 != ierr ){
            fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
            return mpssyntaxerr;
            }
        }

    // Field 5 (Optional): Row identifier
    if( arrayOfTokens[tokIndex] != NULL){
		hasSecondValue = 1;
		strcpy(name3, arrayOfTokens[tokIndex]);
		tokIndex++;
        }

    // Field 6 (Optional): Value of matrix coefficient specified by fields 2 and 5
    if( arrayOfTokens[tokIndex] != NULL){
        *val2 = strtod( arrayOfTokens[tokIndex], &endptr );
	tokIndex++;

        if( endptr[0] != ' ' && endptr[0] != '\0')
            ierr = 1; // This works because we have already tokenized based on space delimiters

		if( 0 != ierr ){
            fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
            return mpssyntaxerr;
            }
        }
    return ierr;
    }

/*
//szhu - extend to character 60
int MpsReader::ParseDataLine( char line[],  char code[],
                              char name1[], char name2[], double * val1,
                              int& hasSecondValue,
                              char name3[], double * val2)
{
  char            valstr1[36], valstr2[36];
  
  int extra_crud = 0, ierr;
  *val1 = 0.0;
  *val2 = 0.0;
  
  assert( line[0] == ' ' );
  this->string_copy(code, &line[1], 2); // characters  2 -  3 to code 
  extra_crud = extra_crud || ' ' != line[3];

  this->string_copy(name1, &line[4], 8);    // characters  5 - 12 to name1
  extra_crud = extra_crud || ' ' != line[12] || ' ' != line[13];
  // name1 is allowed to be empty.

  this->string_copy(name2, &line[14], 8);   // characters 15 - 22 to name2
  extra_crud = extra_crud || ' ' != line[22] || ' ' != line[23];
  if( isOnlySpaces( name2, 0, 7 ) ) { // name2 cannot be empty
    fprintf( stderr, "Empty second name field on line %d.\n", iline );
    return mpssyntaxerr;
  }

  //first check if there is second value
  this->string_copy(valstr1, &line[24], 12); // characters 25 - 36 to valstr 
  extra_crud = extra_crud ||
    ' ' != line[36] || ' ' != line[37] || ' ' != line[38];
  *val1 = asDouble( valstr1, 12, ierr );

  if( !extra_crud ) //may have second value or follows strict character position rule
  {
	  if( !isOnlySpaces( line, 39, 46 ) ) {  //has second value
		hasSecondValue = 1;
		this->string_copy(name3, &line[39], 8); 
		// characters 40 - 47 to name3 
		extra_crud = extra_crud || ' ' != line[47] || ' ' != line[48];
    
		this->string_copy(valstr2, &line[49], 12); 
		// characters 50 - 61 to valstr 
    
		*val2 = asDouble( valstr2, 12, ierr );
		if( 0 != ierr ) {
		  fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
		  return mpssyntaxerr;
		}
	  } else {  //no second value, users follow strict character position rule
		hasSecondValue = 0;
		extra_crud = extra_crud || !isOnlySpaces( line, 47, 60 );
		if( 0 != ierr ) {
		  fprintf( stderr, "Value doesn't parse as number on line %d.\n", iline );
		  return mpssyntaxerr;
		}
	  }
  }
  else  //no second value, the ending character position could be extened to 60
  {
    hasSecondValue = 0;
	extra_crud = 0;
	this->string_copy(valstr1, &line[24], 36); // characters 25 - 60 to valstr 
	*val1 = asDouble( valstr1, 36, ierr );
  }

  if( extra_crud ) {
    fprintf( stderr,
             "Extra characters outside prescribed fields at line %d.\n",
             iline );
    return mpssyntaxerr;
  }
  
  return mpsok;
}

*/

int MpsRowTypeFromCode2( char code )
{

  switch ( code ) {
  case 'N' : case 'n' : return kFreeRow;    break;
  case 'L' : case 'l' : return kLessRow;    break;
  case 'G' : case 'g' : return kGreaterRow; break;
  case 'E' : case 'e' : return kEqualRow;  break;
  default: return kBadRowType; break;
  }
}

/*
int MpsRowTypeFromCode( char code[2] )
{
  char c;

  if ( code[0] != ' ' ) {
    // the first character is not a space. The second
    // character must be.
    c = code[0];
    if ( code[1] != ' ' ) return kBadRowType;
  } else {
    // The first character is a space. The second character must
    // not be.
    c = code[1];
    if ( code[0] != ' ' ) return kBadRowType;
  }
  switch ( c ) {
  case 'N' : case 'n' : return kFreeRow;    break;
  case 'L' : case 'l' : return kLessRow;    break;
  case 'G' : case 'g' : return kGreaterRow; break;
  case 'E' : case 'e' : return kEqualRow;  break;
  default: return kBadRowType; break;
  }
}

  */

/////////////////
// Output section
/////////////////

char * MpsReader::defaultOutputFilename( int& iErr )
{
  if(infilename == NULL) {
    // Why are we here, if there was no input file?
    fprintf( stderr,
	     "Apparently no input file, so can't construct an output file.\n");
    iErr = mpsioerr;
    return 0;
  }
 // Figure out what the output file name should be.
  // Was there a "qps" suffix?
  char * suffix = strstr(infilename,".qps");
  if( suffix == NULL ) { // no qps suffix. How about an .mps suffix
    suffix = strstr(infilename,".mps");
  }
  char * outfilename;
  if(suffix != NULL) {
    // apparently there is a suffix; strip it off
    int len = strlen(infilename) - strlen(suffix);
    outfilename = new char[ len + 5 ];
    strncpy(outfilename, infilename, len);
    outfilename[len] = '\0';
  } else {
    int len = strlen(infilename);
    outfilename = new char[ len + 5 ];
    strcpy(outfilename, infilename);
  }
  // append ".out"
  strcat(outfilename, ".out");

  iErr = mpsok;
  return outfilename;
}

void MpsReader::printSolution( double x[],      int     nx,
			       double xlow[],   char    ixlow[],
			       double xupp[],   char    ixupp[],
			       double gamma[],  double  phi[],
			       double y[],      int     mA,
			       double s[],      int mC,
			       double clow[],   char    iclow[],
			       double cupp[],   char    icupp[],
			       double lambda[], double  pi[],
			       double z[],      double  objective, 
			       int& iErr )
{
  /** output file */
  FILE * outfile;

  iErr = mpsunknownerr;

  {
    char * outfilename = this->defaultOutputFilename( iErr );
    if( 0 != iErr ) return;

    outfile = fopen( outfilename, "w");
    if (outfile == NULL) {
      fprintf(stderr,
	      "Unable to open output file '%s'; solution not printed.\n", 
	      outfilename);
      iErr = mpsfileopenerr;
      return;
    }
    if( gOoqpPrintLevel > 0 ) {
      printf("\nprinting solution:  input file name %s\n", infilename);
      printf("printing solution: output file name %s\n", outfilename);
    }
    delete [] outfilename;
  }
  
  fprintf(outfile, "Solution for '%s'\n\n", problemName); 
  fprintf(outfile, "Rows: %d,  Columns: %d\n\n", totalRows, totalCols);
  
  // print row information
  fprintf(outfile, "  PRIMAL VARIABLES\n\n");
  fprintf(outfile, "       Name      Value "
	  "            Lower Bound      Upper Bound      Multiplier\n\n");

  for(int i=0; i<nx; i++) {
    if(ixlow[i] == 0.0 && ixupp[i] == 0.0) {
      // no lower or upper bounds, just print variable name and value
      fprintf(outfile, " %5d %8s %15.8e\n", i, colInfo[i].name, x[i]);
    } else if(ixlow[i] == 1.0 && ixupp[i] == 0.0) {
      // lower bound but no upper bound
      fprintf(outfile, " %5d %8s %15.8e   %15.8e                    %15.8e\n", 
	      i, colInfo[i].name, x[i], xlow[i], gamma[i]);
    } else if(ixlow[i] == 0.0 && ixupp[i] == 1.0) {
      // upper bound but no lower bound
      fprintf(outfile, " %5d %8s %15.8e "
	      "                    %15.8e  %15.8e\n", 
	      i, colInfo[i].name, x[i], xupp[i], -phi[i]);
    } else {
      // both upper and lower bounds
      fprintf(outfile, " %5d %8s %15.8e   %15.8e  %15.8e  %15.8e\n", 
	      i, colInfo[i].name, x[i], 
	      xlow[i], xupp[i], gamma[i] - phi[i]);
    }
  }
  
  // print column information for equality constrained rows first,
  // then for the rows with upper and lower bounds
  int i_equalities=0, i_inequalities=0;
  
  fprintf(outfile, "\n\n CONSTRAINTS\n\n");

  // first the equality constraints 
    
  if(mA >0) {
    fprintf(outfile, " Equality Constraints: %d\n\n", mA);
    fprintf(outfile, "        Name      Multiplier\n");
    for(int i=0; i<mA; i++) {
      // search for the next equality constraint in the row list
      // (which is sequenced according to their order in the MPS
      // file)
      while(i_equalities < totalRows && 
	    rowInfo[i_equalities].kind != kEqualRow)
	i_equalities++;
      fprintf( outfile,
	       " %5d  %8s  %15.8e\n", i, rowInfo[i_equalities].name, y[i]);
      i_equalities++;
    }
  }

  // now the inequality constraints

  if(mC >0) {
    fprintf(outfile, " Inequality Constraints: %d\n\n", mC);
    fprintf(outfile, "       Name      Value "
	    "            Lower Bound      Upper Bound      Multiplier\n\n");

    for(int i=0; i<mC; i++) {
      // search for the next inequality constraint in the row list
      while(i_inequalities < totalRows && 
	    (rowInfo[i_inequalities].kind == kEqualRow || 
	     rowInfo[i_inequalities].kind == kFreeRow))
	i_inequalities++;
      if(iclow[i] == 1.0 && icupp[i] == 0.0) {
	// lower bound only
	fprintf(outfile,
		" %5d %8s %15.8e   %15.8e                    %15.8e\n", 
		i, rowInfo[i_inequalities].name, s[i], clow[i], lambda[i]);
      } else if(iclow[i] == 0.0 && icupp[i] == 1.0) {
	// upper bound only
      // upper bound but no lower bound
	fprintf(outfile, " %5d %8s %15.8e "
		"                    %15.8e  %15.8e\n", 
		i, rowInfo[i_inequalities].name, s[i], cupp[i], -pi[i]);
      } else if(iclow[i] == 1.0 && icupp[i] == 1.0) {
	// lower and upper bounds both
	fprintf(outfile, " %5d %8s %15.8e   %15.8e  %15.8e  %15.8e\n", 
		i, rowInfo[i_inequalities].name, s[i], 
		clow[i], cupp[i], z[i]);
      }
      i_inequalities++;
    }
  }
  fprintf( outfile, "\nObjective value: %g\n", objective - objminus );
 
  fclose(outfile);

  iErr = mpsok;
  return;
}

