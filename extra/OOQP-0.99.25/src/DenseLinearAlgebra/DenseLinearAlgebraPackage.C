/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DenseLinearAlgebraPackage.h"
#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "DeSymIndefSolver.h"
#include "SimpleVector.h"

DenseLinearAlgebraPackage * DenseLinearAlgebraPackage::soleInstance()
{
  static DenseLinearAlgebraPackage * la = 0;

  if( !la ) {
	la = new DenseLinearAlgebraPackage;
  }

  return la;
}

SymMatrix *
DenseLinearAlgebraPackage::newSymMatrix( int size, int /* nnz */ )
{
  return new DenseSymMatrix( size );
}

GenMatrix *
DenseLinearAlgebraPackage::newGenMatrix( int m, int n, int /* nnz */)
{
  return new DenseGenMatrix( m, n );
}
  
OoqpVector * DenseLinearAlgebraPackage::newVector( int n )
{
  return new SimpleVector(n);
}
 
void DenseLinearAlgebraPackage::whatami( char type[32] )
{
  char type_[] = "DenseLinearAlgebraPackage";

  int i = 0;
  
  type[0] = type_[0];
  while( type[i] != 0 && i < 31 ) {
    ++i;
    type[i] = type_[i];
  }
}

