/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseLinearAlgebraPackage.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "DoubleMatrix.h"

SparseLinearAlgebraPackage * SparseLinearAlgebraPackage::soleInstance()
{
  static SparseLinearAlgebraPackage * la = 0;
  
  if( !la ) la = new SparseLinearAlgebraPackage;

  return la;
}

SymMatrix * SparseLinearAlgebraPackage::newSymMatrix( int size, int nnz)
{
  return new SparseSymMatrix( size, nnz );
}

GenMatrix * SparseLinearAlgebraPackage::newGenMatrix( int m, int n,
							  int nnz )
{
  return new SparseGenMatrix( m, n, nnz );
}
  
OoqpVector * SparseLinearAlgebraPackage::newVector( int n )
{
  return new SimpleVector(n);
}
 
void SparseLinearAlgebraPackage::whatami( char type[32] )
{
  char type_[] = "SparseLinearAlgebraPackage";

  int i = 0;
  
  type[0] = type_[0];
  while( type[i] != 0 && i < 31 ) {
    ++i;
    type[i] = type_[i];
  }
}

