/*
 * OsiOoqpSolverInterface.h
 *
 *  Created on: Sep 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVERINTERFACE_OSIOOQPSOLVERINTERFACE_H_
#define SRC_SOLVERINTERFACE_OSIOOQPSOLVERINTERFACE_H_

#define USE_VECTOR

/** Coin */
#include "OsiSolverInterface.hpp"
/** Ooqp */
#include "QpGen.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"

class OsiOoqpSolverInterface: public OsiSolverInterface {

public:

	/**@name Solve methods */
	//@{
	/// Solve initial LP relaxation
	virtual void initialSolve();
	/// Resolve an LP relaxation after problem modification
	virtual void resolve();
	/// Not supported
	virtual void branchAndBound();
	//@}

	/**@name Methods returning info on how the solution process terminated */
	//@{
	/// Are there numerical difficulties?
	virtual bool isAbandoned() const;
	/// Is optimality proven?
	virtual bool isProvenOptimal() const;
	/// Is primal infeasibility proven?
	virtual bool isProvenPrimalInfeasible() const;
	/// Is dual infeasibility proven?
	virtual bool isProvenDualInfeasible() const;
	/// Iteration limit reached?
	virtual bool isIterationLimitReached() const;
	//@}


	/**@name Warm start methods */
	//@{
	virtual CoinWarmStart *getEmptyWarmStart() const;
	virtual CoinWarmStart* getWarmStart() const;
	virtual bool setWarmStart(const CoinWarmStart* warmstart);
	//@}

	/**@name Problem query methods */
	//@{
    /// Get the number of columns
    virtual int getNumCols() const;
    /// Get the number of rows
    virtual int getNumRows() const;
    /// Get the number of nonzero elements
    virtual int getNumElements() const;
    /// Get a pointer to an array[getNumCols()] of column lower bounds
    virtual const double * getColLower() const;
    /// Get a pointer to an array[getNumCols()] of column upper bounds
    virtual const double * getColUpper() const;
    /*! \brief Get a pointer to an array[getNumRows()] of row constraint senses.
      <ul>
		  <li>'L': <= constraint
		  <li>'E': =  constraint
		  <li>'G': >= constraint
		  <li>'R': ranged constraint
		  <li>'N': free constraint
      </ul>
    */
    virtual const char * getRowSense() const;
	/*! \brief Get a pointer to an array[getNumRows()] of row right-hand sides
	 <ul>
		 <li> if getRowSense()[i] == 'L' then
		 getRightHandSide()[i] == getRowUpper()[i]
		 <li> if getRowSense()[i] == 'G' then
		 getRightHandSide()[i] == getRowLower()[i]
		 <li> if getRowSense()[i] == 'R' then
		 getRightHandSide()[i] == getRowUpper()[i]
		 <li> if getRowSense()[i] == 'N' then
		 getRightHandSide()[i] == 0.0
	 </ul>
	 */
    virtual const double * getRightHandSide() const;
	/*! \brief Get a pointer to an array[getNumRows()] of row ranges.
	 <ul>
		 <li> if getRowSense()[i] == 'R' then
		 getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
		 <li> if getRowSense()[i] != 'R' then
		 getRowRange()[i] is 0.0
	 </ul>
	 */
    virtual const double * getRowRange() const;
    /// Get a pointer to an array[getNumRows()] of row lower bounds
    virtual const double * getRowLower() const;
    /// Get a pointer to an array[getNumRows()] of row upper bounds
    virtual const double * getRowUpper() const;
	/*! \brief Get a pointer to an array[getNumCols()] of objective
	 function coefficients.
	 */
    virtual const double * getObjCoefficients() const;
	/*! \brief Get the objective function sense
	 -  1 for minimisation (default)
	 - -1 for maximisation
	 */
    virtual double getObjSense() const;
    /// Return true if the variable is continuous
    virtual bool isContinuous(int colIndex) const;
    /// Get a pointer to a row-wise copy of the matrix
    virtual const CoinPackedMatrix * getMatrixByRow() const;
    /// Get a pointer to a column-wise copy of the matrix
    virtual const CoinPackedMatrix * getMatrixByCol() const;
    /// Get the solver's value for infinity
    virtual double getInfinity() const;
    //@}

    /**@name Solution query methods */
	//@{
	/// Get a pointer to an array[getNumCols()] of primal variable values
	virtual const double * getColSolution() const;
	/// Get pointer to array[getNumRows()] of dual variable values
	virtual const double * getRowPrice() const;
	/// Get a pointer to an array[getNumCols()] of reduced costs
	virtual const double * getReducedCost() const;
	/** Get a pointer to array[getNumRows()] of row activity levels.
	 The row activity for a row is the left-hand side evaluated at the
	 current solution.
	 */
	virtual const double * getRowActivity() const;
	/// Get the objective function value.
	virtual double getObjValue() const;
	/** Get the number of iterations it took to solve the problem (whatever
	 `iteration' means to the solver).
	 */
	virtual int getIterationCount() const;
	/** Get as many dual rays as the solver can provide. In case of proven
	 primal infeasibility there should (with high probability) be at least
	 one.

	 The first getNumRows() ray components will always be associated with
	 the row duals (as returned by getRowPrice()). If \c fullRay is true,
	 the final getNumCols() entries will correspond to the ray components
	 associated with the nonbasic variables. If the full ray is requested
	 and the method cannot provide it, it will throw an exception.

	 \note
	 Implementors of solver interfaces note that the double pointers in
	 the vector should point to arrays of length getNumRows() (fullRay =
	 false) or (getNumRows()+getNumCols()) (fullRay = true) and they should
	 be allocated with new[].

	 \note
	 Clients of solver interfaces note that it is the client's
	 responsibility to free the double pointers in the vector using
	 delete[]. Clients are reminded that a problem can be dual and primal
	 infeasible.
	 */
	virtual std::vector<double*> getDualRays(int maxNumRays, bool fullRay =
			false) const;
	/** Get as many primal rays as the solver can provide. In case of proven
	 dual infeasibility there should (with high probability) be at least
	 one.

	 \note
	 Implementors of solver interfaces note that the double pointers in
	 the vector should point to arrays of length getNumCols() and they
	 should be allocated with new[].

	 \note
	 Clients of solver interfaces note that it is the client's
	 responsibility to free the double pointers in the vector using
	 delete[]. Clients are reminded that a problem can be dual and primal
	 infeasible.
	 */
	virtual std::vector<double*> getPrimalRays(int maxNumRays) const;
	//@}

	//-------------------------------------------------------------------------
	/**@name Methods to modify the objective, bounds, and solution

	 For functions which take a set of indices as parameters
	 (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
	 \c setRowSetTypes()), the parameters follow the C++ STL iterator
	 convention: \c indexFirst points to the first index in the
	 set, and \c indexLast points to a position one past the last index
	 in the set.

	 */
	//@{
	/** Set an objective function coefficient */
	virtual void setObjCoeff(int elementIndex, double elementValue);
	/** Set the objective function sense.

	 Use 1 for minimisation (default), -1 for maximisation.

	 \note
	 Implementors note that objective function sense is a parameter of
	 the OSI, not a property of the problem. Objective sense can be
	 set prior to problem load and should not be affected by loading a
	 new problem.
	 */
	virtual void setObjSense(double s);

	/** Set a single column lower bound.
	 Use -getInfinity() for -infinity. */
	virtual void setColLower(int elementIndex, double elementValue);
	/** Set a single column upper bound.
	 Use getInfinity() for infinity. */
	virtual void setColUpper(int elementIndex, double elementValue);
	/** Set a single row lower bound.
	 Use -getInfinity() for -infinity. */
	virtual void setRowLower(int elementIndex, double elementValue);
	/** Set a single row upper bound.
	 Use getInfinity() for infinity. */
	virtual void setRowUpper(int elementIndex, double elementValue);
	/** Set the type of a single row */
	virtual void setRowType(int index, char sense, double rightHandSide,
			double range);
	/** Set the primal solution variable values

	 colsol[getNumCols()] is an array of values for the primal variables.
	 These values are copied to memory owned by the solver interface
	 object or the solver.  They will be returned as the result of
	 getColSolution() until changed by another call to setColSolution() or
	 by a call to any solver routine.  Whether the solver makes use of the
	 solution in any way is solver-dependent.
	 */
	virtual void setColSolution(const double *colsol);

	/** Set dual solution variable values

	 rowprice[getNumRows()] is an array of values for the dual variables.
	 These values are copied to memory owned by the solver interface
	 object or the solver.  They will be returned as the result of
	 getRowPrice() until changed by another call to setRowPrice() or by a
	 call to any solver routine.  Whether the solver makes use of the
	 solution in any way is solver-dependent.
	 */
	virtual void setRowPrice(const double * rowprice);
	//@}

	//-------------------------------------------------------------------------
	/**@name Methods to set variable type */
	//@{
	/** Set the index-th variable to be a continuous variable */
	virtual void setContinuous(int index);
	/** Set the index-th variable to be an integer variable */
	virtual void setInteger(int index);
	//@}

	/**@name Methods to modify the constraint system.*/
	//@{
	/** Add a column (primal variable) to the problem. */
	virtual void addCol(const CoinPackedVectorBase& vec, const double collb,
			const double colub, const double obj);
	/** \brief Remove a set of columns (primal variables) from the
	 problem.

	 The solver interface for a basis-oriented solver will maintain valid
	 warm start information if all deleted variables are nonbasic.
	 */
	virtual void deleteCols(const int num, const int * colIndices);
	/** Add a row (constraint) to the problem. */
	virtual void addRow(const CoinPackedVectorBase& vec, const double rowlb,
			const double rowub);
	/** Add a row (constraint) to the problem. */
	virtual void addRow(const CoinPackedVectorBase& vec, const char rowsen,
			const double rowrhs, const double rowrng);
	/** \brief Delete a set of rows (constraints) from the problem.

	 The solver interface for a basis-oriented solver will maintain valid
	 warm start information if all deleted rows are loose.
	 */
	virtual void deleteRows(const int num, const int * rowIndices);
	//@}

	/**@name Methods for problem input and output */
	//@{
	/*! \brief Load in a problem by copying the arguments. The constraints on
	 the rows are given by lower and upper bounds.

	 If a pointer is 0 then the following values are the default:
	 <ul>
		 <li> <code>colub</code>: all columns have upper bound infinity
		 <li> <code>collb</code>: all columns have lower bound 0
		 <li> <code>rowub</code>: all rows have upper bound infinity
		 <li> <code>rowlb</code>: all rows have lower bound -infinity
		 <li> <code>obj</code>: all variables have 0 objective coefficient
	 </ul>

	 Note that the default values for rowub and rowlb produce the
	 constraint -infty <= ax <= infty. This is probably not what you want.
	 */
	virtual void loadProblem(const CoinPackedMatrix& matrix,
			const double* collb, const double* colub, const double* obj,
			const double* rowlb, const double* rowub);

	/*! \brief Load in a problem by assuming ownership of the arguments.
	 The constraints on the rows are given by lower and upper bounds.

	 For default argument values see the matching loadProblem method.

	 \warning
	 The arguments passed to this method will be freed using the
	 C++ <code>delete</code> and <code>delete[]</code> functions.
	 */
	virtual void assignProblem(CoinPackedMatrix*& matrix, double*& collb,
			double*& colub, double*& obj, double*& rowlb, double*& rowub);

	/*! \brief Load in a problem by copying the arguments.
	 The constraints on the rows are given by sense/rhs/range triplets.

	 If a pointer is 0 then the following values are the default:
	 <ul>
		 <li> <code>colub</code>: all columns have upper bound infinity
		 <li> <code>collb</code>: all columns have lower bound 0
		 <li> <code>obj</code>: all variables have 0 objective coefficient
		 <li> <code>rowsen</code>: all rows are >=
		 <li> <code>rowrhs</code>: all right hand sides are 0
		 <li> <code>rowrng</code>: 0 for the ranged rows
	 </ul>

	 Note that the default values for rowsen, rowrhs, and rowrng produce the
	 constraint ax >= 0.
	 */
	virtual void loadProblem(const CoinPackedMatrix& matrix,
			const double* collb, const double* colub, const double* obj,
			const char* rowsen, const double* rowrhs, const double* rowrng);

	/*! \brief Load in a problem by assuming ownership of the arguments.
	 The constraints on the rows are given by sense/rhs/range triplets.

	 For default argument values see the matching loadProblem method.

	 \warning
	 The arguments passed to this method will be freed using the
	 C++ <code>delete</code> and <code>delete[]</code> functions.
	 */
	virtual void assignProblem(CoinPackedMatrix*& matrix, double*& collb,
			double*& colub, double*& obj, char*& rowsen, double*& rowrhs,
			double*& rowrng);

	/*! \brief Load in a problem by copying the arguments. The constraint
	 matrix is is specified with standard column-major
	 column starts / row indices / coefficients vectors.
	 The constraints on the rows are given by lower and upper bounds.

	 The matrix vectors must be gap-free. Note that <code>start</code> must
	 have <code>numcols+1</code> entries so that the length of the last column
	 can be calculated as <code>start[numcols]-start[numcols-1]</code>.

	 See the previous loadProblem method using rowlb and rowub for default
	 argument values.
	 */
	virtual void loadProblem(const int numcols, const int numrows,
			const CoinBigIndex * start, const int* index, const double* value,
			const double* collb, const double* colub, const double* obj,
			const double* rowlb, const double* rowub);

	/*! \brief Load in a problem by copying the arguments. The constraint
	 matrix is is specified with standard column-major
	 column starts / row indices / coefficients vectors.
	 The constraints on the rows are given by sense/rhs/range triplets.

	 The matrix vectors must be gap-free. Note that <code>start</code> must
	 have <code>numcols+1</code> entries so that the length of the last column
	 can be calculated as <code>start[numcols]-start[numcols-1]</code>.

	 See the previous loadProblem method using sense/rhs/range for default
	 argument values.
	 */
	virtual void loadProblem(const int numcols, const int numrows,
			const CoinBigIndex * start, const int* index, const double* value,
			const double* collb, const double* colub, const double* obj,
			const char* rowsen, const double* rowrhs, const double* rowrng);

	/*! \brief Write the problem in MPS format to the specified file.

	 If objSense is non-zero, a value of -1.0 causes the problem to be
	 written with a maximization objective; +1.0 forces a minimization
	 objective. If objSense is zero, the choice is left to the implementation.
	 */
	virtual void writeMps(const char *filename, const char *extension = "mps",
			double objSense=0.0) const;
	//@}

	/**@name Constructors and destructors */
	//@{
	/** default constructor */
	OsiOoqpSolverInterface();

    /** Clone

      The result of calling clone(false) is defined to be equivalent to
      calling the default constructor OsiSolverInterface().
    */
    virtual OsiOoqpSolverInterface * clone(bool copyData = true) const;

    /** Copy constructor */
    OsiOoqpSolverInterface(const OsiOoqpSolverInterface &);

    /** Assignment operator */
    OsiOoqpSolverInterface & operator=(const OsiOoqpSolverInterface& rhs);

	/** default destructor */
	virtual ~OsiOoqpSolverInterface();
	//@}

protected:

	/**@name Protected methods */
	//@{
	/** Apply a row cut (append to the constraint matrix). */
	virtual void applyRowCut(const OsiRowCut & rc);

	/** Apply a column cut (adjust the bounds of one or more variables). */
	virtual void applyColCut(const OsiColCut & cc);

protected:

	/**@name Osi problem data */
	//@{
	double obj_sense_;       /**< objective sense */
	CoinPackedMatrix * mat_; /**< constraint matrix */
	std::vector<double> clbd_;  /**< column lower bounds */
	std::vector<double> cubd_;  /**< column upper bounds */
	std::vector<double> obj_;   /**< objective coefficients */
	std::vector<double> rlbd_;  /**< row lower bounds */
	std::vector<double> rubd_;  /**< row upper bounds */
	std::vector<char> sense_;   /**< row sense */
	std::vector<double> rhs_;   /**< right-hand side */
	std::vector<double> range_; /**< range value = rubd_ - rlbd_ */
	//@}

	/**@name Osi solution data */
	//@{
	int status_;    /**< solution status */
	double objval_; /**< objective value */
	std::vector<double> x_;        /**< primal solution */
	std::vector<double> price_;    /**< row price */
	std::vector<double> reduced_;  /**< reduced cost */
	std::vector<double> activity_; /**< row lhs */
	//@}

	int nIters_; /**< number of iterations */

public:

	/**@name Load problem with a quadratic objective function. The quadratic objective
	 * function is given as a lower triangular matrix.
	 */
	//@{
	virtual void loadProblem(const CoinPackedMatrix& matrix,
			const double* collb, const double* colub, const double* obj,
			int nnzQ, const int* irowQ, const int* jcolQ, const double* dQ,
			const double* rowlb, const double* rowub);

	virtual void assignProblem(CoinPackedMatrix*& matrix, double*& collb,
			double*& colub, double*& obj,
			int nnzQ, int*& irowQ, int*& jcolQ, double*& dQ,
			double*& rowlb, double*& rowub);

	virtual void loadProblem(const CoinPackedMatrix& matrix,
			const double* collb, const double* colub, const double* obj,
			int nnzQ, const int* irowQ, const int* jcolQ, const double* dQ,
			const char* rowsen, const double* rowrhs, const double* rowrng);

	virtual void assignProblem(CoinPackedMatrix*& matrix, double*& collb,
			double*& colub, double*& obj,
			int nnzQ, int*& irowQ, int*& jcolQ, double*& dQ,
			char*& rowsen, double*& rowrhs, double*& rowrng);

	virtual void loadProblem(const int numcols, const int numrows,
			const CoinBigIndex * start, const int* index, const double* value,
			const double* collb, const double* colub, const double* obj,
			int nnzQ, const int* irowQ, const int* jcolQ, const double* dQ,
			const double* rowlb, const double* rowub);

	virtual void loadProblem(const int numcols, const int numrows,
			const CoinBigIndex * start, const int* index, const double* value,
			const double* collb, const double* colub, const double* obj,
			int nnzQ, const int* irowQ, const int* jcolQ, const double* dQ,
			const char* rowsen, const double* rowrhs, const double* rowrng);
	//@}

	/** initialize all from a matrix */
	void initializeProblem(const CoinPackedMatrix& matrix);

	/** free cached result vectors */
	void freeCachedResults();

	/** free cached data */
	void freeCachedData();

	/** free cached Ooqp objects */
	void freeCachedOoqp();

protected:

	/** core part of solve */
	void gutsOfSolve();

	/** set quadratic objective function in a lower triangular matrix */
	void setQuadraticObjective(
			int nnz,            /**< [in] number of nonzeros */
			const int* irow,    /**< [in] row indices */
			const int* jcol,    /**< [in] column indices */
			const double* value /**< [in] matrix elements */);

	/** convert Osi data to Ooqp data */
	void convertOsiToOoqp(QpGen*& qpgen, QpGenData*& prob);

protected:

	/**@name OOQP objects */
	//@{
	QpGen* qpgen_;
	QpGenData* prob_;
	QpGenVars* vars_;
	QpGenResiduals* resid_;
	Solver* solver_;
	//@}

	/**@name Lower triangular matrix of the quadratic objective function */
	//@{
	int nnzQ_;
	std::vector<int> irowQ_;
	std::vector<int> jcolQ_;
	std::vector<double> dQ_;
	//@}

	/**
	 * The original matrix may contain empty or duplicate rows.
	 * So, we remove them and pass the reduced matrix to the solver. Doing so may tighten some row bounds.
	 */
	std::vector<int> ia_to_orglow_; /**< map matrix A row index to the original row index with lower bound */
	std::vector<int> ia_to_orgupp_; /**< map matrix A row index to the original row index with upper bound */
	std::vector<int> iclow_to_org_; /**< map matrix C row index with lower bound to the original row index */
	std::vector<int> icupp_to_org_; /**< map matrix C row index with upper bound to the original row index */

	bool updated_; /**< indicate if the model is updated */
};

#endif /* SRC_SOLVERINTERFACE_OSIOOQPSOLVERINTERFACE_H_ */
