/**
 * DspOsiGrb.h
 *
 */
 
#ifndef SRC_SOLVERINTERFACE_DSPOSIGRB_H_
#define SRC_SOLVERINTERFACE_DSPOSIGRB_H_

#ifdef DSP_HAS_GRB

#include "gurobi_c.h"
#include "OsiGrbSolverInterface.hpp"
#include "SolverInterface/DspOsi.h"

#define debugMessage \
  if (false)         \
  printf

#define GUROBI_CALL(m, x)                                                            \
  do {                                                                               \
    int _retcode;                                                                    \
    if ((_retcode = (x)) != 0) {                                                     \
      char s[1001];                                                                  \
      if (grb_->getEnvironmentPtr()) {                                       \
        sprintf(s, "Error <%d> from GUROBI function call: ", _retcode);              \
        strncat(s, GRBgeterrormsg(grb_->getEnvironmentPtr()), 1000);         \
      } else                                                                         \
        sprintf(s, "Error <%d> from GUROBI function call (no license?).", _retcode); \
      debugMessage("%s:%d: %s", __FILE__, __LINE__, s);                              \
      throw CoinError(s, m, "DspOsiGrb", __FILE__, __LINE__);            \
    }                                                                                \
  } while (false)

class DspOsiGrb : public DspOsi {
public:

	/** default constructor */
	DspOsiGrb() {
		si_ = new OsiGrbSolverInterface();
		grb_ = dynamic_cast<OsiGrbSolverInterface*>(si_);
	}

	/** copy constructor */
	DspOsiGrb(const DspOsiGrb& rhs) {
        si_ = new OsiGrbSolverInterface(*(rhs.grb_));
		grb_ = dynamic_cast<OsiGrbSolverInterface*>(si_);
	}

	/** clone constructor */
	virtual DspOsiGrb* clone() const {
		return new DspOsiGrb(*this);
	}

	/** destructor */
	virtual ~DspOsiGrb() {
		grb_ = NULL;
	}

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
		try{
			// need to delete previous coefficient first, otherwise, the value in mat is added to current model
			GUROBI_CALL("loadQuadraticObjective", GRBdelq(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
			
			if (mat.isColOrdered()) {
				for (int j = 0; j < mat.getMajorDim(); ++j) {
					for (int k = 0; k < mat.getVectorSize(j); ++k) {
						int i = mat.getIndices()[mat.getVectorStarts()[j] + k];
						double v = mat.getElements()[mat.getVectorStarts()[j] + k];
						int row[1]={i};
						int col[1]={j};
						double element[1];
						// if (i!=j){
						// 	element[0]=2*v;
						// }
						// else{
							element[0]=v;
						//}
                    	GUROBI_CALL("loadQuadraticObjective", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
                    	GUROBI_CALL("loadQuadraticObjective", GRBaddqpterms(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL),
                    	1, row, col, element));
					}
				}
			} else {
				GUROBI_CALL("loadQuadraticObjective", GRBdelq(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
				for (int i = 0; i < mat.getMajorDim(); ++i) {
					for (int k = 0; k < mat.getVectorSize(i); ++k) {
						int j = mat.getIndices()[mat.getVectorStarts()[i] + k];
						double v = mat.getElements()[mat.getVectorStarts()[i] + k];
						int row[1]={i};
						int col[1]={j};
						double element[1];
						// if (i!=j){
						// 	element[0]=2*v;
						// }
						// else{
							element[0]=v;
						//}
						
						GUROBI_CALL("loadQuadraticObjective", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
                    	GUROBI_CALL("loadQuadraticObjective", GRBaddqpterms(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL),
                    	1, row, col, element));	
					}
				}
			}
		}
		catch(const CoinError& e){
        	e.print();
		}
    }
    

	/** solve problem */
	virtual void solve() {
		try{
			if (si_->getNumIntegers() > 0)
				si_->branchAndBound();
			else
				si_->initialSolve();
			
		}
		catch(const CoinError& e){
        	e.print();
		}
    }

	virtual void writeMps(const char *filename){
		try{
			std::string f(filename);
  			std::string e("mps");
			std::string fullname = f + "." + e;
			GUROBI_CALL("writeMPS", GRBwrite(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL), const_cast< char * >(fullname.c_str())));
		}
		catch(const CoinError& e){
        	e.print();
		}
	}

	//#############################################################################
	// Methods to input a problem
	//#############################################################################
	/*
	virtual void loadProblem(const CoinPackedMatrix &matrix,
  		const double *collb, const double *colub,
  		const double *obj, const CoinPackedMatrix &qobj,
  		const double *rowlb, const double *rowub)
	{
		debugMessage("OsiGrbSolverInterface::loadProblem(1)(%p, %p, %p, %p, %p, %p)\n", (void *)&matrix, (void *)collb, (void *)colub, (void *)obj, (void *)& qobj, (void *)rowlb, (void *)rowub);
		si_->loadProblem(matrix, collb, colub, obj, rowlb, rowub);
		loadQuadraticObjective(qobj);
	}
	*/
	virtual void use_simplex() {
		try{
			GUROBI_CALL("use simplex", GRBsetintparam(grb_->getEnvironmentPtr(), GRB_INT_PAR_METHOD, 1));
		}
		catch(const CoinError& e){
        	e.print();
		}
	}

	virtual void use_barrier() {
		try{
			GUROBI_CALL("use barrier", GRBsetintparam(grb_->getEnvironmentPtr(), GRB_INT_PAR_METHOD, 2));
			GUROBI_CALL("use barrier", GRBsetintparam(grb_->getEnvironmentPtr(), GRB_INT_PAR_CROSSOVER, 0)); //disabled
			GUROBI_CALL("use barrier", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_BARCONVTOL, 1e-5));
		}
		catch(const CoinError& e){
        	e.print();
		}
	}

	/** solution statue */
	virtual int status() {
		
		int status = DSP_STAT_UNKNOWN;
		int stat;
		try{
        	GUROBI_CALL("status", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
        	GUROBI_CALL("status", GRBgetintattr(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL), "Status", &stat));
		}
		catch(const CoinError& e){
        	e.print();
			return status;
		}

        switch(stat) {
        case GRB_OPTIMAL:
        case GRB_SUBOPTIMAL:
            status = DSP_STAT_OPTIMAL;
            break;
        case GRB_INFEASIBLE:
        case GRB_CUTOFF:
            status = DSP_STAT_PRIM_INFEASIBLE;
            break;
        case GRB_INF_OR_UNBD:
		case GRB_UNBOUNDED:
            status = DSP_STAT_DUAL_INFEASIBLE;
            break;
        case GRB_USER_OBJ_LIMIT:
            status = DSP_STAT_LIM_PRIM_OBJ;
            break;
        case GRB_ITERATION_LIMIT:
            status = DSP_STAT_STOPPED_ITER;
            break;
        case GRB_NODE_LIMIT:
            status = DSP_STAT_STOPPED_NODE;
            break;
        case GRB_TIME_LIMIT:
            status = DSP_STAT_STOPPED_TIME;
            break;

        case GRB_INTERRUPTED:
            status = DSP_STAT_STOPPED_USER;
            break;
        case GRB_LOADED:
        case GRB_INPROGRESS:
        case GRB_NUMERIC:
            status = DSP_STAT_ABORT;
            break;
	    default:
            status = DSP_STAT_UNKNOWN;
            break;
		}
	
        return status;
	}

	/** get dual objective value */
	virtual double getDualObjValue() {
		try{
			double val;
        	GUROBI_CALL("getDualObjVal", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
			GUROBI_CALL("getDualObjVal", GRBgetdblattr(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL), GRB_DBL_ATTR_OBJVAL, &val));
			return val;
		}
		catch(const CoinError& e){
        	e.print();
			exit(1);
    	}
	}

	/** get number of branch-and-bound nodes explored */
	virtual int getNumNodes() {
		try{
			double node;
        	GUROBI_CALL("getNumNodes", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
       		GUROBI_CALL("getNumNodes", GRBgetdblattr(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL), GRB_DBL_ATTR_NODECOUNT, &node));       	
        	return (int)node;
		}
		catch(const CoinError& e){
        	e.print();
			exit(1);
    	}
	}

	/** set number of cores */
	virtual void setNumCores(int num) {
		try{
        	GUROBI_CALL("setNumCores", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
			GUROBI_CALL("setNumCores", GRBsetintparam(grb_->getEnvironmentPtr(), GRB_INT_PAR_THREADS, num));
		}
		catch(const CoinError& e){
        	e.print();
    	}
	}

	/** set time limit */
	virtual void setTimeLimit(double time) {
		try{
        	GUROBI_CALL("setTimeLimit", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
			GUROBI_CALL("setTimeLimit", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_TIMELIMIT, time));
		}
		catch(const CoinError& e){
        	e.print();
    	}
	}

	/** set node limit */
	virtual void setNodeLimit(double num) {
		try{
        	GUROBI_CALL("setNodeLimit", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
			GUROBI_CALL("setNodeLimit", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_NODELIMIT, num));
		}
		catch(const CoinError& e){
        	e.print();
    	}
	}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {
		try{
        	GUROBI_CALL("setRelMipGap", GRBupdatemodel(grb_->getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ALL)));
			GUROBI_CALL("setRelMipGap", GRBsetdblparam(grb_->getEnvironmentPtr(), GRB_DBL_PAR_MIPGAP, tol));
		}
		catch(const CoinError& e){
        	e.print();
    	}
	}

    OsiGrbSolverInterface* grb_;   
};

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSIGRB_H_ */