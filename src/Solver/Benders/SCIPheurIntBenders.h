/*
 * SCIPheurIntBenders.h
 * 
 * This implements a SCIP heuristic handler for integer Benders decomposition.
 * This class is used only when the second-stage has integer variables.
 */

#ifndef SCIPHEURINTBENDERS_H_
#define SCIPHEURINTBENDERS_H_

#include "scip/def.h"
#include "objscip/objheur.h"
#include "Model/DecModel.h"
#include "Solver/Benders/BdSub.h"

class SCIPheurIntBenders : public scip::ObjHeur
{
public:

    /** default constructor */
    SCIPheurIntBenders(SCIP* scip, const char* name) :
    ObjHeur(
        scip,
        name,
        "integer recourse evaluation",
        'R', /**< dispchar */
        1,   /**< high priority */
        1,   /**< freq: always run */
        1,   /**< freqofs: always run */
        -1,  /**< maxdepth: always run */
        SCIP_HEURTIMING_AFTERLPNODE,
        FALSE
    ),
    model_(NULL), bdsub_(NULL) {
        // nothing to do
    }

    /** default desctructor */
    virtual ~SCIPheurIntBenders() {}

    virtual SCIP_DECL_HEUREXEC(scip_exec);

protected:

	DecModel *  model_; /**< DecModel object */
	BdSub *     bdsub_; /**< pointer to cut generator */
};

#endif /* SCIPHEURINTBENDERS_H_ */