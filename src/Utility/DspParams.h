/*
 * DspParams.h
 *
 *  Created on: Oct 20, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_UTILITY_DSPPARAMS_H_
#define SRC_UTILITY_DSPPARAMS_H_

#include <limits>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

enum DSP_DD_MASTER_ALGO
{
	Simplex = 0,
	IPM,
	IPM_Feasible,
	DSBM, /**< doubly stabilized bundle method */
	Subgradient
};

enum DSP_BD_INIT_LB_ALGO
{
	SEPARATE_LP = 0,
	SEPARATE_MILP,
};

enum DSP_EXTERNAL_SOLVER {
	OsiCpx = 0,
	OsiScip,
	OsiOoqp,
	OsiClp,
	OsiGrb
};

enum DSP_DW_BRANCH {
	BRANCH_INT = 0,
	BRANCH_NONANT,
	BRANCH_NONANT2,
	BRANCH_DISJUNCTION_TEST,
};

/**
 * This class create, set and get parameters.
 */
template <class T>
class DspParam {
public:

	/** default constructor */
//	DspParam();

	/** default destructor */
//	virtual ~DspParam();

	/** create parameter */
	void createParam(string name, T const & value);

	/** delete parameter */
	void deleteParam(string name);

	/** set parameter */
	void setParam(string name, T const & value);

	/** get parameter */
	T getParam(string name) const;

private:
	unordered_map<string,T> params_;
};

/** Array parameter */
template <class T>
class DspPtrParam {

	typedef unordered_map<string,int> MapSize;
	typedef unordered_map<string,T*> MapParams;

public:

	/** default constructor */
//	DspPtrParam();

	/** default destructor */
//	virtual ~DspPtrParam();

	/** create parameter */
	void createParam(string name, int size = 0);

	/** set parameter size */
	void setParamSize(string name, int size);

	/** delete parameter */
	void deleteParam(string name);

	/** set parameter */
	void setParam(string name, int index, T const & value);

	/** get parameter array size */
	int getParamSize(string name) const;

	/** get parameter */
	T * getParam(string name) const;

private:
	MapSize   size_; /**< size of the parameter array */
	MapParams params_;
};

/** create parameter */
template<class T>
void DspParam<T>::createParam(string name, T const & value)
{
	if (params_.find(name) != params_.end())
		printf("WARNING: The parameter <%s> already exists.\n", name.c_str());
	else
		params_[name] = value;
}

/** delete parameter */
template<class T>
void DspParam<T>::deleteParam(string name)
{
	if (params_.find(name) != params_.end())
		params_.erase(name);
	else
		printf("WARNING: There is no parameter <%s>.\n", name.c_str());
}

/** set parameter */
template<class T>
void DspParam<T>::setParam(string name, T const & value)
{
	if (params_.find(name) != params_.end())
		params_[name] = value;
	else
		printf("WARNING: There is no parameter <%s>.\n", name.c_str());
}

/** get parameter */
template<class T>
T DspParam<T>::getParam(string name) const
{
	auto found = params_.find(name);
	if (found != params_.end())
		return found->second;
	else
	{
		printf("WARNING: There is no parameter <%s>.\n", name.c_str());
		return T();
	}
}

/** create parameter */
template<class T>
void DspPtrParam<T>::createParam(string name, int size)
{
	if (size_.find(name) != size_.end())
		printf("WARNING: The parameter <%s> already exists.\n", name.c_str());
	else
	{
		size_[name] = size;
		if (size > 0)
			params_[name] = new T [size];
		else
			params_[name] = NULL;
	}
}

/** set parameter size */
template<class T>
void DspPtrParam<T>::setParamSize(string name, int size)
{
	if (params_.find(name) != params_.end())
	{
		deleteParam(name);
		createParam(name, size);
	}
	else
		printf("WARNING: There is no parameter <%s>\n.", name.c_str());
}

template<class T>
void DspPtrParam<T>::deleteParam(string name)
{
	auto found1 = params_.find(name);
	auto found2 = size_.find(name);
	if (found1 != params_.end() &&
		found2 != size_.end())
	{
		if (found1->second != NULL)
		{
			delete [] found1->second;
			found1->second = NULL;
		}
		params_.erase(name);
		size_.erase(name);
	}
	else
		printf("WARNING: There is no parameter <%s>\n.", name.c_str());
}

/** set parameter */
template<class T>
void DspPtrParam<T>::setParam(string name, int index, T const & value)
{
	auto found = params_.find(name);
	if (found != params_.end())
		found->second[index] = value;
	else
		printf("WARNING: There is no parameter <%s>\n.", name.c_str());
}

/** get parameter array size */
template<class T>
int DspPtrParam<T>::getParamSize(string name) const
{
	auto found = size_.find(name);
	if (found != size_.end())
		return found->second;
	else
	{
		printf("WARNING: There is no parameter <%s>\n.", name.c_str());
		return -1;
	}
}

/** get parameter */
template<class T>
T * DspPtrParam<T>::getParam(string name) const
{
	auto found = params_.find(name);
	if (found != params_.end())
		return found->second;
	else
	{
		printf("WARNING: There is no parameter <%s>\n.", name.c_str());
		return NULL;
	}
}

class DspParams {
public:

	/** default constructor */
	DspParams();

	/** default destructor */
	virtual ~DspParams();

	/** read parameter file */
	void readParamFile(const char * param_file);

private:
	/** INITIALIZE */
	void initBoolParams();
	void initIntParams();
	void initDblParams();
	void initStrParams();
	void initBoolPtrParams();
	void initIntPtrParams();
	void initDblPtrParams();

public:
	/** SET */

	/** set boolean type parameter */
	void setBoolParam(string name, bool value)
	{
		BoolParams_.setParam(name, value);
	}

	/** set double type parameter */
	void setDblParam(string name, double value)
	{
		DblParams_.setParam(name, value);
	}

	/** set integer type parameter */
	void setIntParam(string name, int value)
	{
		IntParams_.setParam(name, value);
	}

	/** set string type parameter */
	void setStrParam(string name, string value)
	{
		StrParams_.setParam(name, value);
	}

	/** set boolean pointer parameter size */
	void setBoolPtrParamSize(string name, int size)
	{
		BoolPtrParams_.setParamSize(name, size);
	}

	/** set integer pointer parameter size */
	void setIntPtrParamSize(string name, int size)
	{
		IntPtrParams_.setParamSize(name, size);
	}

	/** set double pointer parameter size */
	void setDblPtrParamSize(string name, int size)
	{
		DblPtrParams_.setParamSize(name, size);
	}

	/** set boolean pointer type parameter */
	void setBoolPtrParam(string name, int index, bool value)
	{
		BoolPtrParams_.setParam(name, index, value);
	}

	/** set integer pointer type parameter */
	void setIntPtrParam(string name, int index, int value)
	{
		IntPtrParams_.setParam(name, index, value);
	}

	/** set double pointer type parameter */
	void setDblPtrParam(string name, int index, double value)
	{
		DblPtrParams_.setParam(name, index, value);
	}

	/** GET */

	/** get boolean type parameter */
	bool getBoolParam(string name) const {return BoolParams_.getParam(name);}

	/** get double type parameter */
	double getDblParam(string name) const {return DblParams_.getParam(name);}

	/** get integer type parameter */
	int getIntParam(string name) const {return IntParams_.getParam(name);}

	/** get string type parameter */
	string getStrParam(string name) const {return StrParams_.getParam(name);}

	/** get boolean pointer parameter size */
	int getBoolPtrParamSize(string name) {return BoolPtrParams_.getParamSize(name);}

	/** get integer pointer parameter size */
	int getIntPtrParamSize(string name) {return IntPtrParams_.getParamSize(name);}

	/** get double pointer parameter size */
	int getDblPtrParamSize(string name) {return DblPtrParams_.getParamSize(name);}

	/** get boolean pointer type parameter */
	bool * getBoolPtrParam(string name) const {return BoolPtrParams_.getParam(name);}

	/** get integer pointer type parameter */
	int * getIntPtrParam(string name) const {return IntPtrParams_.getParam(name);}

	/** get double pointer type parameter */
	double * getDblPtrParam(string name) const {return DblPtrParams_.getParam(name);}

private:

	DspParam<bool>    BoolParams_;
	DspParam<double>  DblParams_;
	DspParam<int>     IntParams_;
	DspParam<string>  StrParams_;
	DspPtrParam<bool>   BoolPtrParams_;
	DspPtrParam<int>    IntPtrParams_;
	DspPtrParam<double> DblPtrParams_;
};

#endif /* SRC_UTILITY_DSPPARAMS_H_ */
