/*
 * StoSmpsScenario.h
 *
 *  Created on: Oct 9, 2014
 *      Author: kibaekkim
 */

#ifndef STOSMPSSCENARIO_H_
#define STOSMPSSCENARIO_H_

#include <string>
#include <vector>

using namespace std;

class StoSmpsScenario
{
public:
	string name_;
	string parent_;
	string period_;
	double probability_;
	vector<string> row_;
	vector<string> col_;
	vector<double> val_;
	vector<string> bound_code_;
	vector<string> bound_name_;
	vector<string> bound_col_;
	vector<double> bound_val_;
public:
	StoSmpsScenario(string name, string parent, string period, double probability) :
		name_(name), parent_(parent), period_(period), probability_(probability) {}
	StoSmpsScenario(const char * name, const char * parent, const char * period, double probability) :
		name_(name), parent_(parent), period_(period), probability_(probability) {}
};


#endif /* STOSMPSSCENARIO_H_ */
