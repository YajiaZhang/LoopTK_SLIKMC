/*
 * Prior.cc
 *
 *  Created on: Apr 29, 2013
 *      Author: Yajia
 */

#include "Prior.h"
#include <iostream>
using namespace std;

Prior::~Prior() {
}

double Prior::evaluate(PChain* protein) {
	cerr << "Warning: No self-defined evaluate function" << endl;
	return 0;
}
