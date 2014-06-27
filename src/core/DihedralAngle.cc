/*
 * DihedralAngle.cc
 *
 *  Created on: Jun 12, 2012
 *      Author: Yajia
 */

#include "DihedralAngle.h"
#include <iostream>
using namespace std;

void DihedralAngle::print() {
	cout << this->phi << "\t" << this->psi << endl;
}

void DihedralAngle::set( const DihedralAngle& da) {
	this->phi = da.phi;
	this->psi = da.psi;
	return;
}


