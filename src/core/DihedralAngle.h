/*
 * DihedralAngle.h
 *
 *  Created on: Jun 12, 2012
 *      Author: Yajia
 */

#ifndef DIHEDRALANGLE_H_
#define DIHEDRALANGLE_H_

#include <string.h>
#include <vector.h>
using namespace std;

class DihedralAngle {
public:
	double phi;
	double psi;
	void print();
	void set( const DihedralAngle& da);
};

#endif /* DIHEDRALANGLE_H_ */
