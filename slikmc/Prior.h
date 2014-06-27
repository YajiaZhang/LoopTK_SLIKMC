/*
 * Prior.h
 *
 *  Created on: Apr 29, 2013
 *      Author: Yajia
 */

#ifndef PRIOR_H_
#define PRIOR_H_

#include <vector.h>
#include <PChain.h>
/**
 * @brief Interface for user defined priors.
 */
class Prior {
public:
	virtual ~Prior();
	/**
	 * @brief Override this function to define custom energy/prior function.
	 * @return Probability density in logarithm.
	 */
	virtual double evaluate( PChain* protein);
};

#endif /* PRIOR_H_ */
