/*
 * RamachandranPlot.h
 *
 *  Created on: Jun 6, 2012
 *      Author: Yajia
 */

#ifndef RAMACHANDRANPLOT_H_
#define RAMACHANDRANPLOT_H_

#include <string.h>
#include <vector.h>
#include "PChain.h"
#include "DihedralAngle.h"
using namespace std;

/**
 * @brief The class RamachandranPlot is used to initialize and store database for Ramachandran plot, sample backbone dihedral angles according to residue type
 * and evaluate the probability of one conformation according to backbone structure.
 */
class RamachandranPlot {
public:
	/**
	 * @brief Constructor
	 * @param size number of grids to discretize the angle space [-180.0, 180.0]. By default, 10-degree is the width of a grid.
	 * @param file_pro library location for Proline
	 * @param file_pre_pro library location for pre-Prolines
	 * @param file_gly library location for Glycine
	 * @param file_generic library location for generic residues
	 */
	RamachandranPlot( int size = 36, char* file_pro = "../Data/RamachandranPlot/pro.txt", char* file_pre_pro = "../Data/RamachandranPlot/pre_pro.txt",
			char* file_gly = "../Data/RamachandranPlot/gly.txt", char* file_generic = "../Data/RamachandranPlot/gen.txt");

	/**
	 * @brief Destructor
	 */
	virtual ~RamachandranPlot();

	/**
	 * @brief Randomly generate a pair of backbone dihedral angles of specific type.
	 * @param type residue type {PRO, PRE_PRO, GLY, GENERIC}
	 * @return a pair of dihedral angles
	 */
	DihedralAngle* getRandomDihedralAngle( int type);

	/**
	 * @brief Randomly generate a pair of backbone dihedral angles given residue name.
	 * @param name name of the residue
	 * @param name_next name of the succeeding residue (for pre-Pro residues)
	 * @return a pair of dihedral angles
	 */
	DihedralAngle* getRandomDihedralAngle( string name, string name_next);

	/**
	 * @brief Evaluate the probability of a residue's structure according to Ramachandran plot
	 * @param chain the chain which contains the residue
	 * @param index the index of the residue in the chain
	 * @return probability of the residue structure
	 */
	double getResidueAngleProbability( PChain* chain, int index);

	static const int PRO = 0;
	static const int PRE_PRO = 1;
	static const int GLY = 2;
	static const int GENERIC = 3;

private:
	//accumulative probabilities
	double* pro_acc;
	double* pre_pro_acc;
	double* gly_acc;
	double* gen_acc;
	//pure probabilities
	double* pro;
	double* pre_pro;
	double* gly;
	double* gen;
	int size;
	int gridnum;
	double gridlength;

	/**
	 * @brief Construct the Ramachandran plot
	 * @param filename data file
	 * @return a database
	 */
	double* construct( char* filename);
};

#endif /* RAMACHANDRANPLOT_H_ */
