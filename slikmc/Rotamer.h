/*
 * Rotamer.h
 *
 *  Created on: Jan 24, 2013
 *      Author: Yajia
 */

#ifndef ROTAMER_H_
#define ROTAMER_H_
#include <vector>
#include <string.h>
#include <iostream>
#include <map>
#include <PProtein.h>
using namespace std;

class ChiDistribution;
class RotamerGrid;
class RotamerGridSpecial;
class RotamerGridKey;
class Rotamer;

/**
 * @brief An auxiliary class for class Rotamer. A Comparator to compare two RotamerGridKeys.
 */
struct RotamerGridKeyComparator {
     bool operator()(RotamerGridKey lhv, RotamerGridKey rhv) const;
};

typedef map<RotamerGridKey, RotamerGrid, RotamerGridKeyComparator> GridMap;
typedef map<RotamerGridKey, RotamerGridSpecial, RotamerGridKeyComparator> GridMapSpecial;
typedef map<RotamerGridKey, RotamerGrid, RotamerGridKeyComparator>::iterator GridMapIter;
typedef map<RotamerGridKey, RotamerGridSpecial, RotamerGridKeyComparator>::iterator GridMapIterSpecial;

/**
 * @brief An auxiliary class for class Rotamer. This class is a data structure for storing a specific rotamer in one dihedral angle grid.
 */
class ChiDistribution {
public:
	double prob;
	vector< double> means;
	vector< double> stds;
};

/**
 * @brief An auxilliary class for class Rotamer. This class is a data structure for storing all rotamers in one dihedral angle grid.
 */
class RotamerGrid{
	friend class Rotamer;
public:
	//number of effective chi angles
	int n_chi;
	vector<double> acc_rotamer;
	//accumulated distribution for each rotamer
	vector<ChiDistribution> dist_rotamer;
};

/**
 * @brief An auxilliary class for class Rotamer. This class extends RotamerGrid and it adds additional structures for storing non-rotameric residues {ASN, ASP, PHE, TRP, HIS, TYR, GLN, GLU}.
 */
class RotamerGridSpecial: public RotamerGrid {
	friend class Rotamer;
public:
	int degree_begin;
	int degree_end;
	int stepLength;
	bool check() const;
	vector<double> dist_terminal;
	vector< vector<double> > acc_terminal;
	vector< vector<double> > chi_terminal;
};

/**
 * @brief An auxilliary class for class Rotamer. This class used as a key to access one RotamerGrid in the map structure.
 */
class RotamerGridKey {
	friend class Rotamer;
	friend struct RotamerGridKeyComparator;
public:
	/**
	 * @brief Constructor
	 * @param name residue name
	 * @param phi backbone phi angle in 10 base
	 * @param psi backebone psi angle in 10 base
	 */
	RotamerGridKey( const string& name, const int& phi, const int& psi);

	/**
	 * @brief Get a string containing residue name, backbone phi and psi angles.
	 * @return a string
	 */
	string toString() const;

	/**
	 * @brief Print the residue name, backbone phi,psi angles.
	 */
	void print() const;

	/**
	 * @brief "==" operator override to compare two RotemerGridKeys.
	 * @param key an instance of RotamerGridKey
	 * @return true, same; otherwise, false
	 */
    bool operator==(const RotamerGridKey& key) const;
private:
	string name;
	int phi;
	int psi;
};

/**
 * @brief The class Rotamer is used to initialize side-chain database by parsing Dunbrack's Backbone Dependent Rotamer Library,
 * sample one side-chain conformation for a given residue, evaluate the probability given one side-chain conformation.
 */
class Rotamer {
public:
	/**
	 * @Constructor
	 * @param protein protein chain to be sampled
	 */
	Rotamer(PProtein* protein);

	/**
	 * @brief Randomly sample side-chain chi angles given residue name and backbone dihedral angle pair
	 * @param residue_name name of the residue
	 * @param phi backbone phi angle
	 * @param psi backbone psi angle
	 * @param rotamer_grid stores the grid index of rotamer in the rotamer database (for evaluation purpose)
	 * @param dAngles stores the side-chain chi angles
	 */
	void sample( const string& residue_name, const double& phi, const double& psi, int& rotamer_grid, vector<double>& dAngles);

	/**
	 * @brief Randomly sample side-chain chi angles given rotameric residue name and backbone dihedral angle pair
	 * @param residue_name name of the residue
	 * @param phi backbone phi angle
	 * @param psi backbone psi angle
	 * @param rotamer_grid stores the grid index of rotamer in the rotamer database (for evaluation purpose)
	 * @param dAngles stores the side-chain chi angles
	 */
	void sample_common( const string& residue_name, const int& phi, const int& psi, int& rotamer_grid, vector<double>& dAngles);

	/**
	 * @brief Randomly sample side-chain chi angles given non-rotameric residue name and backbone dihedral angle pair
	 * @param residue_name name of the residue
	 * @param phi backbone phi angle
	 * @param psi backbone psi angle
	 * @param rotamer_grid stores the grid index of rotamer in the rotamer database (for evaluation purpose)
	 * @param dAngles stores the side-chain chi angles
	 */
	void sample_Special( const string& residue_name, const int& phi, const int& psi, int& rotamer_grid, vector<double>& dAngles);

	/**
	 * @brief Output the database to a file
	 * @param filename output filename
	 */
	void output( char* filename);

	//Evaluate the conformation of sidechain and return log(prob);
	/**
	 * @brief Evaluate side-chain structures given one protein chain
	 * @param chain protein chain to be evaluated
	 * @param s index of starting residue in protein chain
	 * @param e index of ending residue in protein chain
	 * @return probability density in logarithm
	 */
	double evalSidechain_log(PChain* chain, const int s, const int e);
//	double evalResidue(const string& residue_name, const double& phi, const double& psi, const int& rotamer_index, const vector<double>& angles);

	/**
	 * @brief Initialize sidechain database according to residues the protein contains
	 * @param protein protein chain to be sampled
	 */
	void initSidechainDatabase( PProtein* protein);
private:
	void init();
	/**
	 * @brief Read sidechain library given residue name and construct corresponding database
	 * @param res_name residue name
	 */
	void read( const string& res_name);

	/**
	 * @brief Read sidechain library given the name of a rotameric residue and construct corresponding database
	 * @param res_name residue name
	 */
	void read_common( const string& res_name);

	/**
	 * @brief Read sidechain library given the name of a non-rotameric residue and construct corresponding database
	 * @param res_name residue name
	 */
	void read_special( const string& res_name);

	GridMap gridMap;
	GridMapSpecial gridMapSpecial;

	//For special amino acids.
	vector<int> degree_start;
	vector<int> degree_end;
	vector<int> stepLength;

	static const int ASN = 0;
	static const int ASP = 1;
	static const int PHE = 2;
	static const int TRP = 3;
	static const int HIS = 4;
	static const int TYR = 5;
	static const int GLN = 6;
	static const int GLU = 7;

	/**
	 * @brief Return the type of residue given residue name
	 * @return -1: rotameric residue; 0-7: ASN-GLU
	 */
	int type( const string& res_name);

	/**
	 * @brief Check if the residue is rotameric or not
	 * @param res_name residue name
	 * @return true, rotemeric; false, non-rotameric
	 */
	bool isSpecial( const string& res_name);
	map<string, int> typeMap;
};

#endif /* ROTAMER_H_ */
