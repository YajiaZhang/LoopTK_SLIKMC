/*
 * RamachandranPlot.cc
 *
 *  Created on: Jun 6, 2012
 *      Author: Yajia
 */

#include "RamachandranPlot.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "PConstants.h"

using namespace std;

RamachandranPlot::RamachandranPlot(int gridnum, char* file_pro, char* file_pre_pro, char* file_gly, char* file_generic) {
	// TODO Auto-generated constructor stub
	this->size = gridnum * gridnum;
	this->gridnum = gridnum;
	this->gridlength = 360.0 / gridnum;
//	cout << "Constructing " << file_pro << endl;
	this->pro = this->construct( file_pro);

//	cout << "Constructing " << file_pre_pro << endl;
	this->pre_pro = this->construct( file_pre_pro);

//	cout << "Constructing " << file_gly << endl;
	this->gly = this->construct( file_gly);

//	cout << "Constructing " << file_generic << endl;
	this->gen = this->construct( file_generic);

	this->pro_acc = new double[ size + 1];
	this->pre_pro_acc = new double[ size + 1];
	this->gly_acc = new double[ size + 1];
	this->gen_acc = new double[ size + 1];

	double sum_pro = 0; this->pro_acc[0] = 0;
	double sum_pre_pro = 0; this->pre_pro_acc[0] = 0;
	double sum_gly = 0; this->gly_acc[0] = 0;
	double sum_gen = 0; this->gen_acc[0] = 0;

	for( int i = 0; i < size; i++) {
		sum_pro += this->pro[i];
		sum_pre_pro += this->pre_pro[i];
		sum_gly += this->gly[i];
		sum_gen += this->gen[i];

		this->pro_acc[i + 1] = sum_pro;
		this->pre_pro_acc[ i + 1] = sum_pre_pro;
		this->gly_acc[ i + 1] = sum_gly;
		this->gen_acc[ i + 1] = sum_gen;

		if( this->pro[i] == 0.0) {
			this->pro[i] = 0.000000000001;
		}
		if( this->pre_pro[i] == 0.0) {
			this->pre_pro[i] = 0.000000000001;
		}
		if( this->gly[i] == 0.0) {
			this->gly[i] = 0.000000000001;
		}
		if( this->gen[i] == 0.0) {
			this->gen[i] = 0.000000000001;
		}
	}
	return;
}

RamachandranPlot::~RamachandranPlot() {
	if( this->gen_acc != NULL)
		delete this->gen_acc;
	if( this->pre_pro_acc != NULL)
		delete this->pre_pro_acc;
	if( this->gly_acc != NULL)
		delete this->gly_acc;
	if( this->pro_acc != NULL)
		delete this->pro_acc;
	return;
}

//Perform a binary search and return a pair of angles
DihedralAngle* RamachandranPlot::getRandomDihedralAngle( int type) {
	double p = (double)rand()/(double)RAND_MAX;
	assert( p >= 0 && p <= 1);

	double* data = NULL;
	switch( type)
	{
	case PRO:
		data = this->pro_acc;
		break;
	case PRE_PRO:
		data = this->pre_pro_acc;
		break;
	case GLY:
		data = this->gly_acc;
		break;
	case GENERIC:
		data = this->gen_acc;
		break;
	default:
		cout << "Error" << endl;
		exit(-1);
 	}
	/*NOTE: Actually the data size is size + 1 with an extra 0 added.
	This makes sure that data[start] <=p <= data[end] and
	therefore we can perform the binary search.
	*/
	int start = 0;
	int end = this->size;
	while (start < end - 1)
	{
		int mid = (int) ((start + end) / 2);
		if (data[mid] > p)
		{
			end = mid;
		}
		else {
			start = mid;
		}

	}
//	cout << data[start] << "\t" << data[end] << "\t" << p << endl;
	assert( start == (end - 1));
	if( !(data[start] <= p && p <= data[end]))
	{
		cout << p << endl;
		cout << data[start] << "\t" << data[end] << endl;
		cout << start << "\t" << end << endl;
		cout << type << endl;
		exit(-1);
	}

	//NOTE: start ranges from 0 to 399
	int x = (int)(start / this->gridnum);
	int y = start % this->gridnum;
	DihedralAngle* da = new DihedralAngle();
	da->phi = this->gridlength * x - 180 + ((double)rand()/(double)RAND_MAX) * this->gridlength;
	da->psi = this->gridlength * y - 180 + ((double)rand()/(double)RAND_MAX) * this->gridlength;
	return da;
}

double* RamachandranPlot::construct(char* filename) {
	double* data = new double[size];
	ifstream in;
	in.open( filename);
	if( in.is_open())
	{
		string line;
		getline( in, line);
		stringstream ss (stringstream::in | stringstream::out);
		ss << line;
		for( int i = 0; i < size; i++)
		{
			double temp = 0;
			ss >> temp;
			data[i] = temp;
		}
	}
	else
		cout << "File IO error" << endl;
	in.close();
	return data;
}

DihedralAngle* RamachandranPlot::getRandomDihedralAngle(string name, string name_next) {
	/*TODO: currently, doesn't consider pre_proline figure in other code.
	 * Therefore, only call this function by given ONE parameter.
	 */
	int type = -1;
	if( name_next == "PRO")
		type = PRE_PRO;
	else if( name == "GLY")
		type = GLY;
	else if( name == "PRO")
		type = PRO;
	else
		type = GENERIC;
	assert( type != -1);
	return this->getRandomDihedralAngle(type);
}

double RamachandranPlot::getResidueAngleProbability(PChain* chain, int index) {

	const DihedralAngle* angle_pair = chain->getDihedralAngleAtResidue(index);
	string name = chain->getResidue(index)->getName();
	string name_next = "";
	if( index < chain->size() - 1) {
		name_next = chain->getResidue( index + 1)->getName();
	}
	else {
		//Not the top level chain, the next residue is on other subchain
		PChain* chain_top = chain->getTopLevelChain();
		if( chain != chain_top) {
			int end = chain->getTopLevelIndices().second;
			if( end != chain_top->size() - 1)
				name_next = chain_top->getResidue( end + 1)->getName();
			//else it is the end of the top chain
		}
		//else, it is the end of the chain (the chain is also the top chain).
	}
	double* data = NULL;
	if( name_next == "PRO")
		data = this->pre_pro;
	else if( name == "PRO")
		data = this->pro;
	else if( name == "GLY")
		data = this->gly;
	else
		data = this->gen;

	double prob = 0.0;
	//Kind like a marginalization process
	//NOTE: The dihedral angle pair of the first residue in the top level chain.
	if( angle_pair->phi > 300.0) {
			int y = (int)((angle_pair->psi  + 180.0)/ this->gridlength);
			for( int x = 0; x < this->gridnum; x++)
			{
				int data_index = x * this->gridnum + y;
				if(!(data_index >= 0 && data_index < this->size))
				{
					cerr << "Ramachandran Plot error!" << endl;
					chain->printBackbonePosition();
					chain->printDihedralAngles();
					exit(-1);
				}
				prob += data[ data_index];
			}
		}
		//NOTE: The dihedral angle pair of the last residue in the top level chain.
		else if( angle_pair->psi > 300.0)
		{
			int x = (int)((angle_pair->phi + 180.0) / this->gridlength);
			for( int y = 0; y < this->gridnum; y++)
			{
				int data_index = x * this->gridnum + y;
				if(!(data_index >= 0 && data_index < this->size))
				{
					cerr << "Ramachandran Plot error!" << endl;
					chain->printBackbonePosition();
					chain->printDihedralAngles();
					exit(-1);
				}
				prob += data[ data_index];
			}
		}
		//NOTE: Not the first nor the last residue in the top level chain.
		else
		{
			int x = (int)((angle_pair->phi + 180.0) / this->gridlength);
			int y = (int)((angle_pair->psi  + 180.0)/ this->gridlength);

			int data_index = x * this->gridnum + y;
			if(!(data_index >= 0 && data_index < this->size))
			{
				cerr << "Ramachandran Plot error!" << endl;
				chain->printBackbonePosition();
				chain->printDihedralAngles();
				exit(-1);
			}
			prob = data[ data_index];
		}
	return prob;
}
