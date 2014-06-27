/*
 * Rotamer.cc
 *
 *  Created on: Jan 24, 2013
 *      Author: Yajia
 */

#include "Rotamer.h"
#include <fstream>
#include <iostream>
#include "Utility.h"
#include <PResources.h>
#include <PExtension.h>

using namespace std;

Rotamer::Rotamer(PProtein* protein) {
	this->init();
	//Since the reading process is slow, just read necessary amino acids
	map<string, bool> mapTemp;
	for( int i = 0; i < protein->size(); i++) {
		string res_name = protein->getResidue(i)->getName();
		map<string, bool>::iterator mapTempIter = mapTemp.find( res_name);
		if( mapTempIter == mapTemp.end()) {
			cout << "Read new residue " << res_name;
			this->read( res_name);
			cout << " done." << endl;
			mapTemp.insert( std::make_pair( res_name, true));
		}
		else {
			cout << res_name << " is already in the lib" << endl;
		}
	}
	return;
}

void Rotamer::sample(const string& residue_name, const double& phi_d, const double& psi_d, int& rotamerIndex, vector<double>& dAngles) {
	//Process the dihedral angles to find the grid point
	int phi = (int)(phi_d / 10.0) * 10;
	int psi = (int)(psi_d / 10.0) * 10;

	bool isSpecial = this->isSpecial( residue_name);
	if( !isSpecial) {
		this->sample_common( residue_name, phi, psi, rotamerIndex, dAngles);
	}
	else {
		this->sample_Special( residue_name, phi, psi, rotamerIndex, dAngles);
	}
	return;
}

void Rotamer::output(char* filename) {
	ofstream out;
	out.open( filename);
	for(GridMapIter gridIter = this->gridMap.begin(); gridIter != this->gridMap.end(); gridIter++) {
		const RotamerGridKey* key = &gridIter->first;
		RotamerGrid* grid = &gridIter->second;

		out << key->name << "\t" << key->phi << "\t" << key->psi << endl;
		for (int i = 0; i < grid->dist_rotamer.size(); i++) {
			out << grid->dist_rotamer[i].prob << "\t";
			for (int j = 0; j < grid->n_chi; j++) {
				out << grid->dist_rotamer[i].means[j] << "\t";
			}
			for (int j = 0; j < grid->n_chi; j++) {
				out << grid->dist_rotamer[i].stds[j] << "\t";
			}
			out << endl;
		}
	}

	for(GridMapIterSpecial gridIter = this->gridMapSpecial.begin(); gridIter != this->gridMapSpecial.end(); gridIter++) {
		const RotamerGridKey* key = &gridIter->first;
		RotamerGridSpecial* grid = &gridIter->second;

		out << key->name << "\t" << key->phi << "\t" << key->psi << endl;
		for( int i = 0; i < grid->dist_rotamer.size(); i++) {
			out << grid->dist_rotamer[i].prob << "\t";
			for (int j = 0; j < grid->n_chi - 1; j++) {
				out << grid->dist_rotamer[i].means[j] << "\t";
			}
			for (int j = 0; j < grid->n_chi - 1; j++) {
				out << grid->dist_rotamer[i].stds[j] << "\t";
			}
			for( int j = 0; j < grid->chi_terminal[i].size(); j++) {
				out << grid->chi_terminal[i][j] << "\t";
			}
			out << endl;
		}
	}
	out.flush();
	out.close();
	return;
}

double Rotamer::evalSidechain_log(PChain* chain, const int start, const int end) {
	vector<DihedralAngle> da_backbone;
	chain->getDihedralAngles( da_backbone);
	double log_prob = 0;
	for( int i = start; i <= end; i++) {
		PResidue* residue = chain->getResidue(i);
		string name = residue->getName();
		if( name == "ALA" || name == "GLY") continue;
		int phi = (int)(da_backbone[i].phi / 10) * 10;
		int psi = (int)(da_backbone[i].psi / 10) * 10;

		int index = 0;
		vector<double> da_sidechain;
		residue->getSideChainAngles( index, da_sidechain);
		int res_index = this->type( name);
		if( res_index < 0){
			GridMapIter gridIter = this->gridMap.find( RotamerGridKey(name, phi, psi));
			if (gridIter == this->gridMap.end()) {
				cout << "Cannot find the corresponding key" << endl;
				cout << name << "\t" << phi << "\t" << psi << endl;
				abort();
			}
			const RotamerGrid* grid = &gridIter->second;
			log_prob += log( grid->dist_rotamer[index].prob);
			}
		else {
			GridMapIterSpecial gridIter = this->gridMapSpecial.find( RotamerGridKey(name, phi, psi));
			if (gridIter == this->gridMapSpecial.end()) {
				cout << "Cannot find the corresponding key" << endl;
				cout << name << "\t" << phi << "\t" << psi << endl;
				abort();
			}
			const RotamerGridSpecial* grid = &gridIter->second;
			log_prob += log( grid->dist_rotamer[index].prob);
			int terminalIndex = (int)((da_sidechain.back() - this->degree_start[res_index]) / this->stepLength[res_index]);
			log_prob += log( grid->chi_terminal[index][terminalIndex]);
		}
	}
	return log_prob;
}

void Rotamer::read( const string& res_name_U) {
	string res_name = String::toLower( res_name_U);
	if( res_name == "ala" || res_name == "gly") {
		cout << res_name << " does not have a rotamer! Skip reading" << endl;
		return;
	}
	if( this->isSpecial( res_name_U)) {
		this->read_special( res_name);
	}
	else {
		this->read_common( res_name);
	}
	return;
}

void Rotamer::read_common( const string& res_name) {

	string filename = "../Data/ExtendedOpt1-5/" + res_name + ".bbdep.rotamers.lib";
	ifstream in;
	in.open( filename.c_str());
	if( !in.is_open()) {
		cerr << "File Open error:" << filename << endl;
		abort();
	}

	int line_index = 1;
	while( true) {
		string line;
		getline( in, line);
		if( line == "") {
			break;
		}

		//Skip comments
		if( line.c_str()[0] == '#') {
			continue;
		}
		vector<string> tokens;
		StringTokenizer stn( line);
		while( stn.hasMoreTokens()) {
			tokens.push_back( stn.nextToken());
		}
		assert( tokens.size() == 17);

		RotamerGridKey key( tokens[0], atoi( tokens[1].c_str()), atoi( tokens[2].c_str()));

		GridMapIter iter_grid = this->gridMap.find(key);
		//already in the map
		if( iter_grid != this->gridMap.end()) {
			RotamerGrid* grid = &iter_grid->second;
			double acc = atof( tokens[8].c_str());
			grid->acc_rotamer.push_back( grid->acc_rotamer.back() + acc);

			ChiDistribution chis;
			for( int i = 0; i < grid->n_chi; i++) {
				chis.means.push_back( atof( tokens[9 + i].c_str()));
				chis.stds.push_back( atof( tokens[13 + i].c_str()));
			}
			chis.prob = atof( tokens[8].c_str());

			grid->dist_rotamer.push_back( chis);


		}
		else {
			RotamerGrid grid;
			grid.n_chi = 0;
			//count how many chi angles to represent the side-chain conformation.
			for( int i = 4; i <= 7; i++) {
				if( tokens[i] != "0") {
					grid.n_chi += 1;
				}
				else {
					for( ; i <=7; i++)
						assert( tokens[i] == "0");
				}
			}
			grid.acc_rotamer.push_back(0);
			grid.acc_rotamer.push_back( atof( tokens[8].c_str()));

			ChiDistribution chis;
			for( int i = 0; i < grid.n_chi; i++) {
				chis.means.push_back( atof( tokens[9 + i].c_str()));
				chis.stds.push_back( atof( tokens[13 + i].c_str()));
			}
			chis.prob = atof( tokens[8].c_str());

			grid.dist_rotamer.push_back( chis);
			this->gridMap.insert( std::pair< RotamerGridKey, RotamerGrid>( key, grid));
		}
	}
	in.close();
	return;
}

void Rotamer::read_special( const string& res_name) {
	string libname = "../Data/ExtendedOpt1-5/" + res_name + ".bbdep.densities.lib";
	ifstream libin;
	libin.open( libname.c_str());
	if( !libin.is_open()) {
		cerr << "File Open error:" << libname << endl;
		abort();
	}

	int index_t = 8;
	int index_std = 7;
	int index_prob = 5;
	int index_mean = 6;

	int res_index = this->type( String::toUpper( res_name));
	if( res_index > 5) {
		index_t = 11;
		index_std = 9;
		index_prob = 6;
		index_mean = 7;
	}
	int size_terminal = (this->degree_end[res_index] - this->degree_start[res_index]) / this->stepLength[res_index] + 1;

	int line_index = 1;
	while( true) {
		string line; getline( libin, line);
		if( line == "") {
			break;
		}
		//Skip comments
		if( line.c_str()[0] == '#') {
			continue;
		}
		StringTokenizer stn( line);
		vector<string> tokens; stn.getAllTokens( tokens);

		//calculate the size for one record
		int size_tokens = index_t + size_terminal;
		assert( tokens.size() == size_tokens);

		RotamerGridKey key( tokens[0], atoi( tokens[1].c_str()), atoi( tokens[2].c_str()));

		GridMapIterSpecial iter_grid = this->gridMapSpecial.find( key);
		//already in the map
		if( iter_grid != this->gridMapSpecial.end()) {
			RotamerGridSpecial* grid = &iter_grid->second;

			ChiDistribution chis;
			for( int i = 0; i < grid->n_chi - 1; i++) {
				chis.means.push_back( atof( tokens[index_mean + i].c_str()));
				chis.stds.push_back( atof( tokens[index_std + i].c_str()));
			}
			chis.prob = atof( tokens[index_prob].c_str());
			grid->acc_rotamer.push_back( grid->acc_rotamer.back() + chis.prob);
			grid->dist_rotamer.push_back( chis);

			vector<double> chi2_dist;
			vector<double> chi2_dist_acc; chi2_dist_acc.push_back(0);
			for( int i = index_t; i < size_tokens; i++) {
				chi2_dist.push_back( atof( tokens[i].c_str()));
				chi2_dist_acc.push_back( chi2_dist_acc.back() + chi2_dist.back());
			}
			chi2_dist_acc[ chi2_dist_acc.size() - 1] = 1; //NOTE: Awkward method, atof lost precision 10(-4);
			grid->chi_terminal.push_back( chi2_dist);
			grid->acc_terminal.push_back( chi2_dist_acc);
		}
		else {
			RotamerGridSpecial grid;
			grid.n_chi = (index_prob == 5 ? 2: 3);

			grid.degree_begin = this->degree_start[res_index];
			grid.degree_end = this->degree_end[res_index];
			grid.stepLength = this->stepLength[res_index];

			ChiDistribution chis;
			for( int i = 0; i < grid.n_chi - 1; i++) {
				chis.means.push_back( atof( tokens[index_mean + i].c_str()));
				chis.stds.push_back( atof( tokens[index_std + i].c_str()));
			}
			chis.prob = atof( tokens[index_prob].c_str());
			grid.acc_rotamer.push_back(0);
			grid.acc_rotamer.push_back( chis.prob);
			grid.dist_rotamer.push_back( chis);

			vector<double> chi2_dist_acc; chi2_dist_acc.push_back(0);
			vector<double> chi2_dist;

			for( int i = index_t; i < size_tokens; i++) {
				chi2_dist.push_back( atof( tokens[i].c_str()));
				grid.dist_terminal.push_back( this->degree_start[res_index] + this->stepLength[res_index] * (i - index_t));
				chi2_dist_acc.push_back( chi2_dist.back() + chi2_dist_acc.back());
			}
//			NOTE:
			chi2_dist_acc[ chi2_dist_acc.size() - 1] = 1; //NOTE: Awkward method, atof lost precision
			grid.chi_terminal.push_back( chi2_dist);
			grid.acc_terminal.push_back( chi2_dist_acc);

			this->gridMapSpecial.insert( std::pair< RotamerGridKey, RotamerGridSpecial>( key, grid));
		}
	}
	libin.close();
	return;
}

RotamerGridKey::RotamerGridKey(const string& name, const int& phi, const int& psi) {
	this->name = name;
	this->phi = phi;
	this->psi = psi;
	return;
}

string RotamerGridKey::toString() const {
	string s = "";
	s += this->name;
	s += "\t";
	s += this->phi;
	s += "\t";
	s += this->psi;
	return s;
}

void RotamerGridKey::print() const {
	cout << this->name << "\t" << this->phi << "\t" << this->psi;
	return;
}

bool RotamerGridKey::operator ==(const RotamerGridKey& rhv) const {
	if( (this->name == rhv.name) && ( this->phi == rhv.phi) && ( this->psi == rhv.psi)) {
		return true;
	}
	return false;
}


bool RotamerGridKeyComparator::operator ()(RotamerGridKey k1, RotamerGridKey k2) const {
	if( k1.name < k2.name) {
		return true;
	}
	else if( k1.phi < k2.phi && k1.name == k2.name) {
		return true;
	}
	else if( k1.psi < k2.psi && k1.name == k2.name && k1.phi == k2.phi) {
		return true;
	}
	return false;
}

void Rotamer::sample_common(const string& residue_name, const int& phi, const int& psi, int& rotamerIndex, vector<double>& dAngles) {
		RotamerGridKey key( residue_name, phi, psi);
		GridMapIter iter = this->gridMap.find( key);
		if( iter == this->gridMap.end()) {
			cerr << "No matched item? There must be something wrong! sample_common" << endl;
			cout << "key:" << key.name << "\t" << key.phi << "\t" << key.psi << endl;
			abort();
		}

		RotamerGrid* dist = &iter->second;

		double prob = Random::nextDouble( dist->acc_rotamer.back()) - 0.000000001; if( prob < 0) prob = 0.000000001;
		int index = BinarySearch::search( dist->acc_rotamer, prob) - 1;

		assert( dAngles.size() == 0);
		if (!( index <= dist->dist_rotamer.size() && index >= 0)) {
			cout << "i_max:" << dist->dist_rotamer.size()<< endl;;
			cout << "i:" << index << endl;
			cout << "prob:" << prob << endl;
			for( int i = 0; i < dist->acc_rotamer.size(); i++) {
				cout << dist->acc_rotamer[i] << "\t";
			}
			cout << endl;
			abort();
		}
		const ChiDistribution* rotamer_dist = &dist->dist_rotamer[index];
		cout << endl;
		for( int i = 0; i < dist->n_chi; i++) {
			double mean = rotamer_dist->means[i];
			double std = rotamer_dist->stds[i];
			dAngles.push_back( Random::nextNormal( mean, std));
		}
		rotamerIndex = index;
		return;
}

void Rotamer::sample_Special(const string& residue_name, const int& phi, const int& psi, int& rotamerIndex, vector<double>& dAngles) {
	RotamerGridKey key( residue_name, phi, psi);
	GridMapIterSpecial iter = this->gridMapSpecial.find( key);
	if( iter == this->gridMapSpecial.end()) {
		cerr << "No matched item? There must be something wrong! sample_special" << endl;
		cout << "key:" << key.name << "\t" << key.phi << "\t" << key.psi << endl;
		abort();
	}
	RotamerGridSpecial* dist = &iter->second;

	double prob_rotamer = Random::nextDouble( dist->acc_rotamer.back()) - 0.00001; if( prob_rotamer < 0) prob_rotamer = 0.00001;
	int index = BinarySearch::search( dist->acc_rotamer, prob_rotamer) - 1;

	assert( dAngles.size() == 0);

	const ChiDistribution* rotamer_dist = &dist->dist_rotamer[index];
	for( int i = 0; i < dist->n_chi - 1; i++) {
		double mean = rotamer_dist->means[i];
		double std = rotamer_dist->stds[i];
		dAngles.push_back( Random::nextNormal( mean, std));
	}
	//handle the terminal chi
	double prob_terminal = Random::nextDouble(1) - 0.000001; if( prob_terminal < 0) prob_terminal = 0.00001;
	int index_terminal = BinarySearch::search( dist->acc_terminal[index], prob_terminal) - 1;
	double chi_terminal = dist->degree_begin + index_terminal * dist->stepLength;
	dAngles.push_back( chi_terminal);

	rotamerIndex = index;
}

int Rotamer::type(const string& res_name) {
	map<string, int>::iterator iter = this->typeMap.find( res_name);
	if( iter != this->typeMap.end()) {
		return iter->second;
	}
	cout << res_name << endl;
	cout << "no luck?" << endl;
	abort();
	return -2;
}

void Rotamer::init() {
	this->degree_start.resize( 8);
	this->degree_end.resize( 8);
	this->stepLength.resize(8);
	//Initialize the starting degree and ending degree for all 8 special amino acids
	this->degree_start[ASN] = -180;
	this->degree_start[ASP] = -90;
	this->degree_start[PHE] = -30;
	this->degree_start[TRP] = -180;
	this->degree_start[HIS] = -180;
	this->degree_start[TYR] = -30;
	this->degree_start[GLN] = -180;
	this->degree_start[GLU] = -90;

	this->degree_end[ASN] = 170;
	this->degree_end[ASP] = 85;
	this->degree_end[PHE] = 145;
	this->degree_end[TRP] = 170;
	this->degree_end[HIS] = 170;
	this->degree_end[TYR] = 145;
	this->degree_end[GLN] = 170;
	this->degree_end[GLU] = 85;

	this->stepLength[ASN] = 10;
	this->stepLength[ASP] = 5;
	this->stepLength[PHE] = 5;
	this->stepLength[TRP] = 10;
	this->stepLength[HIS] = 10;
	this->stepLength[TYR] = 5;
	this->stepLength[GLN] = 10;
	this->stepLength[GLU] = 5;

	//Amino acids with rotamers
	this->typeMap.insert( std::pair< string, int>( "ARG", -1));
	this->typeMap.insert( std::pair< string, int>( "CYS", -1));
	this->typeMap.insert( std::pair< string, int>( "ILE", -1));
	this->typeMap.insert( std::pair< string, int>( "LEU", -1));
	this->typeMap.insert( std::pair< string, int>( "LYS", -1));
	this->typeMap.insert( std::pair< string, int>( "MET", -1));
	this->typeMap.insert( std::pair< string, int>( "PRO", -1));
	this->typeMap.insert( std::pair< string, int>( "SER", -1));
	this->typeMap.insert( std::pair< string, int>( "THR", -1));
	this->typeMap.insert( std::pair< string, int>( "VAL", -1));
	//Non-rotameric amino acids with 2 chis
	this->typeMap.insert( std::pair< string, int>( "ASN", 0));
	this->typeMap.insert( std::pair< string, int>( "ASP", 1));
	this->typeMap.insert( std::pair< string, int>( "PHE", 2));
	this->typeMap.insert( std::pair< string, int>( "TRP", 3));
	this->typeMap.insert( std::pair< string, int>( "HIS", 4));
	this->typeMap.insert( std::pair< string, int>( "TYR", 5));
	//Non-rotameric amino acids with 3 chis
	this->typeMap.insert( std::pair< string, int>( "GLN", 6));
	this->typeMap.insert( std::pair< string, int>( "GLU", 7));
}

void Rotamer::initSidechainDatabase(PProtein* protein) {

	assert( protein->getTopLevelChain() == protein);
	int start = 1;
	int end = protein->size() - 2;

	vector<DihedralAngle> da_backbone;
	protein->getDihedralAngles( da_backbone);
	for( int j = start; j <= end; j++) {
		PResidue* residue = protein->getResidue(j);
		string name = residue->getName();
		if( name == "GLY" || name == "ALA") continue;

		unsigned chiMax = PResources::numChiIndices( residue->getName());

		residue->angles_sidechain.resize( chiMax);
		vector<string> rotList;
		vector<double> curChi;
		//apply rotamer angles
		for (unsigned i = 1; i <= chiMax; i++) {
			//get the atoms used to define dihedral angles
			rotList = PResources::GetChiIndex( residue->getName(), i);
			//calculate current dihedral angle (chi angle) and calculate the delta needed to achive the goal.
			double curDihedral = PMath::AngleBetweenPlanes(
					residue->getAtomPosition(rotList[0]),
					residue->getAtomPosition(rotList[1]),
					residue->getAtomPosition(rotList[2]),
					residue->getAtomPosition(rotList[3]));
			curChi.push_back( curDihedral);
		}
		int phi = (int)(da_backbone[j].phi / 10) * 10;
		int psi = (int)(da_backbone[j].psi / 10) * 10;
		RotamerGridKey key( name, phi, psi);
		if( !this->isSpecial( name)) {
			GridMapIter iter = this->gridMap.find( key);
			assert( iter != this->gridMap.end());

			double error = numeric_limits<double>::max( );
			int rotamerIndex = -1;
			for( int k = 0; k < iter->second.dist_rotamer.size(); k++) {
				double temp = Utility::dist(curChi, iter->second.dist_rotamer[k].means, chiMax);
				if( error > temp) {
					error = temp;
					rotamerIndex = k;
				}
			}
			residue->type_sidechain = rotamerIndex;
		}
		else {
			GridMapIterSpecial iter = this->gridMapSpecial.find( key);
			assert( iter != this->gridMapSpecial.end());
			double error = numeric_limits<double>::max( );
			int rotamerIndex = -1;
			for( int k = 0; k < iter->second.dist_rotamer.size(); k++) {
				double temp = Utility::dist(curChi, iter->second.dist_rotamer[k].means,  chiMax - 1);
				if( error > temp) {
					error = temp;
					rotamerIndex = k;
				}
			}
			residue->type_sidechain = rotamerIndex;
		}
	}

}

bool Rotamer::isSpecial(const string& res_name) {
	int type = this->type( res_name);
	if( type == -1)
		return false;
	else
		return true;
}
