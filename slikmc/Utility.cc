/*
 * Utility.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: Yajia
 */

#include "Utility.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
using namespace std;

//void Utility::print_Double2D(int m, int n, double** values, int start, string name) {
//	if( name != "")
//		cout << name << ":" << endl;
//	for( int i = start; i < m + start; i++)
//	{
//		for( int j = start; j < n + start; j++)
//		{
//			cout << values[i][j] << "\t";
//		}
//		cout << endl;
//	}
//}

void Utility::matrix_multiply_trash(int m, int n, int p, double** val1, double** val2, double** res, double cov) {
	for( int i = 1; i <= m; i++)
	{
		for( int j = 1; j <=p; j++)
		{
			double temp = 0;
			for( int k = 1; k <=n; k++)
			{
				temp += val1[i][k] * val2[k][j];
			}
			res[i][j] = temp * cov;
		}
	}
	return;
}

void Utility::matrix_square_transposeA(int m, int n, double** A, double** res) {
	//result should be n * n
	for( int i = 1; i <= n; i++)
	{
		for( int j = 1; j <=n; j++)
		{
			double temp = 0;
			for( int k = 1; k <= m; k++)
			{
				temp += A[k][i] * A[k][j];
			}
			res[i][j] = temp;
		}
	}
	return;
}

double** Utility::new_Double2D(int m, int n) {
	double** values = new double*[m];
	for( int i = 0; i < m; i++)
	{
		values[i] = new double[n];
	}
	return values;
}

void Utility::delete_Double2D(double** values, int m) {
	for( int i = 0; i < m; i++)
	{
		delete[] values[i];
	}
	delete[] values;
}

void Utility::setZero_Double2D(int m, int n, double** values) {
	for( int i = 0; i < m; i++)
	{
		for( int j = 0; j < n; j++)
			values[i][j] = 0;
	}
	return;
}

//void Utility::print_vector(vector<double>& v, string delim) {
//	for( int i = 0; i < v.size(); i++)
//	{
//		cout << v[i] << delim;
//	}
//
//	cout << endl;
//}

StringTokenizer::StringTokenizer( string s, string delims)
{
	this->s = s;
	for( int i = 0; i < delims.length(); i++)
	{
		this->delims.push_back( delims[i]);
	}
	this->index = 0;
	return;
}

string StringTokenizer::nextToken()
{
	string s = "";
//	bool end = false;
	for(; this->index < this->s.size(); this->index++)
	{
		char temp = this->s[ this->index];
		bool isIn = false;
		for( int j = 0; j < this->delims.size(); j++)
		{
			if( temp == this->delims[j])
			{
				isIn = true;
				break;
			}
		}
		if( isIn == true && s.size() == 0)
			continue;
		else if( isIn == true && s.size() > 0)
			break;
		else
		{
			s += temp;
			continue;
		}
	}

	if( s == "")
	{
		return NULL;
	}
	else
		return s;
}

bool StringTokenizer::hasMoreTokens() const {
	for( int i = this->index; i < this->s.size(); i++)
	{
		char temp = this->s[i];
		bool isIn = false;
		for( int j = 0; j < this->delims.size(); j++)
		{
			if( temp == this->delims[j])
			{
				isIn = true;
				break;
			}
		}
		if( isIn == false)
			return true;
	}
	return false;
}

//void Utility::record_vector(vector<double>& v, string filename, string delim) {
//	ofstream out;
//	out.open( filename.c_str());
//	if( !out.is_open())
//	{
//		cerr << "File open error!" << endl;
//		exit(-1);
//	}
//
//	for( int i = 0; i < v.size(); i++)
//	{
//		out << v[i] << delim;
//	}
//
//	out.flush();
//	out.close();
//	return;
//}

void StringTokenizer::getAllTokens(vector<string>& tokens) {
	while( this->hasMoreTokens()) {
		tokens.push_back( this->nextToken());
	}
	return;
}

string String::toUpper(const string& temp) {
	string s = "";
	for( int i = 0; i < temp.length(); i++) {
		s += std::toupper( temp[i]);
	}
	return s;
}

string String::toLower(const string& temp) {
	string s = "";
	for( int i = 0; i < temp.length(); i++) {
		s += std::tolower( temp[i]);
	}
	return s;
}

//string Utility::itoa(const int& i) {
//	   stringstream ss;//create a stringstream
//	   ss << i;//add number to the stream
//	   return ss.str();//return a string with the contents of the stream
//}

double Utility::dist(const vector<double>& a, const vector<double>& b, const int minsize) {
	assert( minsize > 0);
	if( !(a.size() >= minsize && b.size() >= minsize)) {
		cerr << "a.size:" << a.size() << endl;
		cout << "b.size:" << b.size() << endl;
		cout << "minsize:" << minsize << endl;
		abort();
	}

	double sum = 0;
	for( int i = 0; i < minsize; i++) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt( sum);
}

void Debug::check(const bool assertion, string comment) {
	if( assertion == false) {
		cout << comment << endl;
		abort();
	}
}

int Random::nextInt(const int range) {
	assert( range > 0);
	return rand() % range;
}

int Random::nextInt(const int r1, const int r2) {
	assert( r1 < r2);
	return rand() % ( r2 - r1 + 1) + r1;
}

double Random::nextDouble(const double range) {
	assert( range > 0);
	return ( rand() / ((double)INT_MAX)) * range;
}

double Random::nextDouble(const double r1, const double r2) {
	if( r1 >= r2)
	{
//		cout << "Random::nextDouble - warning: the lower bound is larger than or equal to upper bound!" << endl;
		return r1;
	}
	return nextDouble( r2 - r1) + r1;
}

//Box¨CMuller transform
double Random::nextNormal() {
	double U = Random::nextDouble(1);
	double V = Random::nextDouble(1);

	double z1 = sqrt( -2 * log(U)) * cos( 2 * PI * V);
	double z2 = sqrt( -2 * log(U)) * sin( 2 * PI * V);
	return U > 0.5? z1: z2;
}

double Random::nextNormal(const double mean, const double std) {
	return std * Random::nextNormal() + mean;
}

//Assume the dataset is in increasing order;
//also d must satisfy data.front <= d <= data.back
//return i if data[i-1] < d <= data[i]
int BinarySearch::search(const vector<double>& data, const double& d) {
//	assert( d >= data[0] && d <= data.back());
	if( !(d >= data[0] && d <= data.back())) {
		cout << "d0:" << data[0] << endl;
		cout << "dback:" << data.back() << endl;
		cout << "d:" << d << endl;
		abort();
	}
	return BinarySearch::search( data, d, 0, data.size()-1);
}

int BinarySearch::search( const vector<double>& data, const double& d, const int s, const int e) {

	if( e <= (s + 1))
		return e;

	int index = (s + e) / 2;
	if( data[index] >= d) {
		return BinarySearch::search( data, d, s, index);
	}
	else {
		return BinarySearch::search( data, d, index, e);
	}
}
