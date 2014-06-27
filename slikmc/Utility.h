/*
 * Utility.h
 *
 *  Created on: Jul 11, 2012
 *      Author: Yajia
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <string.h>
#include <vector.h>
using namespace std;

#define PI 3.14159265

/**
 * @brief An auxiliary class containing vector and matrix operations.
 */
class Utility {
public:

//	static void print_Double2D( int m, int n, double** values, int start = 0, string name = "");

	/**
	 * @brief
	 * @param m
	 * @param n
	 * @param p
	 * @param val1
	 * @param val2
	 * @param res
	 * @param cov
	 */
	static void matrix_multiply_trash( int m, int n, int p, double** val1, double** val2, double** res, double cov = 1);

	static void matrix_square_transposeA( int m, int n, double** A, double** res);

	static double** new_Double2D( int m, int n);

	static void delete_Double2D( double** values, int m);

	static void setZero_Double2D( int m, int n, double** values);

//	static void print_vector( vector<double>& v, string delim = " ");

//	static void record_vector( vector<double>& v, string filename, string delim = "\n");

//	static string itoa( const int& i);

	static double dist( const vector<double>& a, const vector<double>& b, const int minsize);

};

/**
 * @brief An auxiliary class. The class StringTokenzier is a java style utility class for retrieving tokens separated by delims.
 */
class StringTokenizer{
public:
	/**
	 * @brief Construct a StringTokenizer for specific string
	 * @param s string
	 * @param delims characters considered as separators.
	 */
	StringTokenizer( string s, string delims = "\t\n ");

	/**
	 * @brief Get the next avaiable token.
	 * @return a token
	 */
	string nextToken();

	/**
	 * @brief Check if there're more tokens left.
	 * @return true if have more tokens; otherwise false
	 */
	bool hasMoreTokens() const;

	/**
	 * @brief Get all tokens in the string.
	 * @param tokens stores the tokens
	 */
	void getAllTokens( vector<string>& tokens);
private:
	string s;
	vector<char> delims;
	int index;
};

class String {
public:
	static string toUpper( const string&);
	static string toLower( const string&);
};


class Debug {
public:
	static void check( const bool assertion, string comment);
};

/**
 * @brief An auxiliary class. The class is used for generating random values
 */
class Random
{
public:
	/**
	 * @brief Get a random integer from 0 to r - 1
	 * @param r upper bound of integer value (exclusive)
	 * @return an integer
	 */
	static int nextInt( const int r);

	/**
	 * @brief Get a random integer from r1 to r2
	 * @param r1 tight lower bound
	 * @param r2 tight upper bound
	 * @return an integer
	 */
	static int nextInt( const int r1, const int r2);

	/**
	 * @brief Get a random double from 0 to range.
	 * @param range upper bound
	 * @return a double value
	 */
	static double nextDouble( const double range);

	/**
	 * @brief Get a random double from r1 to r2
	 * @param r1 lower bound
	 * @param r2 upper bound
	 * @return a double value
	 */
	static double nextDouble( const double r1, const double r2);

	/**
	 * @brief Draw a random value from normal distribution N( 0, 1)
	 * @return a double value
	 */
	static double nextNormal();

	/**
	 * @brief Draw a random value from normal distribution N( mean, std * std)
	 * @param mean mean of normal distribution
	 * @param std standard deviation of normal distribution
	 * @return a double value
	 */
	static double nextNormal( const double mean, const double std);
};


/**
 * @brief An auxiliary class performing binary search in a given sorted vector
 */
class BinarySearch {
public:
	/**
	 * @brief Find nearest value in dataset less than or equal to v.
	 */
	static int search( const vector<double>& dataset, const double& v);
private:

	/**
	 * @brief Find nearest value in dataset less than or equal to v.
	 */
	static int search( const vector<double>& dataset, const double& v, const int s, const int e);
};

#endif /* UTILITY_H_ */
