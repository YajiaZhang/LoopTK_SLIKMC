#ifndef PBVHBUILDER_H
#define PBVHBUILDER_H
#include "PLibraries.h"
#include "PBNode.h"    
#include "PPCA.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

class PBVHBuilder{

    public:
    const static int MIN_POINT_NUM = 1;
    //given a set of points to build a hierarchical boulding sphere
    static void buildBVH(PBNode *rootNode, Vector3 *points, int number);
    static void buildBVH(PBNode *rootNode, vector<Vector3*> *points);
    static void buildBVH(PBNode *parentNode, PBNode *currentNode, vector<Vector3*> *points);
    //given a set of points to build a bounding sphere
    static void buildNode(PBNode *targetNode, Vector3 *points, int number);
    static void buildNode(PBNode *targetNode, vector<Vector3*> *points);
    static void buildSingleNode(PBNode *targetNode, vector<Vector3*> *points);
    static void buildSingleNode(PBNode *targetNode, Vector3 *points, int number);
    static Vector3 getMassPoint(Vector3 *points, int number); 
    static Vector3 getMassPoint(vector<Vector3*> *points);
    static double getDistance(Vector3 p1, Vector3 p2);
    static void pointsToArray(Vector3 *points, double *values, int number);
    static Vector3 getLargestEigenvector(Vector3 *points, int number);
    static void separatePoints(gsl_vector *normal, Vector3 refPoint, const vector<Vector3*> *points, 
            vector<Vector3*> *underPts, vector<Vector3*> *abovePts);  
    
    private: 
    inline static double getPlaneDistance(gsl_vector *normal, Vector3 p){
        return gsl_vector_get(normal, 0)*p.x + gsl_vector_get(normal,1)*p.y+
                    gsl_vector_get(normal, 2)*p.z;
    }    
};

#endif
