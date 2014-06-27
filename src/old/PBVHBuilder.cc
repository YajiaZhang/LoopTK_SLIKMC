#include "PBVHBuilder.h"

void PBVHBuilder::buildBVH(PBNode *rootNode, Vector3 *points, int number){
    vector<Vector3*> *pts = new vector<Vector3*>();
    for(int i=0; i<number; i++) pts->push_back(&(points[i]));
    PBVHBuilder::buildBVH(NULL, rootNode, pts);
    pts->clear();
}

void PBVHBuilder::buildSingleNode(PBNode *rootNode, Vector3 *points, int number){
    vector<Vector3*> *pts = new vector<Vector3*>();
    for(int i=0; i<number; i++) pts->push_back(&(points[i]));
    PBVHBuilder::buildSingleNode(rootNode, pts);
    pts->clear();
    delete pts;
}

void PBVHBuilder::buildBVH(PBNode* rootNode, vector<Vector3*>*points){
    PBVHBuilder::buildBVH(NULL, rootNode, points);
}

void PBVHBuilder::buildBVH(PBNode *parentNode, PBNode *currentNode, vector<Vector3*> *points){
    //1. build nodes by given points
    PBVHBuilder::buildNode(currentNode, points);
    currentNode->storePoints(points);
    //2. set parent node
    currentNode->setParent(parentNode);
    //2.1 resursive end -> if the nodes contain mininum number of points
    if (points->size() <= PBVHBuilder::MIN_POINT_NUM) return;

    //3. get principal conponent
    gsl_vector *pc = gsl_vector_alloc(3);
    PPCA::getPrincipalComponent(pc, points);

    //4. get mass points
    Vector3 massP = PBVHBuilder::getMassPoint(points);
    
    //5. separate points into two set
    vector<Vector3*> *underPoints = new vector<Vector3*>;
    vector<Vector3*> *abovePoints = new vector<Vector3*>;
    PBVHBuilder::separatePoints(pc, massP, points, underPoints, abovePoints);

    //6. recursively build bounding spheres 
    if (underPoints->size()>0){
        PBNode* childNode = new PBNode();
        PBVHBuilder::buildBVH(currentNode, childNode, underPoints);
        currentNode->addChild(childNode); 
    }
    if (abovePoints->size()>0){
        PBNode* childNode = new PBNode();
        PBVHBuilder::buildBVH(currentNode, childNode, abovePoints);
        currentNode->addChild(childNode);
    }

    //F: release resources
    underPoints->clear();
    abovePoints->clear();
    delete pc;
}

void PBVHBuilder::buildSingleNode(PBNode *targetNode, vector<Vector3*> *points){
    //1. build nodes by given points
    PBVHBuilder::buildNode(targetNode, points);
    targetNode->storePoints(points);
}

void PBVHBuilder::separatePoints(gsl_vector *normal, Vector3 refPoint, const vector<Vector3*> *points,
            vector<Vector3*> *underPts, vector<Vector3*> *abovePts){
    double c = PBVHBuilder::getPlaneDistance(normal, refPoint);    
    vector<Vector3*> *equalPts = new vector<Vector3*>();
    double distance;
    //printf("refP:%f, %f, %f\n", refPoint.x, refPoint.y, refPoint.z);
    for(int i=0; i<points->size(); i++){

        distance = PBVHBuilder::getPlaneDistance(normal, *(*points)[i]);
      //  printf("P:%f, %f, %f\n", (*(*points)[i]).x,(*(*points)[i]).y,(*(*points)[i]).z );
  
        // cout << "c:"<<c<<" Dis:"<< distance<<endl;
   
        if (distance>c)
            abovePts->push_back((*points)[i]);
        else if (distance<c)
            underPts->push_back((*points)[i]);
        else 
            equalPts->push_back((*points)[i]);
    }
 
     //printf("Before:asize:%d, bsize:%d\n", abovePts->size(), underPts->size());

    if (equalPts->size()>0){
        if (abovePts->size()==0 && underPts->size()==0){
            for(int i=0; i<ceil(equalPts->size()/2.0); i++) abovePts->push_back((*equalPts)[i]);
            for(int i=(int)ceil(equalPts->size()/2.0); i<equalPts->size(); i++) underPts->push_back((*equalPts)[i]);
        }else if (abovePts->size()==0)
            for(int i=0; i<equalPts->size(); i++) abovePts->push_back((*equalPts)[i]);
        else if (underPts->size()==0)
            for(int i=0; i<equalPts->size(); i++) underPts->push_back((*equalPts)[i]);
        else if (abovePts->size() > underPts->size())
            for(int i=0; i<equalPts->size(); i++) abovePts->push_back((*equalPts)[i]);
        else        
            for(int i=0; i<equalPts->size(); i++) underPts->push_back((*equalPts)[i]);
    }
    equalPts->clear();

    //printf("asize:%d, bsize:%d\n", abovePts->size(), underPts->size());

}

void PBVHBuilder::buildNode(PBNode *targetNode, Vector3 *points, int number){
    vector<Vector3*> *pts;
    for(int i=0; i<number; i++){
        pts->push_back(&points[i]);
    }
    PBVHBuilder::buildNode(targetNode, pts);
    pts->clear();
}

void PBVHBuilder::buildNode(PBNode *targetNode, vector<Vector3*> *points){
    double radius =0;
    Vector3 massP = getMassPoint(points);
    double tempRadius = 0;
    for(int i=0; i<points->size(); i++){
        tempRadius = PBVHBuilder::getDistance(massP, *(*points)[i]);
        if (tempRadius > radius) radius = tempRadius;
    }
    PBSphere *sphere = new PBSphere(massP, radius);
    targetNode->setSphere(sphere);

}

Vector3 PBVHBuilder::getMassPoint(Vector3 *points, int number){
    Vector3 pt = Vector3(0, 0, 0);
    for(int i=0; i<number; i++){
        pt.x += points[i].x/(double)number;
        pt.y += points[i].y/(double)number;
        pt.z += points[i].z/(double)number;
    } 
    return pt;  
}

Vector3 PBVHBuilder::getMassPoint(vector<Vector3*> *points){
    Vector3 pt = Vector3(0, 0, 0);
    int number = points->size();
    for(int i=0; i<number; i++){
        pt.x += (*(*points)[i]).x/(double)number;
        pt.y += (*(*points)[i]).y/(double)number;
        pt.z += (*(*points)[i]).z/(double)number;
    }
    return pt;
}


double PBVHBuilder::getDistance(Vector3 p1, Vector3 p2){
    return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z));

}

void PBVHBuilder::pointsToArray(Vector3 *points, double *values, int number){
    for(int i =0; i<number; i++){
        values[i*3] = points[i].x;
        values[i*3+1] = points[i].y;
        values[i*3+2] = points[i].z;
    }
}

Vector3 PBVHBuilder::getLargestEigenvector(Vector3 *points, int number){
    double values[number*3];
    PBVHBuilder::pointsToArray(points, values, number);
    gsl_matrix_view m  = gsl_matrix_view_array (values, number, 3);
    gsl_vector *eval = gsl_vector_alloc (3);
    gsl_matrix *evec = gsl_matrix_alloc (3, 3);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);                     

    gsl_vector_view evec_i  = gsl_matrix_column (evec, 0);
    gsl_vector_get(&(evec_i.vector),0);
    Vector3 eigenvector = Vector3(gsl_vector_get(&(evec_i.vector),0), 
            gsl_vector_get(&(evec_i.vector),1),gsl_vector_get(&(evec_i.vector),2));
    return eigenvector;    
}
