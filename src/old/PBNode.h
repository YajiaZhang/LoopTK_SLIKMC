#ifndef PBNODE_H
#define PBNODE_H
#include "PBSphere.h"

class PBNode{
  public:
    static const int PRINT_ALL_NODE = 1;
    //from Johannes Kepler
    static const double PACKING_DENSITY = 0.74048;
    static const double NEAR_SURFACE_THRESHOLD = 0.90;
    static const double NEAR_THRESHOLD = 0.35;

    //static const double NEAR_THRESHOLD = 1.5;

    PBNode(){
        points = new vector<Vector3*>;
        childNum = 0;
    }

    PBNode(PBSphere *sphere){
        childNum = 0;
        this->sphere = sphere;
    }

    ~PBNode(){
        delete sphere;
        delete child;
        //delete parent;
        points->clear();
        delete points;
        childNum = 0;
    }

    void reset(){
        this->sphere = NULL;
        this->child[0] = NULL;
        this->child[1] = NULL;
        this->parent = NULL;
        this->points->clear();
        delete points;        
        delete sphere;
        this->childNum = 0;
    }    
    
    void storePoints(vector<Vector3*> *pts){
       points = new vector<Vector3*>;
        for(int i=0; i<pts->size(); i++){
            points->push_back((*pts)[i]);
        }
    }

    vector<Vector3*> getPointsNearSurface(){
        vector<Vector3*> vec;
        Vector3 center = sphere->getPosition();
        for(int i=0; i<points->size(); i++){
            if ( getDistance((*points)[i], &center )>sphere->getRadius()*NEAR_SURFACE_THRESHOLD){
                vec.push_back( (*points)[i]);
            }
        }
        return vec;
    }

    vector<Vector3*> getPointsNear(Vector3 pt, double threshold){
        vector<Vector3*> vec;
        for(int i=0; i<points->size(); i++){
            if ( getDistance((*points)[i], &pt ) <= threshold)
                vec.push_back( (*points)[i]);
        }
        return vec;
    }

    vector<Vector3*> getPointsNear(Vector3 pt){
        vector<Vector3*> vec;
        double threshold = sphere->getRadius()*NEAR_THRESHOLD;
        for(int i=0; i<points->size(); i++){
            if ( getDistance((*points)[i], &pt ) <= threshold)
                vec.push_back( (*points)[i]);
        }
        return vec;
    }

    double getMaxPackingNumber(double maxRadius){
       if (this->sphere!=NULL){
         double sphereRadius = this->sphere->getRadius()+maxRadius;
         double maxNum = sphereRadius*sphereRadius*sphereRadius*PACKING_DENSITY/(maxRadius*maxRadius*maxRadius);
         return maxNum;
       }else{
         return -1;
       }
    }     

    double getMaxVolumeNumber(double maxRadius){
       if (this->sphere!=NULL){
         double sphereRadius = this->sphere->getRadius()+maxRadius;
         double maxNum = sphereRadius*sphereRadius*sphereRadius/(maxRadius*maxRadius*maxRadius);
         return maxNum;
       }else{
         return -1;
       }
	}

    double getPackingRatio(double maxRadius){
       return this->getPointNumber()/this->getMaxPackingNumber(maxRadius);
    }

    double getVolumeRatio(double maxRadius){
       if (this->sphere!=NULL){
           double sphereRadius = this->sphere->getRadius()+maxRadius;
           double maxNum = sphereRadius*sphereRadius*sphereRadius/(maxRadius*maxRadius*maxRadius);
           return this->getPointNumber()/maxNum;
       }else{
           return -1;
       }
    }
    
    vector<Vector3 *> *getPoints(){
        return points;
    }

    int getPointNumber(){
        return points->size();
    }

    PBSphere* getSphere(){
        return this->sphere;
    }

    void setSphere(PBSphere *sphere){
        this->sphere = sphere;
    }

    void addChild(PBNode *node){
        if (childNum>=2) return;
        child[childNum]= node;
        childNum++;
    }

    void removeChild(int index){
        if (index<0 || index>childNum-1) return;
        child[index] = NULL;
        childNum--;
    }

    PBNode *getChild(int index){
        if (index<childNum && index >=0){
            return child[index];
        }
        return NULL;
    }

    void setParent(PBNode *parent){
        this->parent = parent;
    }

    PBNode * getParent(){
        return this->parent;
    }

    void removeParent(){
        this->parent = NULL;
    }

    void toString(char *str){
        Vector3 pos = sphere->getPosition();
        sprintf(str, "PBNode(%f,%f,%f,%f) childNum:%d, isRoot:%d, pointsNum:%d", pos.x, pos.y, pos.z, 
                sphere->getRadius(), childNum, (parent==NULL), getPointNumber());
    }

    void printHierarchy(){
        int level = 1;
        PBNode* currentNode = this;
        cout << "level:"<<level<<endl;
        if (currentNode== NULL) return;
        else printNode(currentNode);
        PBNode* child1 = currentNode->getChild(0); 
        PBNode* child2 = currentNode->getChild(1);
        if (child1!=NULL)
            printHierarchy(child1, level+1);
        if (child2!=NULL)
            printHierarchy(child2, level+1);
    }

    void printHierarchy(PBNode* currentNode, int l){
        cout << "level:"<<l<<endl;
        if (currentNode== NULL) return;
        else printNode(currentNode);
        PBNode* child1 = currentNode->getChild(0);
        PBNode* child2 = currentNode->getChild(1);
        if (child1!=NULL)
           printHierarchy(child1, l+1);
        if (child2!=NULL)
           printHierarchy(child2, l+1);
    }
    
    void printPackingNumber(double maxRadius){
        int level = 1;
        PBNode* currentNode = this;
        cout <<"level:"<<level ;
        if (currentNode == NULL) return;
        else cout<<" MaxPacking:"<< currentNode->getMaxPackingNumber(maxRadius) << 
            " RealPacking:"<< currentNode->getPointNumber() << 
            " Ratio:"<< currentNode->getPackingRatio(maxRadius)<<endl;
        PBNode* child1 = currentNode->getChild(0);
        PBNode* child2 = currentNode->getChild(1);
        if (child1 != NULL)
            printPackingNumber(child1, level+1, maxRadius);
        if (child2 != NULL)
            printPackingNumber(child2, level+1, maxRadius);
    }

    void printPackingNumber(PBNode* currentNode, int l, double maxRadius){
        cout <<"level:"<<l ;
        if (currentNode == NULL) return;
        else cout<<" MaxPacking:"<< currentNode->getMaxPackingNumber(maxRadius) <<
            " RealPacking:"<< currentNode->getPointNumber() <<
            " Ratio:"<< currentNode->getPackingRatio(maxRadius)<<endl;
        PBNode* child1 = currentNode->getChild(0);
        PBNode* child2 = currentNode->getChild(1);
        if (child1 != NULL)
            printPackingNumber(child1, l+1, maxRadius);
        if (child2 != NULL)
            printPackingNumber(child2, l+1, maxRadius);                                                                         
    }

    void printVolumeRatio(double maxRadius){
        int level = 1;
        PBNode* currentNode = this;
        cout <<"level:"<<level ;
        if (currentNode == NULL) return;
        else cout<<" MaxVolumes:"<< currentNode->getMaxVolumeNumber(maxRadius) <<
            " RealPacking:"<< currentNode->getPointNumber() <<
            " Ratio:"<< currentNode->getVolumeRatio(maxRadius)<<endl;
        PBNode* child1 = currentNode->getChild(0);
        PBNode* child2 = currentNode->getChild(1);
        if (child1 != NULL)
            printVolumeRatio(child1, level+1, maxRadius);
        if (child2 != NULL)
            printVolumeRatio(child2, level+1, maxRadius);
    }

    void printVolumeRatio(PBNode* currentNode, int l, double maxRadius){
        cout <<"level:"<<l ;
        if (currentNode == NULL) return;
        else cout<<" MaxVolumes:"<< currentNode->getMaxVolumeNumber(maxRadius) <<
            " RealPacking:"<< currentNode->getPointNumber() <<
            " Ratio:"<< currentNode->getVolumeRatio(maxRadius)<<endl;
        PBNode* child1 = currentNode->getChild(0);
        PBNode* child2 = currentNode->getChild(1);
        if (child1 != NULL)
            printVolumeRatio(child1, l+1, maxRadius);
        if (child2 != NULL)
            printVolumeRatio(child2, l+1, maxRadius);                                                           
    }
   
    void printNode(PBNode* node){
        char msg[100];
        node->toString(msg);
        cout << msg<<endl;
        if (PRINT_ALL_NODE)
            node->printPoints();
    }

    void printPoints(){
        for(int i=0; i<points->size(); i++){
            printf("(%f, %f, %f)\n", (*(*points)[i]).x, (*(*points)[i]).y, (*(*points)[i]).z);
        }
    }
   
    double getDistance(Vector3 *p1, Vector3 *p2){
         return sqrt((p1->x-p2->x)*(p1->x-p2->x)+
                 (p1->y-p2->y)*(p1->y-p2->y)+
                 (p1->z-p2->z)*((p1->z-p2->z)));
    }

    double getPointVolumeRatio(){
        double volume = 3.0/4.0 * 3.14159 * 
            (this->sphere->getRadius())*(this->sphere->getRadius())*(this->sphere->getRadius());
        return this->getPointNumber()/volume;
    }

  private:
    PBSphere *sphere;
    PBNode *child[2];
    PBNode *parent;
    vector<Vector3*> *points;
    int childNum;
};

#endif
