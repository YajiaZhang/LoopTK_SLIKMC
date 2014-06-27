#ifndef PBSPHERE_H
#define PBSPHERE_H

#include "PLibraries.h"

class PBSphere{

   
    private:
    Vector3 position;
    double radius;
    
    public:
    
    PBSphere();
    
    PBSphere(Vector3 const &position, double radius);

    PBSphere(double x, double y, double z, double radius);
 
    ~PBSphere();

    inline Vector3 getPosition(){
        return position;
    }

    inline double getRadius(){
        return radius;
    }
}; 

#endif
