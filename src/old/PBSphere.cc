#include "PBSphere.h"

PBSphere::PBSphere(){
    position = Vector3(0,0,0);
    this->radius = 0;
}

PBSphere::PBSphere(Vector3 const &position, double radius){
    this->position = Vector3(position);
    this->radius = radius;
}

PBSphere::PBSphere(double x, double y, double z, double radius){
    position = Vector3(x, y, z);
    this->radius = radius;
}

PBSphere::~PBSphere(){
}


