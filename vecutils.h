#ifndef VECUTILS_H_
#define VECUTILS_H_

#include <OpenCAL/cal3D.h>

typedef CALreal VEC3[3];

void crossProduct(const VEC3, const VEC3, VEC3);

void zero(VEC3);

CALreal magnitude(const VEC3);

void normalise(VEC3);

void diff(const VEC3,const VEC3,VEC3);
	
CALreal dot(const VEC3,const VEC3);

CALreal getDistance(const VEC3, const VEC3);

void scale(const VEC3 v, const CALreal value, VEC3 res);
void mult(VEC3 , VEC3 ,const CALreal );

void add( const VEC3 , const  VEC3 , VEC3);

#endif
