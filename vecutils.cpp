#include <math.h>
#include "vecutils.h"

void zero(VEC3 v)
{
v[0] = v[1] = v[2] = 0;
}

void crossProduct(const VEC3 v1, const VEC3 v2, VEC3 result)
{
    result[0]=(v1[1] * v2[2] - v1[2] * v2[1]);
    result[1]=(v1[2] * v2[0] - v1[0] * v2[2]);
    result[2]=(v1[0] * v2[2] - v1[2] * v2[0]);
}

CALreal magnitude(const VEC3 v){
return sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]));
}


void diff(const VEC3 v1,const VEC3 v2,VEC3 result){
result[0] = v1[0] - v2[0];
result[1] = v1[1] - v2[1];
result[2] = v1[2] - v2[2];
}

CALreal dot(const VEC3 v1,const VEC3 v2){
CALreal result;
for(int i=0;i<3;i++){
result+= v1[i]*v2[i];
}
return result;
}

void normalise(VEC3 v)
{
	// Calculate the magnitude of our vector
	CALreal mag = magnitude(v);

	// As long as the magnitude isn't zero, divide each element by the magnitude
	// to get the normalised value between -1 and +1
	if (mag != 0)
	{
	v[0] /= mag;
	v[1] /= mag;
	v[2] /= mag;
	}
}


CALreal getDistance(const VEC3 v1, const VEC3 v2)
{
	CALreal dx = v2[0] - v1[0];
	CALreal dy = v2[1] - v1[1];
	CALreal dz = v2[2] - v1[2];

	return sqrt(dx * dx + dy * dy + dz * dz);
}	

 void scale(const VEC3 v, const CALreal value, VEC3 res)
{
            res[0]= v[0]*value;
            res[1]= v[1] * value;
            res[2]= v[2] * value;
}

 void add(const VEC3 v, const  VEC3 v2, VEC3 res)
{
            res[0]= v[0]*v2[0];
            res[1]= v[1]*v2[1];
            res[2]= v[2]*v2[2];
}

//void scale(const VEC3& v, const CALreal value, VEC& res)
//{
//            res[0]= v[0]*value;
//            res[1]= v[1] * value;
//            res[2]= v[2] * value;
//}
//        


