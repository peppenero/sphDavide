/*
 * sph.h
 *
 *  Created on: Feb 9, 2017
 *      Author: peppe
 */

#ifndef SRC_SPH_H_
#define SRC_SPH_H_
#include<glm/glm.hpp>
extern"C"{
#include<OpenCAL/calCommon.h>

}

typedef glm::tvec3<int> VEC3i;
typedef glm::tvec3<double> VEC3r;


CALreal WPoly6(const CALreal r2, const CALreal h);

bool isNeigh(VEC3r v1, VEC3r v2);

void calcolaDensita(struct CALModel3D* ca, int i, int j, int k);

void computePressureAcceleration(struct CALModel3D* ca, int i, int j, int k);

VEC3r computeExternalForces(struct CALModel3D* ca, int i, int j, int k,int slot);

void advance(struct CALModel3D* ,const int , const int , const int, const int, const VEC3r& n);


#endif /* SRC_SPH_H_ */
