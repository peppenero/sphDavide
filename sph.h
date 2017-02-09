/*
 * sph.h
 *
 *  Created on: Feb 9, 2017
 *      Author: peppe
 */

#ifndef SRC_SPH_H_
#define SRC_SPH_H_

#include<OpenCAL/calCommon.h>
#include "vecutils.h"

CALreal WPoly6(const CALreal r2, const CALreal h);

bool isNeigh(VEC3 v1, VEC3 v2);

void calcolaDensita(struct CALModel3D* ca, int i, int j, int k);

void computePressureAcceleration(struct CALModel3D* ca, int i, int j, int k);

#endif /* SRC_SPH_H_ */
