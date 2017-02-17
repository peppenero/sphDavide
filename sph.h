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

constexpr const static  double STIFFNESS  = 3.0; // Nm/kg is gas constant of water vapor
constexpr const static  double REST_DENSITY = 998.29; //kg/m^3 is rest density of water particle
constexpr const static  double VISCOSITY = 3.5; // Ns/m^2 or Pa*s viscosity of water

constexpr const static  double  SURFACE_TENSION = 0.0728; // N/m
constexpr const static  double  SURFACE_THRESHOLD = 7.065;

constexpr const static  double GRAVITY_ACCELERATION = 9.80665;


constexpr const static  double WALL_DAMPING =-0.9; // wall damping constant
constexpr const static  double WALL_K = 10000.0; // wall spring constant


//aggiusta i valori di questi sotto
constexpr const static  double RADIUS =0.016; // particle radius
constexpr const static  double MASS =0.0008; // particle mass
constexpr const static  double DT =0.001; // time simulation quantum


static bool stampa=false;

typedef glm::tvec3<int> VEC3i;
typedef glm::tvec3<double> VEC3r;


CALreal WPoly6(const CALreal r2, const CALreal h);

bool isNeigh(VEC3r v1, VEC3r v2);

void calcolaDensita(struct CALModel3D* ca, int i, int j, int k);

void computePressureAcceleration(struct CALModel3D* ca, int i, int j, int k);

VEC3r computeExternalForces(struct CALModel3D* ca, int i, int j, int k,int slot);

void advance(struct CALModel3D* ,const int , const int , const int, const int, const VEC3r& n);


#endif /* SRC_SPH_H_ */
