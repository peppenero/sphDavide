#ifndef PROVA_H_
#define PROVA_H_

extern "C" {
#include <OpenCAL/cal3D.h>
#include <OpenCAL/cal3DRun.h>
#include <OpenCAL/cal3DIO.h>
#include <OpenCAL/cal3DUnsafe.h>
#include <OpenCAL/calCommon.h>
}

#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>


#define NUMBER_OF_PARTICLES 14988

// Domain dimensions in m
#define X 0.10
#define Y 0.10
#define Z 0.05

// Cell side  in m
#define CL 0.001

//cell side a cazzo
#define CL1 0.01

//Cell side for divisions
#define CLD 1000

#define MAX_NUMBER_OF_PARTICLES_PER_CELL 3
#define NODATA -10 // No particle condition (used in px, py and pz)
#define PARTICLE_EDGE   -3
#define PARTICLE_ABSENT  0
#define PARTICLE_PRESENT 1

// Domain dimensions in rows, columns and layers
#define ROWS    (int)((X)*(CLD))+1
#define COLS    (int)((Y)*(CLD))+1
#define SLICES  (int)((Z)*(CLD))+1





//Sottostati

struct Substates
{
	struct CALSubstate3Dr *px[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *py[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *pz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vy[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *vz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    
    struct CALSubstate3Dr *nx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *ny[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *nz[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    
	struct CALSubstate3Dr *density[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Dr *pressure[MAX_NUMBER_OF_PARTICLES_PER_CELL];
	struct CALSubstate3Di *imove[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    
    struct CALSubstate3Db *flag[MAX_NUMBER_OF_PARTICLES_PER_CELL];
};

// Main objcts
extern struct CALModel3D* u_modellu;
extern struct Substates Q;
extern struct CALRun3D* a_simulazioni;

// Computational steps
#define STEPS 30
#define VERBOSE

// Functions
CALbyte ncestiFluiduNtraStuSlot(struct CALModel3D*, int, int, int,int);
CALbyte ncestiArmenuNaParticellaNtraIVicini(struct CALModel3D*, int, int, int, int);
void transition(struct CALModel3D*,int,int,int);
void startModello();
void initFunction();
void finalizeModel();
void danciNaPosizioni(struct CALModel3D* ca, const CALreal, const CALreal, const CALreal,const CALreal,const CALreal,const CALreal,const CALreal,const CALint imove);
void partilu();


#endif /* PROVA_H_ */
