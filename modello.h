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
#define X 0.5
#define Y 1
#define Z 0.5


//Cell side for divisions
#define CLD 1000

#define MAX_NUMBER_OF_PARTICLES_PER_CELL 4
#define NODATA -10 // No particle condition (used in px, py and pz)
#define PARTICLE_EDGE   -3
#define PARTICLE_ABSENT  0
#define PARTICLE_PRESENT 1


// Cell side  in m
constexpr const double CL = 0.009;
constexpr const double wx_min = 0; //world xlestmost coordinate
constexpr const double wx_max = 0.5d; //world xlestmost coordinate

constexpr const double wy_min = 0; //world xlestmost coordinate
constexpr const double wy_max = 1.0d; //world xlestmost coordinate

constexpr const double wz_min = 0; //world xlestmost coordinate
constexpr const double wz_max = 0.5d; //world xlestmost coordinate

// Domain dimensions in rows, columns and layers
constexpr const int ROWS = (wx_max-wx_min)/CL;
constexpr const int COLS = (wy_max-wy_min)/CL;
constexpr const int SLICES = (wz_max-wx_min)/CL;



//Sottostati

struct Substates
{
    struct CALSubstate3Dr *px[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *py[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *pz[MAX_NUMBER_OF_PARTICLES_PER_CELL];

    struct CALSubstate3Dr *vx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *vy[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *vz[MAX_NUMBER_OF_PARTICLES_PER_CELL];

    struct CALSubstate3Dr *ax[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *ay[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *az[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    
    struct CALSubstate3Dr *nx[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *ny[MAX_NUMBER_OF_PARTICLES_PER_CELL];
    struct CALSubstate3Dr *nz[MAX_NUMBER_OF_PARTICLES_PER_CELL];

    struct CALSubstate3Dr *mass[MAX_NUMBER_OF_PARTICLES_PER_CELL];
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
#define STEPS 100
#define VERBOSE

// Functions
void pezziala(int , struct CALModel3D* , int , int , int );
void moviliCazzu(struct CALModel3D* , int , int , int );
CALbyte ncestiFluiduNtraStuSlot(struct CALModel3D*, int, int, int,int);
CALbyte ncestiArmenuNaParticellaNtraIVicini(struct CALModel3D*, int, int, int, int);
void transition(struct CALModel3D*,int,int,int);
void startModello();
void initFunction();
void finalizeModel();
void danciNaPosizioni(struct CALModel3D*, const CALreal, const CALreal, const CALreal,const CALreal,const CALreal,const CALreal,const CALreal,const CALint);
void partilu();


#endif /* PROVA_H_ */
