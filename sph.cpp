#include "modello.h"
#include "sph.h"
#include <math.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr const static  double STIFFNESS  = 3.0; // Nm/kg is gas constant of water vapor
constexpr const static  double REST_DENSITY = 998.29; //kg/m^3 is rest density of water particle
constexpr const static  double VISCOSITY = 3.5; // Ns/m^2 or Pa*s viscosity of water

constexpr const static  double  SURFACE_TENSION = 0.0728; // N/m
constexpr const static  double  SURFACE_THRESHOLD = 7.065;

constexpr const static  double GRAVITY_ACCELERATION = 9.80665;


constexpr const static  double WALL_DAMPING =-0.9; // wall damping constant
constexpr const static  double WALL_K = 10000.0; // wall spring constant



//aggiusta i valori di questi sotto
constexpr const static  double RADIUS =0.000045; // particle radius
constexpr const static  double MASS =0.000026254; // particle mass
constexpr const static  double DT =0.000009; // time simulation quantum




CALreal WPoly6(const CALreal r2, const CALreal h) {
    const CALreal coefficient = 315.0/(64.0*M_PI*pow(h,9));
    const CALreal hSquared = pow(h,2);

    return coefficient * pow(hSquared-r2, 3);
}

void Wpoly6Gradient(VEC3r& diffPosition, CALreal radiusSquared, VEC3r& gradient , const CALreal h) {
    static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
    static double hSquared = pow(h,2);

    gradient =  coefficient * pow(hSquared-radiusSquared, 2) * diffPosition ;

}

inline double Wpoly6Laplacian(const double radiusSquared, const double h) {

    static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
    static double hSquared = pow(h,2);

    return coefficient * (hSquared-radiusSquared) * (3.0*hSquared - 7.0*radiusSquared);
}


void WspikyGradient(VEC3r& diffPosition, CALreal radiusSquared, VEC3r& gradient , const CALreal h) {
    static CALreal coefficient = -45.0/(M_PI*pow(h,6));

    CALreal radius = sqrt(radiusSquared);

    gradient = coefficient * pow(h-radius, 2) * diffPosition / radius;
}


inline CALreal WviscosityLaplacian(CALreal radiusSquared, const CALreal h) {

    static CALreal coefficient = 45.0/(M_PI*pow(h,6));

    CALreal radius = sqrt(radiusSquared);

    return coefficient * (h - radius);
}



bool isNeigh(VEC3r p1, VEC3r p2) {
    double d = glm::distance(p1,p2);
    if(d<=RADIUS)
        return true;
    return false;

    //return magnitude(d)<DISTANCE;
}

void calcolaDensita(struct CALModel3D* ca, int i, int j, int k) {
    for(int slot; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++) {
        CALreal density = 0;


        //pos particle
        const CALreal px = calGet3Dr(ca,Q.px[slot],i,j,k);
        const CALreal py = calGet3Dr(ca,Q.py[slot],i,j,k);
        const CALreal pz = calGet3Dr(ca,Q.pz[slot],i,j,k);
        VEC3r p1 = VEC3r(px,py,pz);

        for (int n=0; n<ca->sizeof_X; n++) {
            for(int slot1; slot1<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot1++) {

                const CALreal npx = calGetX3Dr(ca,Q.px[slot1],i,j,k,n);
                const CALreal npy = calGetX3Dr(ca,Q.py[slot1],i,j,k,n);
                const CALreal npz = calGetX3Dr(ca,Q.pz[slot1],i,j,k,n);
                VEC3r p2= VEC3r(npx,npy,npz);

                //v1 and v2 does not refer to the very same particle
                if(!(n==0 && slot1==slot) && isNeigh(p1,p2)) {
                    VEC3r d = p1-p2;

                    CALreal dist2 = glm::dot(d,d);
                    density += WPoly6(dist2,RADIUS);
                }
            }
        }
        density = density*MASS;

        const CALreal pressure = STIFFNESS* ( density - REST_DENSITY );

        calSet3Dr(ca,Q.density[slot],i,j,k,density);
        calSet3Dr(ca,Q.pressure[slot],i,j,k,pressure);
    }

}


inline VEC3r getPositionX(struct CALModel3D* ca , const int i, const int j, const int k, const int slot, const int n) {
    return VEC3r(calGetX3Dr(ca,Q.px[slot],i,j,k,n),calGetX3Dr(ca,Q.py[slot],i,j,k,n),calGetX3Dr(ca,Q.pz[slot],i,j,k,n));
}

inline  VEC3r getVelocityX(struct CALModel3D * ca , const int i, const int j, const int k, const int slot, const int n) {
    return VEC3r(calGetX3Dr(ca,Q.vx[slot],i,j,k,n),calGetX3Dr(ca,Q.vy[slot],i,j,k,n),calGetX3Dr(ca,Q.vz[slot],i,j,k,n));
}


inline VEC3r getPosition(struct CALModel3D* ca , const int i, const int j, const int k, const int slot) {
    return VEC3r(calGet3Dr(ca,Q.px[slot],i,j,k),calGet3Dr(ca,Q.py[slot],i,j,k),calGet3Dr(ca,Q.pz[slot],i,j,k));
}

inline  VEC3r getVelocity(struct CALModel3D * ca , const int i, const int j, const int k, const int slot) {
    return VEC3r(calGet3Dr(ca,Q.vx[slot],i,j,k),calGet3Dr(ca,Q.vy[slot],i,j,k),calGet3Dr(ca,Q.vz[slot],i,j,k));
}

inline  VEC3r getNormal(struct CALModel3D * ca , const int i, const int j, const int k, const int slot) {
    return VEC3r(calGet3Dr(ca,Q.nx[slot],i,j,k),calGet3Dr(ca,Q.ny[slot],i,j,k),calGet3Dr(ca,Q.nz[slot],i,j,k));
}

VEC3r G = VEC3r(0,0,9.81);

void computePressureAcceleration(struct CALModel3D* ca, int i, int j, int k) {

    for(int slot=0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++) {
        //---------------Messo io l'if------------------
        if(ncestiFluiduNtraStuSlot(ca,i,j,k,slot)) {//solo se sono di fluido

            CALreal pdensity = calGet3Dr(ca,Q.density[slot],i,j,k);
            CALreal ppressure = calGet3Dr(ca,Q.pressure[slot],i,j,k);

            VEC3r f_gravity = G * pdensity;
            VEC3r f_pressure=VEC3r(0,0,0);
            VEC3r f_viscosity=VEC3r(0,0,0);
            VEC3r f_surface=VEC3r(0,0,0);
            VEC3r a_external=VEC3r(0,0,0);
            VEC3r colorFieldNormal=VEC3r(0,0,0);
            double colorFieldLaplacian;


            //pos particle
            VEC3r p1 = getPosition(ca,i,j,k,slot);



            //vel particle

            VEC3r v1 = getVelocity(ca,i,j,k,slot);


            //norm particle
            VEC3r n1 = getNormal(ca,i,j,k,slot);

            int cn=0;
            for (int n=0; n<ca->sizeof_X; n++) {
                for(int slot1; slot1<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot1++) {

                    VEC3r p2 = getPositionX(ca,i,j,k,slot1,n);
                    VEC3r v2 = getVelocityX(ca,i,j,k,slot1,n);
                    CALreal ndensity = calGetX3Dr(ca,Q.density[slot1],i,j,k,n);
                    CALreal npressure = calGetX3Dr(ca,Q.pressure[slot1],i,j,k,n);

                    if(isNeigh(p1,p2)) {
                        cn++;
                        VEC3r diff = p1 - p2;
                        double dist_2 = glm::dot(diff,diff);

                        VEC3r poly6Gradient;
                        Wpoly6Gradient(diff, dist_2, poly6Gradient,RADIUS);

                        VEC3r spikyGradient;
                        WspikyGradient(diff, dist_2, spikyGradient, RADIUS);

                        if(!(n==0 && slot1==slot)) {
                            f_pressure +=(  ppressure/pow(pdensity,2)+
                                            npressure/pow(ndensity,2)
                                         )*spikyGradient;

                            f_viscosity += ( v2 -v1)*
                                           WviscosityLaplacian(dist_2,RADIUS) / ndensity;
                        }

                        colorFieldNormal += poly6Gradient / ndensity;
                        colorFieldLaplacian += Wpoly6Laplacian(dist_2,RADIUS) / ndensity;

                    }

                }
            }

            f_pressure *= -MASS* pdensity;
            f_viscosity *=  VISCOSITY * MASS;

            colorFieldNormal *= MASS;
            //!!!!!!p->normal = -1.0 * colorFieldNormal; AGGIUSTA QUESTO---------------------------------------------------------------
            n1 = -1.0 * colorFieldNormal;
            colorFieldLaplacian *=MASS;

            // surface tension force
            double colorFieldNormalMagnitude = colorFieldNormal.length(); //check that this is actually the magnitude of the vector

            if (colorFieldNormalMagnitude >  SURFACE_THRESHOLD) {

                calSet3Db(ca,Q.flag[slot],i,j,k,1);
                f_surface = - SURFACE_TENSION * colorFieldNormal / colorFieldNormalMagnitude * colorFieldLaplacian;

            }
            else {
                calSet3Db(ca,Q.flag[slot],i,j,k,0);
            }

            VEC3r acc = (f_pressure  +f_viscosity + f_surface + f_gravity ) / pdensity;

            printf("%f\n",f_pressure[0]);
            printf("%f\n",f_viscosity[0]);
            printf("%f\n",f_surface[0]);
            printf("%f\n",f_gravity[0]);
            printf("%f\n",pdensity);



            /**add external forces, collision, user interaction, moving rigid bodys etc.**/
            a_external = computeExternalForces(ca,i,j,k,slot);
            acc+= a_external;
            //acc = VEC3r(0.0,0.0,-0.1);
            advance(ca,i,j,k,slot,acc);
        }
    }

}


template<class T>
T ipow(T base, unsigned int exp) {
    T result = 1;
    while (exp) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

void advance(struct CALModel3D* ca , const int i, const int j, const int k , const int slot, const VEC3r& acc) {

    VEC3r p1,v1;
    p1 = getPosition(ca,i,j,k,slot);
    v1 = getVelocity(ca,i,j,k,slot);
    VEC3r newPosition = p1 + v1*DT + acc * ipow<double>(DT,2);
    VEC3r newVelocity = (newPosition - p1) /DT;
    printf("Vecchia: %f,%f,%f\n",p1[0],p1[1],p1[2]);
    printf("Nuova: %f,%f,%f\n",newPosition[0],newPosition[1],newPosition[2]);

//set  new position
    calSet3Dr(ca,Q.px[slot],i,j,k,newPosition[0]);
    calSet3Dr(ca,Q.py[slot],i,j,k,newPosition[1]);
    calSet3Dr(ca,Q.pz[slot],i,j,k,newPosition[2]);
//set  new velocity
    calSet3Dr(ca,Q.vx[slot],i,j,k,newVelocity[0]);
    calSet3Dr(ca,Q.vy[slot],i,j,k,newVelocity[1]);
    calSet3Dr(ca,Q.vz[slot],i,j,k,newVelocity[2]);


}

//TODO!--------------------------------------------
bool hitWall(const VEC3r& p1,const VEC3r& p2) {
    if(glm::distance(p1,p2)<=0.02)
        return true;
    return false;
}


VEC3r computeExternalForces(struct CALModel3D* ca, int i, int j, int k,int slot) {

    VEC3r _a_extern(0,0,0);
    VEC3r p1 = getPosition(ca,i,j,k,slot);

    VEC3r v1 = getVelocity(ca,i,j,k,slot);

    for(int slot1=0; slot1<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot1++) {
        VEC3r p2 = getPosition(ca,i,j,k,slot1);
        CALint isWall =  calGet3Di(ca,Q.imove[slot1],i,j,k);
        VEC3r n2 =  getNormal(ca,i,j,k,slot);
        if( isWall==-3 && slot!=slot1 && hitWall(p1,p2)) {
            double d = glm::dot((p2 - p1),n2) + RADIUS;
            if(d > 0) {
                _a_extern +=  WALL_K * n2 * d;
                _a_extern +=  WALL_DAMPING * glm::dot(v1,n2) * n2;
            }
        }

    }
    return _a_extern;
}




