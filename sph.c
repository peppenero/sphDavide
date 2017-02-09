#include "prova.h"
#include "vecutils.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DISTANCE    (0.001)
#define MASS        (0)
#define STIFFNESS   (0)
#define RADIUS      (0)
#define RESTDENSITY (0)
#define VISCOSITY   (0)
#define SURFACE_THRESHOLD (0)
#define SURFACE_TENSION (0)

CALreal WPoly6(const CALreal r2, const CALreal h) {
    const CALreal coefficient = 315.0/(64.0*M_PI*pow(h,9));
    const CALreal hSquared = pow(h,2);

    return coefficient * pow(hSquared-r2, 3);
}

void Wpoly6Gradient(VEC3 diffPosition, CALreal radiusSquared, VEC3 gradient , const CALreal h) {
    CALreal coefficient = -945.0/(32.0*M_PI*pow(h,9));
    CALreal hSquared = pow(h,2);

    scale(diffPosition, coefficient * pow(hSquared-radiusSquared, 2) , diffPosition);

}

double Wpoly6Laplacian(const double radiusSquared, const double h) {

    double coefficient = -945.0/(32.0*M_PI*pow(h,9));
    double hSquared = pow(h,2);

    return coefficient * (hSquared-radiusSquared) * (3.0*hSquared - 7.0*radiusSquared);
}


void WspikyGradient(VEC3 diffPosition, CALreal radiusSquared, VEC3 gradient , const CALreal h) {

    CALreal coefficient = -45.0/(M_PI*pow(h,6));
    CALreal radius = sqrt(radiusSquared);

    scale(diffPosition,  coefficient * pow(h-radius, 2), diffPosition );
    scale(diffPosition,  1/radius , diffPosition);
}


CALreal WviscosityLaplacian(CALreal radiusSquared, const CALreal h) {

    CALreal coefficient = 45.0/(M_PI*pow(h,6));

    CALreal radius = sqrt(radiusSquared);

    return coefficient * (h - radius);
}




bool isNeigh(VEC3 v1, VEC3 v2) {
    VEC3 d;
    diff(v1,v2,d);
    return true;
    //return magnitude(d)<DISTANCE;
}

void calcolaDensita(struct CALModel3D* ca, int i, int j, int k) {
    for(int slot; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++) {
        CALreal density = 0;


        //pos particle
        const CALreal px = calGet3Dr(ca,Q.px[slot],i,j,k);
        const CALreal py = calGet3Dr(ca,Q.py[slot],i,j,k);
        const CALreal pz = calGet3Dr(ca,Q.pz[slot],i,j,k);
        VEC3 p1= {px,py,pz};

        for (int n=0; n<ca->sizeof_X; n++) {
            for(int slot1; slot1<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot1++) {

                const CALreal npx = calGetX3Dr(ca,Q.px[slot1],i,j,k,n);
                const CALreal npy = calGetX3Dr(ca,Q.py[slot1],i,j,k,n);
                const CALreal npz = calGetX3Dr(ca,Q.pz[slot1],i,j,k,n);
                VEC3 p2= {npx,npy,npz};

                //v1 and v2 does not refer to the very same particle
                if(!(n==0 && slot1==slot) && isNeigh(p1,p2)) {
                    VEC3 d;
                    diff(p1,p2,d);

                    CALreal dist2 = dot(d,d);
                    density += WPoly6(dist2,RADIUS);
                }
            }
        }
        density = density*MASS;

        const CALreal pressure = STIFFNESS* ( density - RESTDENSITY );

        calSet3Dr(ca,Q.density[slot],i,j,k,density);
        calSet3Dr(ca,Q.pressure[slot],i,j,k,pressure);
    }
}




void computePressureAcceleration(struct CALModel3D* ca, int i, int j, int k) {

    const VEC3 G = {0,-9.81,0};
    int cn=0;
    for(int slot; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++) {



        VEC3 f_gravity= {0,0,0};
        VEC3 f_pressure= {0,0,0};
        VEC3 f_viscosity= {0,0,0};
        VEC3 f_surface= {0,0,0};
        VEC3 a_external= {0,0,0};
        VEC3 colorFieldNormal= {0,0,0};
        CALreal colorFieldLaplacian;

        //pos particle
        const CALreal px = calGet3Dr(ca,Q.px[slot],i,j,k);
        const CALreal py = calGet3Dr(ca,Q.py[slot],i,j,k);
        const CALreal pz = calGet3Dr(ca,Q.pz[slot],i,j,k);
        VEC3 p1= {px,py,pz};


        //vel particle
        const CALreal vx = calGet3Dr(ca,Q.vx[slot],i,j,k);
        const CALreal vy = calGet3Dr(ca,Q.vy[slot],i,j,k);
        const CALreal vz = calGet3Dr(ca,Q.vz[slot],i,j,k);
        VEC3 v1= {vx,vy,vz};

        CALreal pressure = calGet3Dr(ca,Q.pressure[slot],i,j,k);
        CALreal density = calGet3Dr(ca,Q.density[slot],i,j,k);
        for (int n=0; n<ca->sizeof_X; n++) {
            for(int slot1; slot1<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot1++) {


                CALreal density = 2.0;

                const CALreal npx = calGetX3Dr(ca,Q.px[slot1],i,j,k,n);
                const CALreal npy = calGetX3Dr(ca,Q.py[slot1],i,j,k,n);
                const CALreal npz = calGetX3Dr(ca,Q.pz[slot1],i,j,k,n);
                VEC3 p2= {npx,npy,npz};


                //vel particle
                const CALreal vvx = calGetX3Dr(ca,Q.vx[slot],i,j,k,n);
                const CALreal vvy = calGetX3Dr(ca,Q.vy[slot],i,j,k,n);
                const CALreal vvz = calGetX3Dr(ca,Q.vz[slot],i,j,k,n);
                VEC3 v2= {vvx,vvy,vvz};


                if(isNeigh(p1,p2)) {

                    CALreal npressure = calGetX3Dr(ca,Q.pressure[slot1],i,j,k,n);
                    CALreal ndensity = calGetX3Dr(ca,Q.density[slot1],i,j,k,n);

                    VEC3 d;
                    diff(p1,p2,d);
                    CALreal dist_2 = dot(d,d);

                    VEC3 poly6Gradient;
                    Wpoly6Gradient(d, dist_2, poly6Gradient,RADIUS);

                    VEC3 spikyGradient;
                    WspikyGradient(d, dist_2, spikyGradient, RADIUS);

                    if(!(n==0 && slot1==slot)) {
                        VEC3 tmp;
                        scale(spikyGradient, pressure/pow(density,2)+
                              npressure/pow(ndensity,2),tmp );

                        add(f_pressure,tmp,f_pressure);


                        VEC3 vd;
                        diff(v1,v2,vd);
                        scale(vd , WviscosityLaplacian(dist_2,RADIUS) / ndensity , tmp );
                        add(f_viscosity, tmp, f_viscosity);

                    }

                    add(colorFieldNormal,poly6Gradient,colorFieldNormal);
                    scale(colorFieldNormal,1/ndensity,colorFieldNormal);
                    colorFieldLaplacian += Wpoly6Laplacian(dist_2,RADIUS) / ndensity;


                }
            }



        }//neight
        scale(f_pressure, -MASS* density, f_pressure);

        scale(f_viscosity, VISCOSITY * MASS, f_viscosity);

        scale(colorFieldNormal,MASS,colorFieldNormal);

        //NORMAL VECTORY GOES SOMEWHERE HERE
        //p->normal = -1.0 * colorFieldNormal;


        colorFieldLaplacian *=MASS;

        // surface tension force
        double colorFieldNormalMagnitude = magnitude(colorFieldNormal); //check that this is actually the magnitude of the vector
        
        if (colorFieldNormalMagnitude >  SURFACE_THRESHOLD) {

                calSet3Db(ca,Q.flag[slot],i,j,k,1);
                
                VEC3 tmp2;
                scale(colorFieldNormal, -SURFACE_TENSION , tmp2);
                scale(tmp2,1/colorFieldNormalMagnitude * colorFieldLaplacian, f_surface);
                //f_surface = - SURFACE_TENSION * colorFieldNormal / colorFieldNormalMagnitude * colorFieldLaplacian;

            }
            else {
                calSet3Db(ca,Q.flag[slot],i,j,k,0);
            }
        
         const CALreal eps = 0.1;
            VEC3 acc={0,0,0};
            // if(cn >1 && (p->density >= eps || p->density <= -eps))
            //p->acc = (f_pressure  +f_viscosity + f_surface + f_gravity ) / p->density;
            add(acc,f_pressure,acc);
            add(acc,f_viscosity,acc);
            add(acc,f_surface,acc);
            add(acc,f_gravity,acc);
            scale(acc,1/density,acc);
            
            

            // std::cout<<p->acc.x<<" "<<p->acc.y<<std::endl;

            /**add external forces, collision, user interaction, moving rigid bodys etc.**/
            //a_external = computeExternalFoces(p);
            add(acc,a_external,acc);
            
            //advance(i,j,k,slot,acc);
        
    }
}





