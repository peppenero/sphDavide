#include "modello.h"
#include "parser.h"
#include "sph.h"
#include<fstream>
#include<iostream>
#include <cstring>
#include <stdlib.h>

using namespace std;


struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;

void avanza(struct CALModel3D* ca, int i, int j, int k){

        VEC3r acc(0,0,0);
	for(int n=0;n<MAX_NUMBER_OF_PARTICLES_PER_CELL;n++){
		if(ncestiFluiduNtraStuSlot(ca,i,j,k,n)){
                        calcolaDensita(ca,i,j,k);
			calUpdate3D(ca);
                        computePressureAcceleration(ca,i,j,k);
		}
	}
}


//Ritorna true se almeno uno tra gli slot di tutte le celle vicine contiene
//almeno una particella di fluido.
CALbyte ncestiArmenuNaParticellaNtraIVicini(struct CALModel3D* ca, int i, int j, int k, int n)
{
	for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
		if (calGetX3Di(ca, Q.imove[slot],i,j,k,n) == PARTICLE_PRESENT)
			return CAL_TRUE;

	return CAL_FALSE;
}

CALbyte ncestiFluiduNtraStuSlot(struct CALModel3D* ca, int i, int j, int k,int slot){
	if(calGet3Di(ca,Q.imove[slot],i,j,k)==1)
		return CAL_TRUE;
	return CAL_FALSE;
}

//Funzione di prova per simulare il movimento sull'asse Y.
void movili(struct CALModel3D* ca, int i, int j, int k)
{
	if (!ncestiArmenuNaParticellaNtraIVicini(ca, i, j, k, 0))
		return;

	CALreal delta_z = -CL/3;
	CALreal z;
	CALreal z_new;

	for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
		if (calGet3Di(ca, Q.imove[c],i,j,k) == PARTICLE_PRESENT)
		{
			z = calGet3Dr(ca, Q.pz[c],i,j,k);
			z_new = z + delta_z;
			calSet3Dr(ca, Q.pz[c],i,j,k,z_new);

			//            CALreal z_new_test = calGetNext3Dr(ca,Q.pz[c],i,j,k);
			//            int a = 0;
		}
}

//Funzione che dato uno slot appartenente alla cella (i,j,k) lo libera.
void pezziala(int slot, struct CALModel3D* ca, int i, int j, int k)
{
	calSet3Dr(ca, Q.px[slot],   i,j,k,NODATA);
	calSet3Dr(ca, Q.py[slot],   i,j,k,NODATA);
	calSet3Dr(ca, Q.pz[slot],   i,j,k,NODATA);
	calSet3Dr(ca, Q.vx[slot],   i,j,k,0);
	calSet3Dr(ca, Q.vy[slot],   i,j,k,0);
	calSet3Dr(ca, Q.vz[slot],   i,j,k,0);
	calSet3Dr(ca, Q.density[slot],   i,j,k,0);
	calSet3Dr(ca, Q.pressure[slot],   i,j,k,0);
	calSet3Di(ca, Q.imove[slot],i,j,k,PARTICLE_ABSENT);
}

//Funzione che dato uno slot libero muove le particelle che hanno cambiato cella.
void sucala(int slot, int e, struct CALModel3D* ca, int i, int j, int k, int n)
{
	calSet3Dr(ca, Q.px[slot],i,j,k,calGetX3Dr(ca,Q.px[e],i,j,k,n));
	calSet3Dr(ca, Q.py[slot],i,j,k,calGetX3Dr(ca,Q.py[e],i,j,k,n));
	calSet3Dr(ca, Q.pz[slot],i,j,k,calGetX3Dr(ca,Q.pz[e],i,j,k,n));
	calSet3Dr(ca, Q.vx[slot],i,j,k,calGetX3Dr(ca,Q.vx[e],i,j,k,n));
	calSet3Dr(ca, Q.vy[slot],i,j,k,calGetX3Dr(ca,Q.vy[e],i,j,k,n));
	calSet3Dr(ca, Q.vz[slot],i,j,k,calGetX3Dr(ca,Q.vz[e],i,j,k,n));
	calSet3Dr(ca, Q.density[slot],i,j,k,calGetX3Dr(ca,Q.density[e],i,j,k,n));
	calSet3Dr(ca, Q.pressure[slot],i,j,k,calGetX3Dr(ca,Q.pressure[e],i,j,k,n));
	calSet3Di(ca, Q.imove[slot],i,j,k,calGetX3Di(ca,Q.imove[e],i,j,k,n));
}

//Funzione che gestisce il movimento delle particelle tramite la chiamata delle
//funzioni pezziala e sucala
void moviliCazzu(struct CALModel3D* ca, int i, int j, int k)
{
	//    CALreal z_min = k*CL;
	//    CALreal z_max = (k+1)*CL;
	//    CALreal z_min = (SLICES-1 - (k+1))*CL;
	//    CALreal z_max = (SLICES-1 - k)*CL;

	CALreal x,  y,  z;
	CALint _i, _j, _k;

	if (ncestiArmenuNaParticellaNtraIVicini(ca, i, j, k, 0))
	{
		//pezziali
		for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
			if (calGet3Di(ca, Q.imove[slot],i,j,k) == PARTICLE_PRESENT)
			{
				x = calGet3Dr(ca, Q.px[slot],i,j,k);
				y = calGet3Dr(ca, Q.py[slot],i,j,k);
				z = calGet3Dr(ca, Q.pz[slot],i,j,k);

				_i = x/CL;
				_j = y/CL;
				_k = z/CL;

				if ((i != _i) || (j != _j) || (k != _k))
					pezziala(slot, ca,i,j,k);
			}
	}

	//sucali
	for (int n=1; n<ca->sizeof_X; n++)
		for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
			if (calGetX3Di(ca,Q.imove[slot],i,j,k,n) == PARTICLE_PRESENT)
			{
				x = calGetX3Dr(ca, Q.px[slot],i,j,k,n);
				y = calGetX3Dr(ca, Q.py[slot],i,j,k,n);
				z = calGetX3Dr(ca, Q.pz[slot],i,j,k,n);

				_i = x/CL;
				_j = y/CL;
				_k = z/CL;

				int a;
				if (i==25 && j==25 && k==24 && n==18 && slot==0)
					a = 0;

				if ((i == _i) && (j == _j) && (k == _k))
				{
					int c;
					for (c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
						if (calGetNext3Di(ca,Q.imove[c],i,j,k) == PARTICLE_ABSENT)
						{
							slot = c;
							break;
						}

					if (c < MAX_NUMBER_OF_PARTICLES_PER_CELL)
						sucala(c,slot,ca,i,j,k,n);
				}
			}
}


//Funzione che date le posizione reali della particella le converte in indici
//della matrice dell'automa cellulare ed imposta i valori iniziali.
void danciNaPosizioni(struct CALModel3D* ca, const CALreal x, const CALreal y, const CALreal z,const CALreal nx, const CALreal ny, const CALreal nz, const CALreal rho, const CALint imove)
{
	int i = x/CL;
	int j = y/CL;
	int k = z/CL;
	//int k = (SLICES -1) -  z/CL;

	CALint slot = -1;
	for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
		if (calGet3Di(ca,Q.imove[c],i,j,k) == PARTICLE_ABSENT)
		{
			slot = c;
			break;
		}

	if(slot != -1)
	{
		calInit3Dr(ca,Q.px[slot],   i,j,k,x	);
		calInit3Dr(ca,Q.py[slot],   i,j,k,y	);
		calInit3Dr(ca,Q.pz[slot],   i,j,k,z	);
		calInit3Dr(ca,Q.nx[slot],	i,j,k,nx);
		calInit3Dr(ca,Q.ny[slot],	i,j,k,ny);
		calInit3Dr(ca,Q.nz[slot],	i,j,k,nz);
                calInit3Dr(ca,Q.density[slot],	i,j,k,rho);
		calInit3Di(ca,Q.imove[slot],i,j,k,imove);
	}
}

//Funzione che fa muovere e fermare l'automa.
CALbyte caminalu(struct CALModel3D* modello)
{
	CALint step = a_simulazioni->step;

	if (step <= STEPS)
		return CAL_FALSE;

	return CAL_TRUE;
}

//Funzione di prova
void stampaAPosizioni(struct CALModel3D* ca, int i, int j, int k)
{
	CALreal z = 0;
	if (k == 10)
		for (int c=0; c<MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
		{
			z = calGet3Dr(ca,Q.pz[0],i,j,k);
			if (z != NODATA)
                                printf ("[(i,j,k),[c],z]=[(%d,%d,%d),[%d],%f]\n", i,j,k,c,z);
		}
}

//Funzione che lascia solo le particelle di boundary
void sbacanta(struct CALModel3D* ca, int i, int j, int k)
{
	if (ncestiArmenuNaParticellaNtraIVicini(ca, i, j, k, 0))
		for (int slot = 0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
			pezziala(slot, ca, i, j, k);
}

//Funzione di transizione globale
void transizioniGlobali(struct CALModel3D* modello)
{
	//calApplyElementaryProcess3D(modello,movili);
	//calUpdate3D(modello);

	//    printf("z(25,25,25) = %f\n", calGet3Dr(modello,Q.pz[0],25,25,25));
	//    for (int n=0; n<modello->sizeof_X; n++)
	//        printf("z(25,25,24; %d) = %f\n", n, calGetX3Dr(modello,Q.pz[0],25,25,24,n));

	//calApplyElementaryProcess3D(modello,moviliCazzu);
	//calUpdate3D(modello);

	calApplyElementaryProcess3D(modello,avanza);
	calUpdate3D(modello);

	calApplyElementaryProcess3D(modello,moviliCazzu);
	calUpdate3D(modello);

}

void danciNaPosizioni1(struct CALModel3D* ca, const CALreal x, const CALreal y, const CALreal z, const CALint imove){


    int i = x/CL1;
    int j = y/CL1;
    int k = z/CL1;

    //int k = (SLICES -1) -  z/CL;

    CALint slot = -1;
    for (int c = 0; c < MAX_NUMBER_OF_PARTICLES_PER_CELL; c++)
            if (calGet3Di(ca,Q.imove[c],i,j,k) == PARTICLE_ABSENT)
            {
                    slot = c;
                    break;
            }

     if(slot != -1)
     {
           calInit3Dr(ca,Q.px[slot],   i,j,k,x	);
           calInit3Dr(ca,Q.py[slot],   i,j,k,y	);
           calInit3Dr(ca,Q.pz[slot],   i,j,k,z	);
           calInit3Di(ca,Q.imove[slot],i,j,k,imove);
     }

}

int leggiFile(){

    FILE *posizioni;

    char line[100];
    char *token;
    int n=0;
    double s[3];

    posizioni = fopen("posizioni.txt","r");

    while (fgets(line, 100, posizioni) != NULL) {
                    int campo=0;

                    token = strtok(line," ");

                    while(token != NULL){
                            switch (campo) {
                            case 0:
                                    //printf("%f\n",atof(token));
                                    s[0]=atof(token);
                                    //particella.posizione[0]=atof(token);
                                    campo++;
                                    break;
                            case 1:

                                    token = strtok(NULL, " ");
                                    s[1]=atof(token);
                                    //printf("%f\n",atof(token));
                                    //particella.posizione[1]=atof(token);
                                    campo++;
                                    break;
                            case 2:

                                    token = strtok(NULL, " ");
                                    s[2]=atof(token);
                                    //printf("%f\n",atof(token));
                                    //particella.posizione[2]=atof(token);
                                    campo++;
                                    break;

                            case 3:

                                    token = strtok(NULL, " ");
                                    campo=0;
                                    break;
                            }
                    }

                        danciNaPosizioni1(u_modellu, s[0]+0.4,s[1]+0.4,s[2]+0.4,1);

    }
}


//Funzione d'inizializzazione e di start dell'automa
void partilu()
{
	u_modellu = calCADef3D(ROWS,COLS,SLICES,CAL_MOORE_NEIGHBORHOOD_3D,CAL_SPACE_TOROIDAL,CAL_NO_OPT);

	for(int i=0;i<MAX_NUMBER_OF_PARTICLES_PER_CELL;i++)
	{
		Q.px[i]    = calAddSubstate3Dr(u_modellu);
		Q.py[i]    = calAddSubstate3Dr(u_modellu);
		Q.pz[i]    = calAddSubstate3Dr(u_modellu);
		Q.vx[i]    = calAddSubstate3Dr(u_modellu);
		Q.vy[i]    = calAddSubstate3Dr(u_modellu);
		Q.vz[i]    = calAddSubstate3Dr(u_modellu);
		Q.nx[i]    = calAddSubstate3Dr(u_modellu);
		Q.ny[i]    = calAddSubstate3Dr(u_modellu);
		Q.nz[i]    = calAddSubstate3Dr(u_modellu);
		Q.density[i] = calAddSubstate3Dr(u_modellu);
		Q.pressure[i] = calAddSubstate3Dr(u_modellu);
		Q.imove[i] = calAddSubstate3Di(u_modellu);
		Q.flag[i]    = calAddSubstate3Db(u_modellu);

		calInitSubstate3Dr(u_modellu,	Q.px[i],		NODATA);
		calInitSubstate3Dr(u_modellu,	Q.py[i],		NODATA);
		calInitSubstate3Dr(u_modellu,	Q.pz[i],		NODATA);
		calInitSubstate3Dr(u_modellu,	Q.vx[i],		0);
		calInitSubstate3Dr(u_modellu,	Q.vy[i],		0);
		calInitSubstate3Dr(u_modellu,	Q.vz[i],		0);
		calInitSubstate3Dr(u_modellu,	Q.nx[i],		0);
		calInitSubstate3Dr(u_modellu,	Q.ny[i],		0);
		calInitSubstate3Dr(u_modellu,	Q.nz[i],		0);
		calInitSubstate3Dr(u_modellu,	Q.density[i],	0);
		calInitSubstate3Dr(u_modellu,	Q.pressure[i],	0);
		calInitSubstate3Di(u_modellu,	Q.imove[i],		PARTICLE_ABSENT);
		calInitSubstate3Db(u_modellu,	Q.flag[i],		0);
	}

	particellaT particelle[NUMBER_OF_PARTICLES];

	parse(particelle);

	for(int i=0;i<NUMBER_OF_PARTICLES;i++)
		danciNaPosizioni(u_modellu, particelle[i].posizione[0],	particelle[i].posizione[1],	particelle[i].posizione[2],
                                        particelle[i].normale[0],	particelle[i].normale[1],	particelle[i].normale[2],
                                        particelle[i].rho,              particelle[i].imove);

	// calApplyElementaryProcess3D(modello, printPos);

//        calApplyElementaryProcess3D(u_modellu, sbacanta);
//        calUpdate3D(u_modellu);
//



      //  danciNaPosizioni(u_modellu,0.025,0.025,0.025,0,0,0,1000.68373876,1);
      //  calUpdate3D(u_modellu);

//        leggiFile();

        calUpdate3D(u_modellu);


	a_simulazioni = calRunDef3D(u_modellu,1,CAL_RUN_LOOP,CAL_UPDATE_IMPLICIT);

	calRunAddGlobalTransitionFunc3D(a_simulazioni, transizioniGlobali);

	calRunAddStopConditionFunc3D(a_simulazioni, caminalu);
}

