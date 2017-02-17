#include "modello.h"
#include "parser.h"
#include "sph.h"
#include <cstring>


struct CALModel3D* u_modellu;
struct CALRun3D* a_simulazioni;
struct Substates Q;

constexpr const double II=0.205010;
constexpr const double JJ=0.395000;
constexpr const double KK=0.265000;

void computeDensityElementaryProcess(struct CALModel3D* ca, int i, int j, int k){
    for(int n=0;n<MAX_NUMBER_OF_PARTICLES_PER_CELL;n++){
        if(ncestiFluiduNtraStuSlot(ca,i,j,k,n)){
            calcolaDensita(ca,i,j,k);
            break; //messo perche senza chiamerebbe le due funzioni sopra più di una volta
        }
    }
}



void stampaAltroStepElementaryProcess(struct CALModel3D* ca, int i, int j, int k){
    //printf("%i",a_simulazioni->step);
    if(a_simulazioni->step==2){
        for(int a=0;a<MAX_NUMBER_OF_PARTICLES_PER_CELL;a++){
            //if(calGet3Di(ca,Q.imove[a],i,j,k)==PARTICLE_PRESENT){

            printf("density-----,%f\n",calGet3Dr(ca,Q.density[a],i,j,k));
            printf("pressure-----,%f\n",calGet3Dr(ca,Q.pressure[a],i,j,k));

        }
    }
}



void stampaturiElementaryProcess(struct CALModel3D* ca, int i, int j, int k){
    if(a_simulazioni->step<=10){
        for(int a=0;a<MAX_NUMBER_OF_PARTICLES_PER_CELL;a++){
            if(calGet3Di(ca,Q.imove[a],i,j,k)==PARTICLE_PRESENT){
                printf("%f %f %f %f %f %f %f %f %f %f %f %f\n",calGet3Dr(ca,Q.px[a],i,j,k),calGet3Dr(ca,Q.py[a],i,j,k),calGet3Dr(ca,Q.pz[a],i,j,k),calGet3Dr(ca,Q.vx[a],i,j,k),calGet3Dr(ca,Q.vy[a],i,j,k),calGet3Dr(ca,Q.vz[a],i,j,k),calGet3Dr(ca,Q.ax[a],i,j,k),calGet3Dr(ca,Q.ay[a],i,j,k),calGet3Dr(ca,Q.az[a],i,j,k),calGet3Dr(ca,Q.mass[a],i,j,k),calGet3Dr(ca,Q.density[a],i,j,k),calGet3Dr(ca,Q.pressure[a],i,j,k));
            }
        }
    }
    //printf("USCITO--------------------------------------\n");
}



void computePressureAndAccelerationElementaryProcess(struct CALModel3D* ca, int i, int j, int k){
    for(int n=0;n<MAX_NUMBER_OF_PARTICLES_PER_CELL;n++){
        if(ncestiFluiduNtraStuSlot(ca,i,j,k,n)){
            computePressureAcceleration(ca,i,j,k);
            break; //messo perche senza chiamerebbe le due funzioni sopra più di una volta
        }
    }

}

//Funzione che lascia solo le particelle di boundary
void sbacantaElementaryProcess(struct CALModel3D* ca, int i, int j, int k)
{
    if (ncestiArmenuNaParticellaNtraIVicini(ca, i, j, k, 0))
        for (int slot = 0; slot<MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++)
            pezziala(slot, ca, i, j, k);
}

//Funzione di prova per simulare il movimento sull'asse Y.
void moviliElementaryProcess(struct CALModel3D* ca, int i, int j, int k)
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
        }
}


void stampaPosElementaryProcess(struct CALModel3D* ca, int i, int j, int k){
    if(a_simulazioni->step == 2){
        for (int slot = 0; slot < MAX_NUMBER_OF_PARTICLES_PER_CELL; slot++) {
            if(calGet3Di(ca,Q.imove[slot],i,j,k)==PARTICLE_PRESENT){
                printf("%f,%f,%f\n",calGet3Dr(ca,Q.px[slot],i,j,k),calGet3Dr(ca,Q.py[slot],i,j,k),calGet3Dr(ca,Q.pz[slot],i,j,k));
            }
        }
    }
}


//Funzione di transizione globale
void transizioniGlobali(struct CALModel3D* modello)
{
    if(stampa){
        printf("------------pos------------|------------vel------------|------------acc------------|--------norm--------|mass|pradius|density|pressure|\n");
       // calApplyElementaryProcess3D(modello,stampaturiElementaryProcess);
        printf("******************************************STEP******************************************\n");
    }

    calApplyElementaryProcess3D(modello,computeDensityElementaryProcess);
    calUpdate3D(modello);

    //calApplyElementaryProcess3D(modello,stampaAltroStepElementaryProcess);

    calApplyElementaryProcess3D(modello,computePressureAndAccelerationElementaryProcess);
    calUpdate3D(modello);


    //calApplyElementaryProcess3D(modello,stampaPosElementaryProcess);

    calApplyElementaryProcess3D(modello,moviliCazzu);
    calUpdate3D(modello);
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







//Funzione che dato uno slot appartenente alla cella (i,j,k) lo libera.
void pezziala(int slot, struct CALModel3D* ca, int i, int j, int k)
{
    calSet3Dr(ca, Q.px[slot],   i,j,k,NODATA);
    calSet3Dr(ca, Q.py[slot],   i,j,k,NODATA);
    calSet3Dr(ca, Q.pz[slot],   i,j,k,NODATA);
    calSet3Dr(ca, Q.vx[slot],   i,j,k,0);
    calSet3Dr(ca, Q.vy[slot],   i,j,k,0);
    calSet3Dr(ca, Q.vz[slot],   i,j,k,0);
    calSet3Dr(ca, Q.ax[slot],   i,j,k,0);
    calSet3Dr(ca, Q.ay[slot],   i,j,k,0);
    calSet3Dr(ca, Q.az[slot],   i,j,k,0);
    calSet3Dr(ca, Q.density[slot],   i,j,k,0);
    calSet3Dr(ca, Q.pressure[slot],   i,j,k,0);
    calSet3Dr(ca, Q.mass[slot], i,j,k,0);
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
    calSet3Dr(ca, Q.ax[slot],i,j,k,calGetX3Dr(ca,Q.vx[e],i,j,k,n));
    calSet3Dr(ca, Q.ay[slot],i,j,k,calGetX3Dr(ca,Q.vy[e],i,j,k,n));
    calSet3Dr(ca, Q.az[slot],i,j,k,calGetX3Dr(ca,Q.vz[e],i,j,k,n));
    calSet3Dr(ca, Q.density[slot],i,j,k,calGetX3Dr(ca,Q.density[e],i,j,k,n));
    calSet3Dr(ca, Q.pressure[slot],i,j,k,calGetX3Dr(ca,Q.pressure[e],i,j,k,n));
    calSet3Dr(ca, Q.mass[slot],i,j,k,calGetX3Dr(ca,Q.mass[e],i,j,k,n));
    calSet3Di(ca, Q.imove[slot],i,j,k,calGetX3Di(ca,Q.imove[e],i,j,k,n));
}



//Funzione che gestisce il movimento delle particelle tramite la chiamata delle
//funzioni pezziala e sucala
void moviliCazzu(struct CALModel3D* ca, int i, int j, int k)
{

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
                //PERICOLO----XXXXXXXXXX
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
                //PERICOLO----XXXXXXXXXX
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
    fprintf(stderr,"%d step\n", a_simulazioni->step);
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


void danciNaPosizioni1(struct CALModel3D* ca, const CALreal x, const CALreal y, const CALreal z, const CALint imove){

    int i = x/CL;
    int j = y/CL;
    int k = z/CL;

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
        calInit3Dr(ca,Q.mass[slot], i,j,k,MASS);
        calInit3Di(ca,Q.imove[slot],i,j,k,imove);
    }

}

int initParticelle(){
    double px =0;
    double py =0;
    double pz =0;
    for(int i=ROWS/4; i < 3*ROWS/4; i++){
        for (int j = COLS/4; j < 3*COLS/4; ++j) {
            for (int k = SLICES/2-1; k <= SLICES/2+1; ++k) {



                px = CL/2 + (i)*CL;
                py = CL/2 + (j)*CL;
                pz = CL/2 + (k)*CL;


                danciNaPosizioni1(u_modellu, px,py,pz,1);
            }
        }
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
        n++;
        while(token != NULL){
            switch (campo) {
            case 0:
                s[0]=atof(token);
                campo++;
                break;
            case 1:
                token = strtok(NULL, " ");
                s[1]=atof(token);
                campo++;
                break;
            case 2:
                token = strtok(NULL, " ");
                s[2]=atof(token);
                campo++;
                break;
            case 3:
                token = strtok(NULL, " ");
                campo=0;
                break;
            }
        }

        danciNaPosizioni1(u_modellu, s[0],s[1],s[2],1);

    }

    calApplyElementaryProcess3D(u_modellu,stampaturiElementaryProcess);
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

        Q.ax[i]    = calAddSubstate3Dr(u_modellu);
        Q.ay[i]    = calAddSubstate3Dr(u_modellu);
        Q.az[i]    = calAddSubstate3Dr(u_modellu);

        Q.nx[i]    = calAddSubstate3Dr(u_modellu);
        Q.ny[i]    = calAddSubstate3Dr(u_modellu);
        Q.nz[i]    = calAddSubstate3Dr(u_modellu);

        Q.mass[i] = calAddSubstate3Dr(u_modellu);
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

        calInitSubstate3Dr(u_modellu,	Q.ax[i],		0);
        calInitSubstate3Dr(u_modellu,	Q.ay[i],		0);
        calInitSubstate3Dr(u_modellu,	Q.az[i],		0);

        calInitSubstate3Dr(u_modellu,	Q.nx[i],		0);
        calInitSubstate3Dr(u_modellu,	Q.ny[i],		0);
        calInitSubstate3Dr(u_modellu,	Q.nz[i],		0);

        calInitSubstate3Dr(u_modellu,	Q.mass[i],		0);
        calInitSubstate3Dr(u_modellu,	Q.density[i],   	0);
        calInitSubstate3Dr(u_modellu,	Q.pressure[i],          0);
        calInitSubstate3Di(u_modellu,	Q.imove[i],		PARTICLE_ABSENT);
        calInitSubstate3Db(u_modellu,	Q.flag[i],		0);
    }


    particellaT particelle[NUMBER_OF_PARTICLES];
    //    parse(particelle);
    //    for(int i=0;i<NUMBER_OF_PARTICLES;i++)
    //        danciNaPosizioni(u_modellu, particelle[i].posizione[0],	particelle[i].posizione[1],	particelle[i].posizione[2],
    //                particelle[i].normale[0],	particelle[i].normale[1],	particelle[i].normale[2],
    //                particelle[i].rho,              particelle[i].imove);



    //    calApplyElementaryProcess3D(u_modellu, sbacantaElementaryProcess);
    //  calUpdate3D(u_modellu);
    //
    //    //leggiFile();
    initParticelle();
    calUpdate3D(u_modellu);


    a_simulazioni = calRunDef3D(u_modellu,1,1,CAL_UPDATE_EXPLICIT);
    calRunAddGlobalTransitionFunc3D(a_simulazioni, transizioniGlobali);
    calRunAddStopConditionFunc3D(a_simulazioni, caminalu);
}

