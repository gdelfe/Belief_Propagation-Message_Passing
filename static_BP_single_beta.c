// How to run:
// ./a.out N beta FileName
// Results are printed in BP_magn_b%s.dat

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>
//#include <cblas.h>

#define directory "DMP_data" //directory containing the
#define OPEN 0
#define CLOSE 1

#define Tmax 100000
#define EPSILON 1e-9
#define DAMP 0.9

struct variable{ //Structure for MC simulations
    
    // number of neighbours of a given site, incoming and outcoming ones
    int degree;
    int *neigh;         //neighbours of a given site
    
    double *J;     // couplings J_{ij}
};

/* -------------------------------------------  */
struct bp{ // Structure for dynamic message passing algorithm
    
    double *to; //value of the message from the site to another site
};
/* -------------------------------------------  */

void get_parameters(int argc, char **argv);
void allocate_memory(struct variable **site, struct bp **u, struct bp **temp);
void read_ERRG(struct variable *site, char **argv);
void initialize_bp(struct variable *site, struct bp *u,double bias);
void update_bp(struct variable *site, struct bp *u, struct bp *temp, double beta, int *t);
void magnetization(struct variable *site, struct bp *u, double beta,char *argv[]);
void verify_bp(struct variable *site, struct bp *u, double beta);
double get_rand(void);

int N;
double beta;

/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char *argv[]){
    
    chdir(directory); // Move to the directory
  
    struct variable *site;
    struct bp *u,*temp;
    int t;

    
    get_parameters(argc,argv); // get parameters from command line, stdin
    allocate_memory(&site,&u,&temp);
    read_ERRG(site,argv);
    
    initialize_bp(site,u,1.);
    
    for(t=0;t<Tmax;t++){
        update_bp(site,u,temp,beta,&t);
    }
   

    magnetization(site,u,beta,argv);
    verify_bp(site,u,beta);
    
    return 1;
}

/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */

// get parameter from stdout and convert them in the right format
void get_parameters(int argc, char **argv){
    
    if(argc>1){
        
        N=atoi(argv[1]);
        beta=atof(argv[2]);
    }
    
    return;
    
}


/*************************************************************************************/

void allocate_memory(struct variable **site, struct bp **u, struct bp **temp){
    
    int i;
    
    /* allocate memory for the array of structures */
    *site=(struct variable *)malloc(N*sizeof(struct variable));        //all.mem. for an array of structures
    *u=(struct bp *)malloc(N*sizeof(struct bp));        //all.mem. for an array of structures
    *temp=(struct bp *)malloc(N*sizeof(struct bp));        //all.mem. for an array of structures
    
    for(i=0;i<N;i++){
        
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours
        (*site)[i].J=(double *)malloc(N*sizeof(double));              //For each element of the array, all. mem. for its couplings
        (*site)[i].degree=0;                                    //set the initial degree equal to zero
        
        (*u)[i].to=(double *)malloc( N * sizeof(double));
        (*temp)[i].to=(double *)malloc( N * sizeof(double));

    }
    
    
    return;
    
}

/*************************************************************************************/

//Read the graph from a file generated with another program

void read_ERRG(struct variable *site, char **argv){
    
    int i,j,nn;
    double J1,J2;
    FILE *fp_file;
    
    char filename[101];
    
    sprintf(filename,"%s_MC_graph.dat",argv[3]);
    fp_file=fopen(filename,"r");
    
    if(fp_file==NULL){
        fprintf(stderr,"PROBLEM OPENING FILE %s\n\n" ,"XXX_MC_graph.dat");
        exit(errno);
    }
    
    for(i=0;i<N;i++){
        fscanf(fp_file,"%d",&site[i].degree);
  //      printf("degree %d = %d\n",i,site[i].degree);
        for(j=0;j<site[i].degree;j++){
            fscanf(fp_file,"%d%d%lf%lf",&i,&nn,&J1,&J2);
            site[i].neigh[j]=nn;
            site[i].J[nn]=J1;
            site[nn].J[i]=J2;
      //      printf("%d ---> %d \t J[%d][%d] = %lf \t J[%d][%d] = %lf\n",i,nn,i,nn,site[i].J[nn],nn,i,site[nn].J[i]);
        }
        
    }
    
    return;
}

/*************************************************************************************/

// Initialize all the initial fields to random values
void initialize_bp(struct variable *site, struct bp *u, double bias){

    int i,j,nn;
    for(i=1;i<N;i++){
        for(j=0;j<site[i].degree;j++){
            nn=site[i].neigh[j];
            u[i].to[nn]= (10 * get_rand() ) + bias;
        }
    }
    
    return;

}

/*************************************************************************************/

void update_bp(struct variable *site, struct bp *u, struct bp *temp, double beta, int *t){

    int i,j,k,nnj,nnk;
    double sum;
    
    double messChange,check;
    
    if(*t%1000==0) printf("Iterations = %d\n",*t);
    
    for(i=0;i<N;i++){ //for each i
        for(j=0;j<site[i].degree;j++){ // for each outgoing from i
            nnj=site[i].neigh[j];
            temp[i].to[nnj]=u[i].to[nnj]; //old message
            sum=0;
            for(k=0;k<site[i].degree;k++){ //sum over the incoming to i
                if(k!=j){
                    nnk=site[i].neigh[k];
                    sum += atanh(tanh(beta*site[nnk].J[i]) * tanh(beta*u[nnk].to[i]));
                }
            }
            
            u[i].to[nnj]=DAMP*(1/beta * sum)+(1.-DAMP)*temp[i].to[nnj];
        }
    }

    /* ---- CHECK DIFFERENCE ----- */
    
    check=0;
    for(i=0;i<N;i++){
        for(j=0;j<site[i].degree;j++){
            nnj=site[i].neigh[j];
            messChange=fabs(u[i].to[nnj]-temp[i].to[nnj]);
            if(messChange > EPSILON){
                check = 1;
            }
        }
    }
    if(check == 0){
        printf("Numb iterations to convergence = %d\n",*t);
        *t = Tmax;
    }
    
    
    return;
}

/*************************************************************************************/

// Compute the value of the magnetization
// M = \sum_i tanh( \beta * u_{cav_i})
// u_{cav_i} = 1/ \beta \sum_{j \in i}

void magnetization(struct variable *site, struct bp *u, double beta, char *argv[]){
    
    int i,j,nn;
    double u_cav,M;
    char filename[101];

    FILE *fp_m;
    
    M = 0.;
    for(i=0;i<N;i++){
        u_cav = 0.;
        for(j=0;j<site[i].degree;j++){ // sum over all the j neighbours
            nn = site[i].neigh[j];
            u_cav += atanh( tanh(beta*site[nn].J[i]) * tanh(beta * u[nn].to[i])); // sum of all the incoming contributions to node i
                                                            // NOTE: there should be a factor 1/beta but it cancels in the magn= tanh(beta*u_cav)
        }
        M += tanh(u_cav);
    }
    
    
    sprintf(filename,"BP_magn_%s_b%.2lf.dat",argv[3],beta);
    fp_m=fopen(filename,"w");
    fprintf(fp_m,"%lf\t%lf\n",beta,M/((double) N));
    
    return ;
    
}

/***************************************************************/

// OPTIONAL:
// Verify that the value of the converged fields satisfy the BP equation

void verify_bp(struct variable *site, struct bp *u, double beta){

    int i,j,k,nnj,nnk;
    double sum;
    FILE *fp_u;
    
    fp_u=fopen("Field_value.dat","w");
    
    for(i=0;i<N;i++){
        for(j=0;j<site[i].degree;j++){
            nnj=site[i].neigh[j];
            sum=0;
            for(k=0;k<site[i].degree;k++){
                nnk=site[i].neigh[k];
                if(nnk!=nnj){
                    sum += atanh( tanh(beta * site[nnk].J[i]) * tanh(beta * u[nnk].to[i]));
                }
            }
            fprintf(fp_u,"RHS = %lf\tLHS = %lf\tdiff = %lf\n",sum/beta, u[i].to[nnj], sum/beta-u[i].to[nnj]);
        }
        
    }

    return;
}

/*************************************************************************************/

// (Very simple) Random generator

double get_rand(void){
    
    return(1.0 * rand()/(1.0 + RAND_MAX));
}



