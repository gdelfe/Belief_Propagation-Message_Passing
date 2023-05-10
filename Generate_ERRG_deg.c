
// Generate a Erdős–Rényi Random Graph with a given average degree c, with the option of having: 
// 1. asymetric couplings, 2. positive and negative couplings
//  Generate_ERRG.c
// To run the code on linux: $ gcc Generate_ERRG_comm_line.c -lm, then: $./a.out N c J1 J11 J2 J22,
// where J1 and J11 are forward and backward link that connect to nodes i and j. If you want this connection to be
// symmetric and positive choose: J1 = A; J11 = A; with A positive number. You can also choose to have the direct connection
// positive and the indirect connection negative, s.t. J1 = A; J11 = -A. Alternative, you can have that the two connections have different
// magnitude: J1 = A; J11 = B.
// On average, 50% of the links are assigned through J1 and J11, the other 50% via J2 and J22. These latter variables have the same meaning
// they are used just to extend the number of options. 
//
// For instance if you choose J1 = 1, J11 = 1, J2 = -1, J22 = -1, you'll have 50% of the 
// connections which are direct and unitary positive, and 50% which are direct and unitary negative. 

// If you choose J1 = 1, J11 = 2, J2 = -1, J2 = -2, you'll have that 50% of the connections are positive, and asymmetric, and 50% of the 
// connections are negative and asymmetric. 
// 
// To run: 
// ./a.out N c J1 J11 J2 J22, name_file; where N =  size of the network, c = average connectivity, name_file = name of file output. 
// For instance:
// $./a.out 1000 3 1 1 -1 -1 output_file
//
//  Created by Gino Del Ferraro on 1/20/15.
//
// OUTPUT: 
// output_file_ERRgraph.dat -> file to plot the graph using grapviz
// output_file_MC_graph.dat -> file to be loaded for MC simulations or BP simulations or other. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define directory "DMP_data" //directory containing the data


struct variable{ //Structure for MC simulations
    int degree; // number of neighbours of a given site
    int *neigh; //neighbours of a given site
    double *J;     // couplings J_{ij}
};

// get parameters from stdout ------
void get_parameter(int argc, char **argv);

void allocate_memory(struct variable **site); // allocate memory
void generate_ERRG(struct variable *site); // generate ERRG WITHOUT isolated nodes

void print_graph_on_file(struct variable *site, char **);   // print graph on 2 files
// Rand numb generator
double get_rand(void);

int N,c;  // N spin and average connectivity c
double J1,J11,J2,J22; // couplings


/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char *argv[]){
    
    struct variable *site;
    int i,j;

    chdir(directory); // Move to the directory
    
    srand(time(NULL)); //initialize get_rand(), random numb generator
    
    get_parameter(argc,argv);
    allocate_memory(&site);
    generate_ERRG(site);
    print_graph_on_file(site,argv);
    

    return 0;
}




/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */


void get_parameter(int argc, char **argv){
    
    if(argc>1){
        
        N=atoi(argv[1]);
        c=atoi(argv[2]);
        J1=atof(argv[3]);
        J11=atof(argv[4]);
        J2=atof(argv[5]);
        J22=atof(argv[6]);
        
    }

    return;

}

/***************************************************/

void allocate_memory(struct variable **site){
    
    int i;
    
    /* allocate memory for the array of structures */
    *site=(struct variable *)malloc(N*sizeof(struct variable));        //all.mem. for an array of structures
    for(i=0;i<N;i++){
        
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours
        (*site)[i].J=(double *)malloc(N*sizeof(double));              //For each element of the array, all. mem. for its couplings
        (*site)[i].degree=0;                                    //set the initial degree equal to zero
    }
    
    return;
    
}

/***************************************************/

void generate_ERRG(struct variable *site){
    
    int i,j,nodes=0;
    double z;
    
    
    printf("\nER random graph, N = %d, c = %d, ",N,c);
    
    z=((double)c/(double)N); //average connectivity of the graph
    
    for(i=0;i<N;i++){               //note: when calling this function, site[i].degree=0 initially
        for(j=i+1;j<N;j++){
            
            if(get_rand()<z){  // get_rand() gives a random numb between [0,1]
                
                site[i].neigh[site[i].degree]=j; //connect i to j
                site[i].degree++;
                site[j].neigh[site[j].degree]=i; // connect j to i
                site[j].degree++;
                
                if(get_rand()<0.5){ // generate random symmetric couplings J_{ij}=\pm 1
                    site[i].J[j]=(double)J1/c; // J_{ij}
                    site[j].J[i]=(double)J11/c; // J_{ji}
                }
                else{
                    site[i].J[j]=(double)J2/c;  // J_{ij}
                    site[j].J[i]=(double)J22/c;  // J_{ji}
                }
            }
    
        }
        if(site[i].degree==0){  // go back and recompute the nodes if it was isolated
            i--; 
        }
        
    }
    // No isolated nodes at the end of this process
    
    printf("Numb isolated nodes = %d\n\n",nodes);
    return;
}



/***************************************************/

// Print graph on the file ERRgraph.dat to be plotted using graphics and on MC_graph.dat
// to be read by the C program MC_DMP.c which runs MC simulations based on the same graph

void print_graph_on_file(struct variable *site, char **argv){
    
    int i,j;
    
    FILE *fp_graph, *fp_MC, *fp_Amat;
    char filename[101];
    
    sprintf(filename,"%s_ERRgraph.dat",argv[7]);  // file to plot the graph using grapviz
    fp_graph=fopen(filename,"w");
    sprintf(filename,"%s_MC_graph.dat",argv[7]);  // file with graph to be used in MC and other simulations
    fp_MC=fopen(filename,"w");
    
    fprintf(fp_graph,"digraph {\n        graph [rankdir=LR];\n");  //instructions for graphviz
    for(i=0;i<N;i++){
        fprintf(fp_MC,"%d\n",site[i].degree);
        if(site[i].degree==0){
            fprintf(fp_graph,"        %d;\n",i);
        }
        for(j=0;j<site[i].degree;j++){
            if(site[i].J[site[i].neigh[j]]>0.)
                fprintf(fp_graph,"        %d -> %d[color=\"red\"];\n",i,site[i].neigh[j]);  //instructions for graphviz
            if(site[i].J[site[i].neigh[j]]<0.)
                fprintf(fp_graph,"        %d -> %d[color=\"blue\"];\n",i,site[i].neigh[j]);  //instructions for graphviz
                
            fprintf(fp_MC,"%d\t%d\t%lf\t%lf\n",i,site[i].neigh[j],site[i].J[site[i].neigh[j]],site[site[i].neigh[j]].J[i]); // print i --> i-neigh, J[i][nn] J[nn][i]
        }
    }
    fprintf(fp_graph,"}\n");
    
    
    fclose(fp_graph);
    fclose(fp_MC);
    return;
}

/***************************************************/



double get_rand(void){
    
    return(1.0 * rand()/(1.0 + RAND_MAX));
}










