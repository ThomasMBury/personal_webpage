/*******************************************************************************
 
 Simulates an SIRV model coupled to social contagion on a multiplex network
 
 October 2015
 
 Tom  Bury, Chris Bauch
 
 ********************************************************************************/

#define MAX_TS 1000
/* max number of timepoints in the output timeseries */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

/* Construct a function to return a random number between 0 and 1 */
double kyrand()
{
    int maxnumber = 2147483647l;   /* 2^31 - 1 */
    double answer;
    answer=((double)random())/maxnumber;
    return(answer);
}

/* Construct a function to return the maximum of two numbers */

float max(x,y)
float x,y;
{
    float ans;
    
    if( x > y ){
        ans = x;
    }
    
    else{
        ans = y;
    }
    
    return ans;
}

/* MAIN */

int main()
{
    
    int i,j,k,m,p,q,r;				/* assorted counters */
    
    int sus_count;                  /* susceptible counter */
    int inf_count;                  /* infected counter */
    int rec_count;                  /* removed counter */
    int vac_count;                  /* vaccinated counter */
    int pos_count;                  /* positive view counter */
    int neg_count;                  /* negative view counter */
    
    int count_real;                 /* keeps count of the number of realizations */
    float sigma;                    /* Proportion of nodes that take the virtual neighbours of a randomly chosen node */
    

    float beta = 0.2;               /* disease transmission rate */
    float beta_plus = 0.1/4;          /* positive information/view transmission rate */
    float beta_minus = 0.1/4;         /* negative information/view transmission rate */
    float mu = 0.1/2;               /* recovery rate */
    
    int n_x;                        /* x coordinate of physical neighbour */
    int n_y;                        /* y coordinate of physical neighbour */
    int v_x;                        /* x coordinate of virtual neighbour */
    int v_y;                        /* y coordinate of virtual neighbour */
    
    
    
    /* this is a structure which holds information about the status of each node and the identity of its neighbours */
    /* a node is specified by it's x/y integer coordinates */
    typedef struct node_data
    {
        int neighbour_x[4];                /* identity of neighbours on the physical layer*/
        int neighbour_y[4];
        int neighbour_virtual_x[4];        /* identity of neighbours in the virtual layer */
        int neighbour_virtual_y[4];
        char status[2];                    /* infection and awareness status */
    } node_structure;                      /* structure for specifications of individual */
    
    /* file declarations */
    FILE *fout1a;
    FILE *fout1b;
    FILE *fout1c;
    FILE *fout2a;
    FILE *fout2b;
    FILE *fout3a;
    FILE *fout3b;
    
    
    int N = 50;                 /* length of one side of the square lattice. population size = N times N */
    int SAMPLE_INCR =  1500;  	/* frequency with which lattice state is sampled for output files (greater value=less frequent sampling) */
    int MAX_ITER = 600000;      /* number of iterations */
    int NUM_REALIZATIONS = 10;  /* number of simulation runs */
    int INIT_inf = (int)(N*N/25); 	/* initial number of infected nodes */
    int INIT_pos = (int)(N*N/10);   /* initial number of nodes with positive vaccination views */
    int INIT_neg = (int)(N*N/10);  /* initial number of nodes with negative vaccingation views */
                                     
    node_structure node[N][N]; 	/* data structure for the nodes */
    
    float time_array[MAX_TS];       /* stores the timeseries of time points */
    float S_timeseries[MAX_TS];     /* stores the timeseries of number susceptible */
    float I_timeseries[MAX_TS];     /* stores the timeseries of number infected */
    float R_timeseries[MAX_TS];     /* stores the timeseries of number recovered */
    float V_timeseries[MAX_TS];     /* stores the timeseries of number vaccinated */
    float pos_timeseries[MAX_TS];   /* stores the timeseries of number of positive nodes */
    float neg_timeseries[MAX_TS];   /* stores the timeseries of number of negative nodes */
    int ts_count;                   /* keeps track of number of timesteps accumulated */
    
    /* output files */
    fout1a = fopen("ts.txt","w");               /* timeseries sigma = 0 */
    fout1b = fopen("ts2.txt","w");              /* timeseries sigma = 0.5 */
    fout1c = fopen("ts3.txt","w");              /* timeseries sigma = 1 */
    fout2a = fopen("phys_sigma0.txt","w");   /* physical node data for visualisations, sigma = 0, 1 realisation */
    fout2b = fopen("virt_sigma0.txt","w");   /* virtual node data, sigma = 0 */
    fout3a = fopen("phys_sigma1.txt","w");  /* physical node data, sigma = 1 */
    fout3b = fopen("virt_sigma1.txt","w");  /* virtual node data, sigma = 1 */
    
    /* initialize node array on physical layer */
    /* this section defines the lattice.  the lattice maps as a toroid so edges of square are mapped to opposite edges. */
    /* for each node, it's neighbouring node is identified via x and y integer coordinates. */
    /* 0 = west, 1 = south, 2 = east, 3 = north */
    
    for (i = 0 ; i < N ; i ++)
        for (j = 0 ; j < N ; j ++)
        {
            if (i == 0 && j == 0)
            {
                node[i][j].neighbour_x[0] = N-1;
                node[i][j].neighbour_y[0] = 0;
                node[i][j].neighbour_x[1] = 0;
                node[i][j].neighbour_y[1] = N-1;
                node[i][j].neighbour_x[2] = 1;
                node[i][j].neighbour_y[2] = 0;
                node[i][j].neighbour_x[3] = 0;
                node[i][j].neighbour_y[3] = 1;
            }
            else if (i == 0 && j == N-1)
            {
                node[i][j].neighbour_x[0] = N-1;
                node[i][j].neighbour_y[0] = N-1;
                node[i][j].neighbour_x[1] = 0;
                node[i][j].neighbour_y[1] = N-2;
                node[i][j].neighbour_x[2] = 1;
                node[i][j].neighbour_y[2] = N-1;
                node[i][j].neighbour_x[3] = 0;
                node[i][j].neighbour_y[3] = 0;
            }
            else if (i == N-1 && j == N-1)
            {
                node[i][j].neighbour_x[0] = N-2;
                node[i][j].neighbour_y[0] = N-1;
                node[i][j].neighbour_x[1] = N-1;
                node[i][j].neighbour_y[1] = N-2;
                node[i][j].neighbour_x[2] = 0;
                node[i][j].neighbour_y[2] = N-1;
                node[i][j].neighbour_x[3] = N-1;
                node[i][j].neighbour_y[3] = 0;
            }
            else if (i == N-1 && j == 0)
            {
                node[i][j].neighbour_x[0] = N-2;
                node[i][j].neighbour_y[0] = 0;
                node[i][j].neighbour_x[1] = N-1;
                node[i][j].neighbour_y[1] = N-1;
                node[i][j].neighbour_x[2] = 0;
                node[i][j].neighbour_y[2] = 0;
                node[i][j].neighbour_x[3] = N-1;
                node[i][j].neighbour_y[3] = 1;
            }
            else
            {
                if (i == 0)
                {
                    node[i][j].neighbour_x[0] = N-1;
                    node[i][j].neighbour_y[0] = j;
                    node[i][j].neighbour_x[1] = i;
                    node[i][j].neighbour_y[1] = j-1;
                    node[i][j].neighbour_x[2] = i+1;
                    node[i][j].neighbour_y[2] = j;
                    node[i][j].neighbour_x[3] = i;
                    node[i][j].neighbour_y[3] = j+1;
                }
                else if (j == 0)
                {
                    node[i][j].neighbour_x[0] = i-1;
                    node[i][j].neighbour_y[0] = j;
                    node[i][j].neighbour_x[1] = i;
                    node[i][j].neighbour_y[1] = N-1;
                    node[i][j].neighbour_x[2] = i+1;
                    node[i][j].neighbour_y[2] = j;
                    node[i][j].neighbour_x[3] = i;
                    node[i][j].neighbour_y[3] = j+1;
                }
                else if (i == N-1)
                {
                    node[i][j].neighbour_x[0] = i-1;
                    node[i][j].neighbour_y[0] = j;
                    node[i][j].neighbour_x[1] = i;
                    node[i][j].neighbour_y[1] = j-1;
                    node[i][j].neighbour_x[2] = 0;
                    node[i][j].neighbour_y[2] = j;
                    node[i][j].neighbour_x[3] = i;
                    node[i][j].neighbour_y[3] = j+1;
                }
                else if (j == N-1)
                {
                    node[i][j].neighbour_x[0] = i;
                    node[i][j].neighbour_y[0] = N-1;
                    node[i][j].neighbour_x[1] = i;
                    node[i][j].neighbour_y[1] = j-1;
                    node[i][j].neighbour_x[2] = i+1;
                    node[i][j].neighbour_y[2] = j;
                    node[i][j].neighbour_x[3] = i;
                    node[i][j].neighbour_y[3] = j+1;
                }
                else
                {
                    node[i][j].neighbour_x[0] = i-1;
                    node[i][j].neighbour_y[0] = j;
                    node[i][j].neighbour_x[1] = i;
                    node[i][j].neighbour_y[1] = j-1;
                    node[i][j].neighbour_x[2] = i+1;
                    node[i][j].neighbour_y[2] = j;
                    node[i][j].neighbour_x[3] = i;
                    node[i][j].neighbour_y[3] = j+1;
                }
            }
        }
    
    
    /* Initialize node array for the virual layer */
    
    /* Initially set 'virtual neighbours' to being the same as the 'physical neighbours' */
    
    for (i=0; i<N; i++){
        for(j=0; j<N; j++){
            for(m=0; m<4; m++){
                node[i][j].neighbour_virtual_x[m]=node[i][j].neighbour_x[m];
                node[i][j].neighbour_virtual_y[m]=node[i][j].neighbour_y[m];
            }
        }
    }
    
    /* Select sigma percent of nodes to exchange their 'virtual neighbours' for four other randomly chosen nodes */
    
    
    for(r = 0; r<3; r++){
        sigma = (float)r/2;
        
        printf("Sigma = %f\n",sigma);
        
        
        for (i=0; i < (int)(N * N * sigma); i++){
            
            /* Select a random node */
            
            j = (int)(kyrand() * N);
            k = (int)(kyrand() * N);
            
            
            /* Change its virtual neighbours to those of 4 randomly selected nodes */
            
            for(m=0; m<4; m++){
                node[j][k].neighbour_virtual_x[m] = (int)(kyrand()*N);
                node[j][k].neighbour_virtual_y[m] = (int)(kyrand()*N);
            }
            
            
            
        }
        
        /* initialize data arrays */
        for (i = 0 ; i < MAX_TS ; i ++)
        {
            S_timeseries[i] = 0;
            I_timeseries[i] = 0;
            R_timeseries[i] = 0;
            V_timeseries[i] = 0;
            pos_timeseries[i] = 0;
            neg_timeseries[i] = 0;
            
        }
        
        /* this is the loop for each simulation run */
        for (count_real = 0 ; count_real < NUM_REALIZATIONS ; count_real ++)
        {
            
            printf("Realization number %d of %d \n",count_real,NUM_REALIZATIONS);
            ts_count = 0;
            
            /* set each node to be susceptible and neutral at first */
            for (i = 0 ; i < N ; i ++)
            {
                for (j = 0 ; j < N ; j ++)
                {
                    node[i][j].status[0] = 's';
                    node[i][j].status[1] = 'n';
                }
            }
        
            
        
            /* INITIAL INFECTANTS CHOSEN AT RANDOM
            
            for (i = 0 ; i < INIT_inf ; i ++)
            {
                do{
                    j = (int)(N*kyrand());
                    k = (int)(N*kyrand());}
                while (node[j][k].status[0] != 's');
                
                node[j][k].status[0] = 'i';
            
            } */
            
            
            /*  INITIAL INFECTANTS CLUSTERED */
            
            for( i=N/2-2; i<N/2+2; i++){
                for( j=N/2-2; j<N/2+2; j++){
                    node[i][j].status[0]='i';
                }
            }
            
            
            
            /* Those initailly with positive vaccination views chosen at random */
            
            for(i = 0; i < INIT_pos; i++)
            {
                do{
                    j = (int)(N*kyrand());
                    k = (int)(N*kyrand());}
                while (node[j][k].status[1] != 'n');
                
                node[j][k].status[1] = '+';
            }
            
            
            
            /* Those initailly with negative vaccination views chosen at random */
            
            for(i = 0; i < INIT_neg; i++)
            {
                do{
                    j = (int)(N*kyrand());
                    k = (int)(N*kyrand());}
                while (node[j][k].status[1] != 'n');
                
                node[j][k].status[1] = '-';
            }
            
            
            
            /* let the disease spread for MAX_ITER timesteps and record data at intervals (every SAMPLE_INCR timesteps) */
            
            for (i = 0 ; i < MAX_ITER ; i ++)
            {
                /* select a random node and see if an event occurs */
                
                j = (int)(N*kyrand());
                k = (int)(N*kyrand());
                
                
                /* Count the number of positive / negative / infected individuals in neighbourhood */
                
                inf_count = 0;
                pos_count = 0;
                neg_count = 0;
                
                for (m = 0 ; m < 4 ; m ++)
                {
                    n_x = node[j][k].neighbour_x[m];
                    n_y = node[j][k].neighbour_y[m];
                    v_x = node[j][k].neighbour_virtual_x[m];
                    v_y = node[j][k].neighbour_virtual_y[m];
                    
                    if (node[n_x][n_y].status[0] == 'i') inf_count ++;
                    if (node[v_x][v_y].status[1] == '+') pos_count ++;
                    if (node[v_x][v_y].status[1] == '-') neg_count ++;
                }
                
                
                
                /* Susceptible-Positive individual outcome */
                
                if (node[j][k].status[0] == 's' && node[j][k].status[1] == '+')
                {
                    node[j][k].status[0] = 'v';
                }
                
                /* Susceptible-Neutral individual outcome */
                
                if (node[j][k].status[0] == 's' && node[j][k].status[1] == 'n')
                    
                {
                    if (kyrand() < inf_count * beta){
                        node[j][k].status[0] = 'i';
                    }
                    
                    if (kyrand() < max(beta_plus * pos_count - beta_minus * neg_count, 0)){
                        node[j][k].status[1] = '+';
                    }
                    
                    if (kyrand() < max(beta_minus * neg_count - beta_plus * pos_count, 0)){
                        node[j][k].status[1] = '-';
                    }
                    
                }
                    
                /* Susceptible-Negative individual outcome */
                
                if (node[j][k].status[0] == 's' && node[j][k].status[1] == '-')
                    
                {
                    if (kyrand() < inf_count * beta){
                        node[j][k].status[0] = 'i';
                    }
                    
                    if (kyrand() < max(beta_plus * pos_count - beta_minus * neg_count, 0)){
                        node[j][k].status[1] = 'n';
                    }
                        
                }
                
                /* Infected-Positive individual outcome */
        
                if (node[j][k].status[0] == 'i' && node[j][k].status[1] == '+')
                    
                {
                    if (kyrand() < mu){
                        node[j][k].status[0] = 'r';
                    }
                    
                    if (kyrand() < max(beta_minus * neg_count - beta_plus * pos_count, 0)){
                        node[j][k].status[1] = 'n';
                        
                    }
                
                }
                
                /* Infected-Neutral individual outcome */
                
                if (node[j][k].status[0] == 'i' && node[j][k].status[1] == 'n')
                    
                {
                    if (kyrand() < mu){
                        node[j][k].status[0] = 'r';
                    }
                    
                    if (kyrand() < max(beta_plus * pos_count - beta_minus * neg_count, 0)){
                        node[j][k].status[1] = '+';
                    }
                    
                    if (kyrand() < max(beta_minus * neg_count - beta_plus * pos_count, 0)){
                        node[j][k].status[1] = '-';
                        
                    }
                }
                
                /* Infected-Negative individual outcome */
                
                if (node[j][k].status[0] == 'i' && node[j][k].status[1] == '-')
                    
                {
                    if (kyrand() < mu){
                        node[j][k].status[0] = 'r';
                    }
                    
                    if (kyrand() < max(beta_plus * pos_count - beta_minus * neg_count, 0)){
                        node[j][k].status[1] = 'n';
                    }
                }
                
                /* Recovered-Positive individual outcome */
                
                if (node[j][k].status[0] == 'r' && node[j][k].status[1] == '+')
                    
                {
                    if (kyrand() < max(beta_minus * neg_count - beta_plus * pos_count, 0)){
                        node[j][k].status[1] = 'n';
                    }
                }
                
                
                /* Recovered-Neutral individual outcome */
                
                if (node[j][k].status[0] == 'r' && node[j][k].status[1] == 'n')
                    
                {
                    if (kyrand() < max(beta_plus * pos_count - beta_minus * neg_count, 0)){
                        node[j][k].status[1] = '+';
                    }
                    
                    if (kyrand() < max(beta_minus * neg_count - beta_plus * pos_count, 0)){
                        node[j][k].status[1] = '-';
                    }
                }
                
                /* Recovered-Negative individual outcome */
                
                if (node[j][k].status[0] == 'r' && node[j][k].status[1] == '-')
                    
                {
                    if (kyrand() < max(beta_plus * pos_count - beta_minus * neg_count, 0)){
                        node[j][k].status[1] = 'n';
                    }
                        
                }
                        
                /* Vaccinated-Positive individual stays as they are */
                
                
                /* every SAMPLE_INCR timesteps, update statistics and record data */
                if (i%SAMPLE_INCR == 0)
                {
                    
                    
                    
                    /* calculate totals of each state */
                    
                    sus_count = 0;
                    inf_count = 0;
                    rec_count = 0;
                    vac_count = 0;
                    pos_count = 0;
                    neg_count = 0;
                    
                    for (j = 0 ; j < N ; j ++){
                        for (k = 0 ; k < N ; k ++)
                        {
                            if (node[j][k].status[0] == 's') sus_count ++;
                            if (node[j][k].status[0] == 'i') inf_count ++;
                            if (node[j][k].status[0] == 'r') rec_count ++;
                            if (node[j][k].status[0] == 'v') vac_count ++;
                            if (node[j][k].status[1] == '+') pos_count ++;
                            if (node[j][k].status[1] == '-') neg_count ++;


                        }
                    }
                    
                    /* Print the node data at each sample increment time step for sigma = 0 */
                    /* This data is used for a visualisation - the second realisation is used */
                    if( count_real==1 && sigma==0 ){
                        
                        for(p=0; p<N; p++){
                            for(q=0; q<N; q++){
                                
                                fprintf(fout2a,"%c ",node[p][q].status[0]);
                                fprintf(fout2b,"%c ",node[p][q].status[1]);
                                
                            }
                            fprintf(fout2a,"\n");
                            fprintf(fout2b,"\n");
                        }
                        
                    }
                    
                    /* Print node data in the case sigma = 1 */
                    if( count_real==1 && sigma==1 ){
                        
                        for(p=0; p<N; p++){
                            for(q=0; q<N; q++){
                                
                                fprintf(fout3a,"%c ",node[p][q].status[0]);
                                fprintf(fout3b,"%c ",node[p][q].status[1]);
                                
                            }
                            fprintf(fout3a,"\n");
                            fprintf(fout3b,"\n");
                        }
                        
                    }
                    
                    /* Update time series */
                    
                    S_timeseries[ts_count] += sus_count;
                    I_timeseries[ts_count] += inf_count;
                    R_timeseries[ts_count] += rec_count;
                    V_timeseries[ts_count] += vac_count;
                    pos_timeseries[ts_count] += pos_count;
                    neg_timeseries[ts_count] += neg_count;
       
                    time_array[ts_count] = (float)i/N/N;  /* array that outputs the time */
                    ts_count ++; /* timestep counter */
                    
                    
                }
            }
        }
        
        /* determine average quantities across all realizations and print results to data files */
        
        
        for (i = 0 ; i < ts_count ; i ++)
        {
            S_timeseries[i] = S_timeseries[i]/NUM_REALIZATIONS;
            I_timeseries[i] = I_timeseries[i]/NUM_REALIZATIONS;
            R_timeseries[i] = R_timeseries[i]/NUM_REALIZATIONS;
            V_timeseries[i] = V_timeseries[i]/NUM_REALIZATIONS;
            pos_timeseries[i] = pos_timeseries[i]/NUM_REALIZATIONS;
            neg_timeseries[i] = neg_timeseries[i]/NUM_REALIZATIONS;
            
            if(r==0){
                fprintf(fout1a,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",time_array[i],S_timeseries[i],I_timeseries[i],R_timeseries[i],V_timeseries[i],pos_timeseries[i],neg_timeseries[i]);
                
            }
            if(r==1) {
                fprintf(fout1b,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",time_array[i],S_timeseries[i],I_timeseries[i],R_timeseries[i],V_timeseries[i],pos_timeseries[i],neg_timeseries[i]);
            }
            if(r==2){
                fprintf(fout1c,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",time_array[i],S_timeseries[i],I_timeseries[i],R_timeseries[i],V_timeseries[i],pos_timeseries[i],neg_timeseries[i]);
                
            }
            
        }
        
        
        
    }
    
 
    
    
    /* close data files */
    
    fclose(fout1a);
    fclose(fout1b);
    fclose(fout1c);
    fclose(fout2a);
    fclose(fout2b);
    fclose(fout3a);
    fclose(fout3b);
    
    
}


