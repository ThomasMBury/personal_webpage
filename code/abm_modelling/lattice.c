/*******************************************************************************
 
 Simulates an SIS (Susceptible-Infected-Susceptible) epidemic
 on a square lattice and outputs results. 4 neighbours per person.
 
 15 June 2006
 
 ********************************************************************************/

#define MAX_TS 1000
/* max number of timepoints in the output timeseries */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double kyrand()   /* this function returns a random number between 0 and 1; note that what random() returns may vary across computing platforms */
{
    int maxnumber = 2147483647l;   /* 2^31 - 1 */
    double answer;
    answer=((double)random())/maxnumber;
    return(answer);
}


int main()
{
    
    
    int i,j,k,m,ii;				/* assorted counters */
    int n_x,n_y; 					/* temporary counter variables */
    int n_inf,n_sus; 				/* temporary counter variables */
    int inf_count;				/* inf_count is the number of infected individuals */
    int phi_I_0,phi_I_1,phi_I_2,phi_I_3,phi_I_4; 	/* phi_I_n is the number of clusters of an infected surrounded by n infected individuals. */
    int phi_S_0,phi_S_1,phi_S_2,phi_S_3,phi_S_4;  /* phi_S_n is the number of clusters of a susceptible surrounded by n infected individuals. */
    int SI_count, II_count; 			/* number of SI, II edges */
    int count_real;			/* keeps count of the number of realizations */
    float CSI;					/* SI correlation function */
    
    int avg_count=0,avgI=0;			/* averages across all realizations */
    int avgSI_count=0,avgII_count=0;
    int avgphiI0=0,avgphiI1=0,avgphiI2=0,avgphiI3=0,avgphiI4=0;
    int avgphiS0=0,avgphiS1=0,avgphiS2=0,avgphiS3=0,avgphiS4=0;
    
    /* this is a structure which holds information about the status of each node and the identity of its neighbours */
    /* a node is specified by it's x/y integer coordinates */
    typedef struct node_data
    {
        int neighbour_x[4];                 /* identity of neighbours */
        int neighbour_y[4]; 		/* identity of neighbours */
        char status;                        /* infection status */
    } node_structure;                     /* structure for specifications of individual */
    
    /* file declarations */
    FILE *fout1;
    FILE *fout2;
    FILE *fout3;
    FILE *fout4;
    FILE *fout5;
    FILE *favg;
    FILE *fCSI;
    
    float STARTAVG = 0; 		/* number between 0 and 1 indicating when to start taking a running time-average */
    int N = 50;   		/* length of one side of the square lattice. population size = N times N */
    int SAMPLE_INCR =  1500;  	/* frequency with which lattice state is sampled for output files (greater value=less frequent sampling) */
    int MAX_ITER = 600000;       	/* number of iterations */
    int NUM_REALIZATIONS = 10;    /* number of simulation runs */
    float CSI_min = 10;		/* minimal CSI value (starting value set to 10 arbitrarily) */
    float lambda = 0.2/5; 	/* transmission rate along a lattice edge */
    float nu = 0.1;		/* recovery rate */
    int INOC = (int)(N*N/100); 	/* number of nodes initially infected */
    node_structure node[N][N]; 	/* data structure for the nodes */
    float CSI_timeseries[MAX_TS]; /* stores the timeseries of SI correlation values */
    float time_array[MAX_TS]; 	/* stores the timeseries of time points */
    float I_timeseries[MAX_TS]; 	/* stores the timeseries of number infected */
    int ts_count; 		/* keeps track of number of timesteps accumulated */
    
    /* output files */
    fout1 = fopen("ts.dat","w");     	/* number infected over time, average across all simulations */
    fCSI = fopen("CSI.dat","w");		/* SI correlation values over time, averaged across all simualation */
    fout2 = fopen("phiI.dat","w"); 	/* numbers of infected clusters over time, for each simulation */
    fout3 = fopen("phiS.dat","w"); 	/* numbers of susceptible clusters over time, for each simulation */
    fout4 = fopen("QSI.dat","w"); 	/* average number of susceptible neighbours per infected, for each simulation */
    fout5 = fopen("QII.dat","w"); 	/* average number of infected neighbours per infected, for each simulation */
    favg = fopen("avg.dat","w"); 		/* various statistics averaged across both time and simulation runs */
    
    /* initialize node array */
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
    
    /* initialize data arrays */
    for (i = 0 ; i < MAX_TS ; i ++)
    {
        CSI_timeseries[i] = 0;
        I_timeseries[i] = 0;
    }
    
    /* this is the loop for each simulation run */
    for (count_real = 0 ; count_real < NUM_REALIZATIONS ; count_real ++)
    {
        
        printf("Realization number %d of %d \n",count_real,NUM_REALIZATIONS);
        ts_count = 0;
        
        /* set each node to be susceptible at first */
        for (i = 0 ; i < N ; i ++)
            for (j = 0 ; j < N ; j ++)
            {
                node[i][j].status = 's';
            }
        
        /* and then choose some of the nodes to inoculate the disease in */
        for (i = 0 ; i < INOC ; i ++)
        {
            do{
                j = (int)(N*kyrand());
                k = (int)(N*kyrand());}
            while (node[j][k].status != 's');
            
            node[j][k].status = 'i';
        }
        
        /* let the disease spread for MAX_ITER timesteps and record data at intervals (every SAMPLE_INCR timesteps) */
        
        for (i = 0 ; i < MAX_ITER ; i ++)
        {
            /* select a random node and see if an event occurs */
            
            j = (int)(N*kyrand());
            k = (int)(N*kyrand());
            
            if (node[j][k].status == 's') /* if susceptible, see if they get infected by a neighbour and update accordingly */
            {
                inf_count = 0;
                for (m = 0 ; m < 4 ; m ++)
                {
                    n_x = node[j][k].neighbour_x[m];
                    n_y = node[j][k].neighbour_y[m];
                    if (node[n_x][n_y].status == 'i') inf_count ++;
                }
                if (kyrand() < inf_count*lambda)
                {
                    node[j][k].status = 'i';
                }
            }
            else /* if infected, see if they recover and update accordingly */
            {
                if (kyrand() < nu)
                {
                    node[j][k].status = 's';
                }
            }
            
            /* every SAMPLE_INCR timesteps, update statistics and record data */
            if (i%SAMPLE_INCR == 0)
            {
                /* statistics calculated for number of infecteds, susceptibles, various cluster types */
                
                inf_count = 0;
                phi_I_0 = phi_I_1 = phi_I_2 = phi_I_3 = phi_I_4 = phi_S_0 = phi_S_1 = phi_S_2 = phi_S_3 = phi_S_4 = 0;
                SI_count = II_count = 0;
                
                for (j = 0 ; j < N ; j ++)
                    for (k = 0 ; k < N ; k ++)
                    {
                        if (node[j][k].status == 'i') inf_count ++;
                        n_inf = 0;
                        for (m = 0 ; m < 4 ; m ++)
                        {
                            n_x = node[j][k].neighbour_x[m];
                            n_y = node[j][k].neighbour_y[m];
                            if (node[n_x][n_y].status == 'i') n_inf++;
                        }
                        if (node[j][k].status == 'i')
                        {
                            if (n_inf == 0) phi_I_0 ++;
                            else if (n_inf == 1)
                            {
                                phi_I_1 ++;
                                II_count ++;
                            }
                            else if (n_inf == 2)
                            {
                                phi_I_2 ++;
                                II_count += 2;
                            }
                            else if (n_inf == 3)
                            {
                                phi_I_3 ++;
                                II_count += 3;
                            }
                            else
                            {
                                phi_I_4 ++;
                                II_count += 4;
                            }
                        }
                        else
                        {
                            if (n_inf == 0) phi_S_0 ++;
                            else if (n_inf == 1)
                            {
                                phi_S_1 ++;
                                SI_count += 1;
                            }
                            else if (n_inf == 2)
                            {
                                phi_S_2 ++;
                                SI_count += 2;
                            }
                            else if (n_inf == 3)
                            {
                                phi_S_3 ++;
                                SI_count += 3;
                            } 
                            else 
                            { 
                                phi_S_4 ++; 
                                SI_count += 4; 
                            } 
                        } 
                    } 
                if (i > (int)(STARTAVG*MAX_ITER)) /* update averages */ 
                { 
                    avg_count ++; 
                    avgI += inf_count; 
                    avgphiI0 += phi_I_0;
                    avgphiI1 += phi_I_1; 
                    avgphiI2 += phi_I_2; 
                    avgphiI3 += phi_I_3; 
                    avgphiI4 += phi_I_4; 
                    avgphiS0 += phi_S_0; 
                    avgphiS1 += phi_S_1; 
                    avgphiS2 += phi_S_2; 
                    avgphiS3 += phi_S_3; 
                    avgphiS4 += phi_S_4; 
                    avgSI_count += SI_count; 
                    avgII_count += II_count; 
                } 
                
                /* print out values for this timestep */ 
                fprintf(fout2,"%d %f %d %d %d %d %d %f %f %f %f %f \n",i,(float)i/N/N,phi_I_0,phi_I_1,phi_I_2,phi_I_3,phi_I_4,(float)phi_I_0/N/N,(float)phi_I_1/N/N,(float)phi_I_2/N/N,(float)phi_I_3/N/N,(float)phi_I_4/N/N);
                fprintf(fout3,"%d %f %d %d %d %d %d %f %f %f %f %f \n",i,(float)i/N/N,phi_S_0,phi_S_1,phi_S_2,phi_S_3,phi_S_4,(float)phi_S_0/N/N,(float)phi_S_1/N/N,(float)phi_S_2/N/N,(float)phi_S_3/N/N,(float)phi_S_4/N/N);
                fprintf(fout4,"%d %f %d %d %f %f \n",i,(float)i/N/N,SI_count,N*N - inf_count,(float)SI_count/inf_count,(float)SI_count/(4*inf_count*(1-(float)inf_count/N/N)));
                fprintf(fout5,"%d %f %d %d %f \n",i,(float)i/N/N,II_count,inf_count,(float)II_count/inf_count); 
                
                if (inf_count != 0)
                    CSI_timeseries[ts_count] += (float)SI_count/(4*inf_count*(1-(float)inf_count/N/N)); 
                
                I_timeseries[ts_count] += inf_count;
                
                time_array[ts_count] = (float)i/N/N;  /* array that outputs the time */ 
                ts_count ++; /* timestep counter */ 
            } 
        } 
    }
    
    
    /* determine average quantities across all realizations and print results to data files */ 
    
    for (i = 0 ; i < ts_count ; i ++)
    { 
        CSI_timeseries[i] = CSI_timeseries[i]/NUM_REALIZATIONS; 
        I_timeseries[i] = I_timeseries[i]/NUM_REALIZATIONS;
        
        fprintf(fCSI,"%f %f \n",time_array[i],CSI_timeseries[i]); 
        fprintf(fout1,"%f %f \n",time_array[i],I_timeseries[i]);
        
        if (CSI_timeseries[i] < CSI_min) 
        { 
            CSI_min = CSI_timeseries[i]; 
        } 
    } 
    
    fprintf(fCSI,"\n%f \n",CSI_min); 
    
    /* print out values averaged across all simulations */ 
    
    fprintf(favg,"averaged values after a proportion %f of the time series is complete: \n",STARTAVG); 
    fprintf(favg,"samples = %d, I = %f, I(prop) = %f, SI = %f, II = %f \n",avg_count,(float)avgI/avg_count,(float)avgI/avg_count/N/N,(float)avgSI_count/avg_count,(float)avgII_count/avg_count); 
    fprintf(favg,"phi_I_0 = %f, phi_I_1 =  %f, phi_I_2 =  %f, phi_I_3 =  %f, phi_I_4 =  %f \n",(float)avgphiI0/avg_count/N/N,(float)avgphiI1/avg_count/N/N,(float)avgphiI2/avg_count/N/N,(float)avgphiI3/avg_count/N/N,(float)avgphiI4/avg_count/N/N); 
    fprintf(favg,"phi_S_0 = %f, phi_S_1 =  %f, phi_S_2 =  %f, phi_S_3 =  %f, phi_S_4 =  %f \n",(float)avgphiS0/avg_count/N/N,(float)avgphiS1/avg_count/N/N,(float)avgphiS2/avg_count/N/N,(float)avgphiS3/avg_count/N/N,(float)avgphiS4/avg_count/N/N);
    
    /* close data files */ 
    
    fclose(fout1);
    fclose(fout2);
    fclose(fout3);
    fclose(fout4);
    fclose(fout5);
    fclose(favg);
    fclose(fCSI);
}
