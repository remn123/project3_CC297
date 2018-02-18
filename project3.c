#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <math.h>
#include "structs.h"
#include "settings.h"
#include "functions.h"

// Debugging
#define _GNU_SOURCE         /* See feature_test_macros(7) */
#include <fenv.h>
int feenableexcept(int excepts);

/*
    ./project3 mesh.dat 0.70
*/

int main(int argc, char ** argv)
{
    // Verifying Inputs
        // if(VerInput(argc,argv) == 0) return -1;
        feenableexcept(FE_INVALID | FE_OVERFLOW);

        IntegerPoints * IP   = calloc(1, sizeof(IntegerPoints ));
        HalfPoints    * HP   = calloc(1, sizeof(HalfPoints    ));
        MeshGrid      * mesh = calloc(1, sizeof(MeshGrid      ));
        ThomStr       * T    = calloc(1, sizeof(ThomStr       ));
        double fi_inf, U_INF;
        
        IP->M_INF=0.0;
        
        sscanf(argv[2], "%lf", &IP->M_INF);
        
    // Control Variables
        int n=0,i=0,j=0;                            // iteration number
        clock_t start_it, end_it;                   // time per iteration
        clock_t start_loop, end_loop;               // loop time
        double cpu_time_used, cpu_time_loop;
            
    // Reading Mesh
        ReadMesh(mesh,argv[1],IP,HP,T);
        printf("Mesh has been read!\n");
    
    // Initial Conditions
        InitialCond(mesh, IP);
        printf("Initial Condition assigned!\n");    
    
    // Constant Initial Settings and Calculations
        ConstSetCal(mesh,HP,IP);
        printf("Initial Settings Calculated!\n");
        printf("\n");

        
    /* Solver - Holst & Ballhaus method of Artificial Density 
                    for Full-Potential Equation Solution              */
        printf("Starting iteration...\n");
        start_loop = clock(); // Starting loop clock
        while(n<NMAX)
        {
            start_it = clock(); // STARTING CLOCK! TIC TAC...TIC TAC
                // printf("Iteration(%d):\n",n);

            // Variable Settings and Calculations
                VarSetCal(HP,IP,mesh);
            
            //  // Boundary Conditions
            
            //     printf("Boundary Conditions Applied\n");

            // Calculate Residue
                Residue(HP,IP,n);
                // printf("Residue(%d) Calculated\n",n);     
            // Convergence Test
                if(IP->resL2[n]<RMIN)
                {
                    PrintConv(n,IP);
                    end_it = clock();
        	    	cpu_time_used = ((double) (end_it-start_it)) / CLOCKS_PER_SEC;
                    break;
                }
                else if(IP->resL2[n]>RMAX)
                {
                    PrintDiverg(n,IP);
                    end_it = clock();
        	    	cpu_time_used = ((double) (end_it-start_it)) / CLOCKS_PER_SEC;
                    break;
                }
                else
                {
                    end_it = clock();
        	    	cpu_time_used = ((double) (end_it-start_it)) / CLOCKS_PER_SEC;
                    
                    PrintIter(n,IP,cpu_time_used);
                    
                }
                
            // Iteration Method
                IterMethod(IP,HP,T,n);
                // printf("Updating fi...\n");
            // Updating Potencial Flow
                for(i=1;i<IMAX;i++)
                {
                    for(j=1;j<JMAX;j++)
                    {
                        // printf("IP->corr[%d][%d]=%lf\n",i,j,IP->corr[i][j]);
                        IP->fi[i][j] = IP->fi[i][j] + IP->corr[i][j];
                    }   
                }
            
                U_INF = pow( (gamma+1.0)/(gamma-1.0+2.0/(IP->M_INF*IP->M_INF)),0.5 );
            
            // Boundary Conditions
               /** At the outer boundary */
            	for (int i = 1 ; i <= IMAX-1; i++ )
            	{
            		fi_inf = U_INF*mesh->x[i][0];
            		
            		IP->fi[i][0] = fi_inf; 
            
            	}
            	
                /** Enforce jump condition */
                for (int j = 0; j <= JMAX-1; j++)
                {
            		IP->fi[0][j] =  IP->fi[IMAX-1][j];
            // 		IP->U[0][j] =  IP->U[IMAX-1][j];
            	}
                // printf("Boundary Conditions Applied\n");

                Output(mesh,IP);
            n++;
        }
        end_loop = clock();
        cpu_time_loop = ((double) (end_loop-start_loop)) / CLOCKS_PER_SEC;
        
    // Output
        Output(mesh,IP);

    // Free up memmory
        FreeUp(IP,HP,mesh,T);
        
    return 0;
}