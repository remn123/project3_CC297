#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "settings.h"
#include "functions.h"

/*
    ./project3 mesh.dat 0.70
*/

int main(int argc, char ** argv)
{
    // Verifying Inputs
        // if(VerInput(argc,argv) == 0) return -1;
    
        IntegerPoints * IP   = calloc(1, sizeof(IntegerPoints ));
        HalfPoints    * HP   = calloc(1, sizeof(HalfPoints    ));
        MeshGrid      * mesh = calloc(1, sizeof(MeshGrid      ));
        ThomStr       * T    = calloc(1, sizeof(ThomStr       ));
        
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

        // Boundary Conditions
        BCs(IP,HP,mesh);
        printf("Boundary Conditions Applied\n");

    /* Solver - Holst & Ballhaus method of Artificial Density 
                    for Full-Potential Equation Solution              */
        printf("Starting iteration...\n");
        start_loop = clock(); // Starting loop clock
        while(n<NMAX)
        {
            start_it = clock(); // STARTING CLOCK! TIC TAC...TIC TAC
                printf("Iteration(%d):\n",n);

            // Variable Settings and Calculations
                VarSetCal(HP,IP);
            
             // Boundary Conditions
                BCs(IP,HP,mesh);
                printf("Boundary Conditions Applied\n");

            // Calculate Residue
                Residue(HP,IP,n);
                printf("Residue(%d) Calculated\n",n);     
            // Convergence Test
                if(IP->resL2[n]<RMIN)
                {
                    PrintConv(n,IP);
                    break;
                }
                else if(IP->resL2[n]>RMAX)
                {
                    PrintDiverg(n,IP);
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
                printf("Updating fi...\n");
            // Updating Potencial Flow
                for(i=1;i<IMAX;i++)
                {
                    for(j=1;j<JMAX;j++)
                    {
                        printf("IP->corr[%d][%d]=%lf\n",i,j,IP->corr[i][j]);
                        IP->fi[i][j] = IP->fi[i][j] + IP->corr[i][j];
                    }   
                }
                

            // Boundary Conditions
                BCs(IP,HP,mesh);
                printf("Boundary Conditions Applied\n");

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