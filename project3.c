#include <stdio.h>
#include <time.h>

#include "structs.h"
#include "settings.h"
#include "functions.c"

int main(int argc, char ** argv)
{
    // Verifying Inputs
        // if(VerInput(argc,argv) == 0) return -1;
    
        IntegerPoints * IP   = calloc(1, sizeof(IntegerPoints *));
        HalfPoints    * HP   = calloc(1, sizeof(HalfPoints    *));
        MeshGrid      * mesh = calloc(1, sizeof(MeshGrid      *));
        ThomStr       * T    = calloc(1, sizeof(ThomStr       *));
        
        IP->M_INF=0.0;
        
        sscanf(argv[2], "%f", &IP->M_INF);
        
    // Control Variables
        int n=0;                                    // iteration number
        clock_t start_it, end_it;                   // time per iteration
        clock_t start_loop, end_loop;               // loop time
        double cpu_time_used, cpu_time_loop;
        
    // Allocating Memory
        AllocateAll(mesh,HP,IP,T);
        
    // Reading Mesh
        ReadMesh(mesh,argv[1]);
    
    // Initial Conditions
        InitialCond(mesh, IP);
        
    // Constant Initial Settings and Calculations
        ConstSetCal(mesh,HP,IP);
        
    /* Solver - Holst & Ballhaus method of Artificial Density 
                    for Full-Potential Equation Solution              */
        
        start_loop = clock(); // Starting loop clock
        while(n<NMAX)
        {
            start_it = clock(); // STARTING CLOCK! TIC TAC...TIC TAC
            
            // Variable Settings and Calculations
                VarSetCal(HP,IP);
                
            // Calculate Residue
                Residue(HP,IP);
                    
            // Convergence Test
                if(IP->resMax<RMIN)
                {
                    PrintConv(n,IP);
                    break;
                }
                else if(IP->resMax>RMAX)
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
                
            // Boundary Conditions
                BCs(IP,mesh);
                
            n++;
        }
        end_loop = clock();
        cpu_time_loop = ((double) (end_loop-start_loop)) / CLOCKS_PER_SEC;
        
    // Output
        Output(mesh,IP);
        
    return 0;
}