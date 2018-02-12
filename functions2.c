#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"
#include "settings.h"

/* This file contains all functions that were used on project 3 */
/* ------------------------------------------------------------ */

/* ----------------------- Verifying Input -------------------- */
// ------------------------------------------------------------ //
    // int VerInput( int argc, char ** argv)
    // {
        // bla bla bla
    // }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* --------------------- Allocating Memory -------------------- */
// ------------------------------------------------------------ //
    void AllocateAll(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP, ThomStr * T)
    {
        int i;
        
        // mesh
        mesh->x     = calloc(IMAX, sizeof(double *));
        mesh->y     = calloc(IMAX, sizeof(double *));
        // IP
        IP->resMax  = calloc(1   , sizeof(double ));
        IP->beta    = calloc(1   , sizeof(double ));
        IP->resL1   = calloc(1   , sizeof(double ));
        IP->resL2   = calloc(1   , sizeof(double ));
        IP->residue = calloc(IMAX, sizeof(double *));
        IP->fi      = calloc(IMAX, sizeof(double *));
        IP->corr    = calloc(IMAX, sizeof(double *));
        IP->rho     = calloc(IMAX, sizeof(double *));
        IP->x_csi   = calloc(IMAX, sizeof(double *));
        IP->x_eta   = calloc(IMAX, sizeof(double *));
        IP->y_csi   = calloc(IMAX, sizeof(double *));
        IP->y_eta   = calloc(IMAX, sizeof(double *));
        
        IP->csi_x   = calloc(IMAX, sizeof(double *));
        IP->csi_y   = calloc(IMAX, sizeof(double *));
        IP->eta_x   = calloc(IMAX, sizeof(double *));
        IP->eta_y   = calloc(IMAX, sizeof(double *));
        IP->J       = calloc(IMAX, sizeof(double *));
        IP->A1      = calloc(IMAX, sizeof(double *));
        IP->A2      = calloc(IMAX, sizeof(double *));
        IP->A3      = calloc(IMAX, sizeof(double *));
        
        IP->fi_csi  = calloc(IMAX, sizeof(double *));
        IP->fi_eta  = calloc(IMAX, sizeof(double *));
        IP->U       = calloc(IMAX, sizeof(double *));
        IP->V       = calloc(IMAX, sizeof(double *));
        
        // HP
        HP->xhi_csi    = calloc(IMAX, sizeof(double *));
        HP->yhi_csi    = calloc(IMAX, sizeof(double *));
        HP->xhj_eta    = calloc(IMAX, sizeof(double *));
        HP->yhj_eta    = calloc(IMAX, sizeof(double *));
          
        HP->xhj_csi    = calloc(IMAX, sizeof(double *));
        HP->yhj_csi    = calloc(IMAX, sizeof(double *));
        HP->xhi_eta    = calloc(IMAX, sizeof(double *));
        HP->yhi_eta    = calloc(IMAX, sizeof(double *));
          
        HP->csi_xhi    = calloc(IMAX, sizeof(double *));
        HP->csi_yhi    = calloc(IMAX, sizeof(double *));
        HP->eta_xhi    = calloc(IMAX, sizeof(double *));
        HP->eta_yhi    = calloc(IMAX, sizeof(double *));
        
        HP->csi_xhj    = calloc(IMAX, sizeof(double *));
        HP->csi_yhj    = calloc(IMAX, sizeof(double *));
        HP->eta_xhj    = calloc(IMAX, sizeof(double *));
        HP->eta_yhj    = calloc(IMAX, sizeof(double *));
        
        HP->A1_hi      = calloc(IMAX, sizeof(double *));
        HP->A2_hi      = calloc(IMAX, sizeof(double *));
        HP->A3_hi      = calloc(IMAX, sizeof(double *));
        
        HP->A1_hj      = calloc(IMAX, sizeof(double *));
        HP->A2_hj      = calloc(IMAX, sizeof(double *));
        HP->A3_hj      = calloc(IMAX, sizeof(double *));
    
        HP->J_hi       = calloc(IMAX, sizeof(double *));
        HP->J_hj       = calloc(IMAX, sizeof(double *));
        
        HP->fi_csi_hi  = calloc(IMAX, sizeof(double *));
        HP->fi_csi_hj  = calloc(IMAX, sizeof(double *));
        HP->fi_eta_hi  = calloc(IMAX, sizeof(double *));
        HP->fi_eta_hj  = calloc(IMAX, sizeof(double *));
        
        HP->U_hi       = calloc(IMAX, sizeof(double *));
        HP->U_hj       = calloc(IMAX, sizeof(double *));
        HP->V_hi       = calloc(IMAX, sizeof(double *));
        HP->V_hj       = calloc(IMAX, sizeof(double *));
        
        HP->rho_hi     = calloc(IMAX, sizeof(double *));
        HP->rho_hj     = calloc(IMAX, sizeof(double *));
        
        HP->rho_til_hi = calloc(IMAX, sizeof(double *));
        HP->rho_bar_hj = calloc(IMAX, sizeof(double *));
        
        HP->nu_hi      = calloc(IMAX, sizeof(double *));
        HP->nu_hj      = calloc(IMAX, sizeof(double *));
        
        // T
        T->a           = calloc(IMAX, sizeof(double  ));
        T->b           = calloc(IMAX, sizeof(double  ));
        T->c           = calloc(IMAX, sizeof(double  ));
        T->RHS         = calloc(IMAX, sizeof(double  ));
        T->LHS         = calloc(IMAX, sizeof(double  ));
          
        for(i = 0; i < IMAX; i++) { 
            mesh->x[i]        = calloc(JMAX, sizeof(double));
            mesh->y[i]        = calloc(JMAX, sizeof(double));
            IP->residue[i]    = calloc(JMAX, sizeof(double));
            IP->fi[i]         = calloc(JMAX, sizeof(double));
            IP->corr[i]       = calloc(JMAX, sizeof(double));
            IP->rho[i]        = calloc(JMAX, sizeof(double));
            IP->x_csi[i]      = calloc(JMAX, sizeof(double));
            IP->x_eta[i]      = calloc(JMAX, sizeof(double));
            IP->y_csi[i]      = calloc(JMAX, sizeof(double));
            IP->y_eta[i]      = calloc(JMAX, sizeof(double));
            
            IP->csi_x[i]      = calloc(JMAX, sizeof(double));
            IP->csi_y[i]      = calloc(JMAX, sizeof(double));
            IP->eta_x[i]      = calloc(JMAX, sizeof(double));
            IP->eta_y[i]      = calloc(JMAX, sizeof(double));
            IP->J[i]          = calloc(JMAX, sizeof(double));
            IP->A1[i]         = calloc(JMAX, sizeof(double));
            IP->A2[i]         = calloc(JMAX, sizeof(double));
            IP->A3[i]         = calloc(JMAX, sizeof(double));
              
            IP->fi_csi[i]     = calloc(JMAX, sizeof(double));
            IP->fi_eta[i]     = calloc(JMAX, sizeof(double));
            IP->U[i]          = calloc(JMAX, sizeof(double));
            IP->V[i]          = calloc(JMAX, sizeof(double));
            
            HP->xhi_csi[i]    = calloc(JMAX, sizeof(double));
            HP->yhi_csi[i]    = calloc(JMAX, sizeof(double));
            HP->xhj_eta[i]    = calloc(JMAX, sizeof(double));
            HP->yhj_eta[i]    = calloc(JMAX, sizeof(double));
                
            HP->xhj_csi[i]    = calloc(JMAX, sizeof(double));
            HP->yhj_csi[i]    = calloc(JMAX, sizeof(double));
            HP->xhi_eta[i]    = calloc(JMAX, sizeof(double));
            HP->yhi_eta[i]    = calloc(JMAX, sizeof(double));
              
            HP->csi_xhi[i]    = calloc(JMAX, sizeof(double));
            HP->csi_yhi[i]    = calloc(JMAX, sizeof(double));
            HP->eta_xhi[i]    = calloc(JMAX, sizeof(double));
            HP->eta_yhi[i]    = calloc(JMAX, sizeof(double));
             
            HP->csi_xhj[i]    = calloc(JMAX, sizeof(double));
            HP->csi_yhj[i]    = calloc(JMAX, sizeof(double));
            HP->eta_xhj[i]    = calloc(JMAX, sizeof(double));
            HP->eta_yhj[i]    = calloc(JMAX, sizeof(double));
            
            HP->A1_hi[i]      = calloc(JMAX, sizeof(double));
            HP->A2_hi[i]      = calloc(JMAX, sizeof(double));
            HP->A3_hi[i]      = calloc(JMAX, sizeof(double));
             
            HP->A1_hj[i]      = calloc(JMAX, sizeof(double));
            HP->A2_hj[i]      = calloc(JMAX, sizeof(double));
            HP->A3_hj[i]      = calloc(JMAX, sizeof(double));
        
            HP->J_hi[i]       = calloc(JMAX, sizeof(double));
            HP->J_hj[i]       = calloc(JMAX, sizeof(double));
            
            HP->fi_csi_hi[i]  = calloc(JMAX, sizeof(double));
            HP->fi_csi_hj[i]  = calloc(JMAX, sizeof(double));
            HP->fi_eta_hi[i]  = calloc(JMAX, sizeof(double));
            HP->fi_eta_hj[i]  = calloc(JMAX, sizeof(double));
            
            HP->U_hi[i]       = calloc(JMAX, sizeof(double));
            HP->U_hj[i]       = calloc(JMAX, sizeof(double));
            HP->V_hi[i]       = calloc(JMAX, sizeof(double));
            HP->V_hj[i]       = calloc(JMAX, sizeof(double));
            
            HP->rho_hi[i]     = calloc(JMAX, sizeof(double));
            HP->rho_hj[i]     = calloc(JMAX, sizeof(double));
            
            HP->rho_til_hi[i] = calloc(JMAX, sizeof(double));
            HP->rho_bar_hj[i] = calloc(JMAX, sizeof(double));
            
            HP->nu_hi[i]      = calloc(JMAX, sizeof(double));
            HP->nu_hj[i]      = calloc(JMAX, sizeof(double));

        }
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Initial Condition --------------------- */
// ------------------------------------------------------------ //
    void InitialCond(MeshGrid * mesh, IntegerPoints * IP)
    {
        int i=0,j=0;
        double fi_inf, U_INF;
    
        U_INF = pow((gamma+1.0)/(gamma-1.0+2.0/(IP->M_INF*IP->M_INF)),0.5);
    	printf("M_INF= %lf; U_INF= %lf\n",IP->M_INF,U_INF);
        for(i=0;i<IMAX;i++)
        {
        	// External Boundary Condition
	        fi_inf = U_INF*mesh->x[i][0];
	            
        	for(j=0;j<JMAX;j++)
        	{
	            IP->fi[i][j] = fi_inf;
	        }
        }   
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* --------------------- Read Mesh ---------------------------- */
// ------------------------------------------------------------ //
    void ReadMesh(MeshGrid * mesh, char * filename, IntegerPoints * IP, HalfPoints * HP, ThomStr * T)
    {
        int nelm=0, i=0, j=0;
        double x=0.0, y=0.0;
        
        FILE * fr;
        
        fr = fopen(filename,"r");
        fscanf(fr,"%d %d",&IMAX,&JMAX); // reading the mesh size on global variables
        printf("IMAX: %d \n JMAX: %d \n",IMAX,JMAX);
        
        // Allocating Memory
        AllocateAll(mesh,HP,IP,T);
        printf("Memmory Allocated!\n");

        while(fscanf(fr,"%lf %lf",&x,&y)!=EOF)
        {	
        	mesh->x[i][j] = x;
            mesh->y[i][j] = y;
            
            if(j==JMAX-1)
            {
                j=0;
                i++;
            }
            else
            {
                j++;
            }
            
            nelm++;
        }
        fclose(fr);
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* ------- Constant Initial Settings and Calculations --------- */
// ------------------------------------------------------------ //
    void InverseMesh(MeshGrid * mesh)
    {
        double auX=0.0, auY=0.0;
        int i=0, j=0;
        
        // Inverting csi direction
        for(i=0;i<(IMAX+1)/2;i++)
        {
            for(j=0;j<JMAX;j++)
            {
                auX = mesh->x[i][j];
                mesh->x[i][j] = mesh->x[IMAX-1-i][j];
                mesh->x[IMAX-1-i][j] = auX;
                
                auY = mesh->y[i][j];
                mesh->y[i][j] = mesh->y[IMAX-1-i][j];
                mesh->y[IMAX-1-i][j] = auY;
            }
        }
        
        // Inverting eta direction
        for(i=0;i<IMAX;i++)
        {
            for(j=0;j<(JMAX+1)/2;j++)
            {
                auX = mesh->x[i][j];
                mesh->x[i][j] = mesh->x[i][JMAX-1-j];
                mesh->x[i][JMAX-1-j] = auX;
                
                auY = mesh->y[i][j];
                mesh->y[i][j] = mesh->y[i][JMAX-1-j];
                mesh->y[i][JMAX-1-j] = auY;
            }
        }
    }
    
    void HPderivPM(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {	
    	if(i==0 && j>0 && j<JMAX-1)
    	{
    		IP->x_csi[i][j] = 0.5*(mesh->x[i+1][j] - mesh->x[IMAX-2][j]);  // x_csi(i,j) 
	        IP->y_csi[i][j] = 0.5*(mesh->y[i+1][j] - mesh->y[IMAX-2][j]);  // y_csi(i,j)
	        IP->x_eta[i][j] = 0.5*(mesh->x[i][j+1] - mesh->x[i][j-1]);  // x_eta(i,j)
	        IP->y_eta[i][j] = 0.5*(mesh->y[i][j+1] - mesh->y[i][j-1]);  // y_eta(i,j) 

	        // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[i+1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[i+1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)

	        HP->xhj_csi[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[IMAX-2][j+1]) + (mesh->x[i+1][j] - mesh->x[IMAX-2][j]));  // x_csi(i,j+1/2)
	        HP->yhj_csi[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[IMAX-2][j+1]) + (mesh->y[i+1][j] - mesh->y[IMAX-2][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[i+1][j-1]) + (mesh->x[i][j+1] - mesh->x[i][j-1]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[i+1][j-1]) + (mesh->y[i][j+1] - mesh->y[i][j-1]));  // y_eta(i+1/2,j)	
    	}
    	else if(i>0 && i<IMAX-1 && j==0)
    	{
    		IP->x_csi[i][j] = 0.5*(mesh->x[i+1][j] - mesh->x[i-1][j]);  // x_csi(i,j) 
	        IP->y_csi[i][j] = 0.5*(mesh->y[i+1][j] - mesh->y[i-1][j]);  // y_csi(i,j)
	       	IP->x_eta[i][j] = 0.5*(-3.0*mesh->x[i][j] + 4.0*mesh->x[i][j+1] - mesh->x[i][j+2]);  // x_eta(i,j)
	        IP->y_eta[i][j] = 0.5*(-3.0*mesh->y[i][j] + 4.0*mesh->y[i][j+1] - mesh->y[i][j+2]);  // y_eta(i,j) 

	        // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[i+1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[i+1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)

	        HP->xhj_csi[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[i-1][j+1]) + (mesh->x[i+1][j] - mesh->x[i-1][j]));  // x_csi(i,j+1/2)
	        HP->yhj_csi[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[i-1][j+1]) + (mesh->y[i+1][j] - mesh->y[i-1][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*(-3.0*(mesh->x[i+1][j]+mesh->x[i][j]) +4.0*(mesh->x[i+1][j+1]+mesh->x[i][j+1]) -(mesh->x[i+1][j+2] + mesh->x[i][j+2]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*(-3.0*(mesh->y[i+1][j]+mesh->y[i][j]) +4.0*(mesh->y[i+1][j+1]+mesh->y[i][j+1]) -(mesh->y[i+1][j+2] + mesh->y[i][j+2]));  // y_eta(i+1/2,j)	
    	}
    	else if(i==0 && j==0)
    	{
    		IP->x_csi[i][j] = 0.5*(mesh->x[i+1][j] - mesh->x[IMAX-2][j]);  // x_csi(i,j) 
	        IP->y_csi[i][j] = 0.5*(mesh->y[i+1][j] - mesh->y[IMAX-2][j]);  // y_csi(i,j)
	        IP->x_eta[i][j] = 0.5*(-3.0*mesh->x[i][j] + 4.0*mesh->x[i][j+1] - mesh->x[i][j+2]);  // x_eta(i,j)
	        IP->y_eta[i][j] = 0.5*(-3.0*mesh->y[i][j] + 4.0*mesh->y[i][j+1] - mesh->y[i][j+2]);  // y_eta(i,j) 

	        // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[i+1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[i+1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)

	        HP->xhj_csi[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[IMAX-2][j+1]) + (mesh->x[i+1][j] - mesh->x[IMAX-2][j]));  // x_csi(i,j+1/2)
	        HP->yhj_csi[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[IMAX-2][j+1]) + (mesh->y[i+1][j] - mesh->y[IMAX-2][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*(-3.0*(mesh->x[i+1][j]+mesh->x[i][j]) +4.0*(mesh->x[i+1][j+1]+mesh->x[i][j+1]) -(mesh->x[i+1][j+2] + mesh->x[i][j+2]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*(-3.0*(mesh->y[i+1][j]+mesh->y[i][j]) +4.0*(mesh->y[i+1][j+1]+mesh->y[i][j+1]) -(mesh->y[i+1][j+2] + mesh->y[i][j+2]));  // y_eta(i+1/2,j)		
    	}
    	else if(i==IMAX-1 && j <JMAX-1)
    	{
			IP->x_csi[i][j] = 0.5*(mesh->x[1][j] - mesh->x[i-1][j]);  // x_csi(i,j) 
	        IP->y_csi[i][j] = 0.5*(mesh->y[1][j] - mesh->y[i-1][j]);  // y_csi(i,j)
    
    		if(j==0)
    		{
    			IP->x_eta[i][j] = 0.5*(-3.0*mesh->x[i][j] + 4.0*mesh->x[i][j+1] - mesh->x[i][j+2]);  // x_eta(i,j)
		        IP->y_eta[i][j] = 0.5*(-3.0*mesh->y[i][j] + 4.0*mesh->y[i][j+1] - mesh->y[i][j+2]);  // y_eta(i,j) 

    		}
    		else
    		{
    			IP->x_eta[i][j] = 0.5*(mesh->x[i][j+1] - mesh->x[i][j-1]);  // x_eta(i,j)
		        IP->y_eta[i][j] = 0.5*(mesh->y[i][j+1] - mesh->y[i][j-1]);  // y_eta(i,j) 

    		}
            // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)

	        HP->xhj_csi[i][j] = 0.25*((mesh->x[1][j+1] - mesh->x[i-1][j+1]) + (mesh->x[1][j] - mesh->x[i-1][j]));  // x_csi(i,j+1/2)
	        HP->yhj_csi[i][j] = 0.25*((mesh->y[1][j+1] - mesh->y[i-1][j+1]) + (mesh->y[1][j] - mesh->y[i-1][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*((mesh->x[1][j+1] - mesh->x[1][j-1]) + (mesh->x[i][j+1] - mesh->x[i][j-1]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*((mesh->y[1][j+1] - mesh->y[1][j-1]) + (mesh->y[i][j+1] - mesh->y[i][j-1]));  // y_eta(i+1/2,j)	
    	}
    	else if(i<IMAX-1 && j==JMAX-1)
    	{
    		if(i==0)
    		{
    			IP->x_csi[i][j] = 0.5*(mesh->x[i+1][j] - mesh->x[IMAX-2][j]);  // x_csi(i,j) 
	        	IP->y_csi[i][j] = 0.5*(mesh->y[i+1][j] - mesh->y[IMAX-2][j]);  // y_csi(i,j)
    		}
    		else
    		{
    			IP->x_csi[i][j] = 0.5*(mesh->x[i+1][j] - mesh->x[i-1][j]);  // x_csi(i,j) 
	        	IP->y_csi[i][j] = 0.5*(mesh->y[i+1][j] - mesh->y[i-1][j]);  // y_csi(i,j)
    		}
	        
	        IP->x_eta[i][j] = 0.5*(3.0*mesh->x[i][j] - 4.0*mesh->x[i][j-1] + mesh->x[i][j-2]);  // x_eta(i,j)
	        IP->y_eta[i][j] = 0.5*(3.0*mesh->y[i][j] - 4.0*mesh->y[i][j-1] + mesh->y[i][j-2]);  // y_eta(i,j) 

	        // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[i+1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[i+1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        // HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        // HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)
	        
	        // HP->xhj_csi[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[i-1][j+1]) + (mesh->x[i+1][j] - mesh->x[i-1][j]));  // x_csi(i,j+1/2)
	        // HP->yhj_csi[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[i-1][j+1]) + (mesh->y[i+1][j] - mesh->y[i-1][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*(3.0*(mesh->x[i+1][j]+mesh->x[i][j]) -4.0*(mesh->x[i+1][j-1]+mesh->x[i][j-1]) +(mesh->x[i+1][j-2] + mesh->x[i][j-2]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*(3.0*(mesh->y[i+1][j]+mesh->y[i][j]) -4.0*(mesh->y[i+1][j-1]+mesh->y[i][j-1]) +(mesh->y[i+1][j-2] + mesh->y[i][j-2]));  // y_eta(i+1/2,j)		
    	}
    	else if(i==IMAX-1 && j==JMAX-1)
    	{
    		IP->x_csi[i][j] = 0.5*(mesh->x[1][j] - mesh->x[i-1][j]);  // x_csi(i,j) 
	        IP->y_csi[i][j] = 0.5*(mesh->y[1][j] - mesh->y[i-1][j]);  // y_csi(i,j)
	        IP->x_eta[i][j] = 0.5*(3.0*mesh->x[i][j] - 4.0*mesh->x[i][j-1] + mesh->x[i][j-2]);  // x_eta(i,j)
	        IP->y_eta[i][j] = 0.5*(3.0*mesh->y[i][j] - 4.0*mesh->y[i][j-1] + mesh->y[i][j-2]);  // y_eta(i,j) 

	        // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        // HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        // HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)

	        // HP->xhj_csi[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[i-1][j+1]) + (mesh->x[i+1][j] - mesh->x[i-1][j]));  // x_csi(i,j+1/2)
	        // HP->yhj_csi[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[i-1][j+1]) + (mesh->y[i+1][j] - mesh->y[i-1][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*(3.0*(mesh->x[1][j]+mesh->x[i][j]) -4.0*(mesh->x[1][j-1]+mesh->x[i][j-1]) +(mesh->x[1][j-2] + mesh->x[i][j-2]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*(3.0*(mesh->y[1][j]+mesh->y[i][j]) -4.0*(mesh->y[1][j-1]+mesh->y[i][j-1]) +(mesh->y[1][j-2] + mesh->y[i][j-2]));  // y_eta(i+1/2,j)	
    	}
    	else
    	{
    		IP->x_csi[i][j] = 0.5*(mesh->x[i+1][j] - mesh->x[i-1][j]);  // x_csi(i,j) 
	        IP->y_csi[i][j] = 0.5*(mesh->y[i+1][j] - mesh->y[i-1][j]);  // y_csi(i,j)
	        IP->x_eta[i][j] = 0.5*(mesh->x[i][j+1] - mesh->x[i][j-1]);  // x_eta(i,j)
	        IP->y_eta[i][j] = 0.5*(mesh->y[i][j+1] - mesh->y[i][j-1]);  // y_eta(i,j) 

	        // centered in half step is equivalent to a forward in integer step
	        HP->xhi_csi[i][j] = (mesh->x[i+1][j] - mesh->x[i][j]);  // x_csi(i+1/2,j)
	        HP->yhi_csi[i][j] = (mesh->y[i+1][j] - mesh->y[i][j]);  // y_csi(i+1/2,j)
	        HP->xhj_eta[i][j] = (mesh->x[i][j+1] - mesh->x[i][j]);  // x_eta(i,j+1/2)
	        HP->yhj_eta[i][j] = (mesh->y[i][j+1] - mesh->y[i][j]);  // y_eta(i,j+1/2)

	        HP->xhj_csi[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[i-1][j+1]) + (mesh->x[i+1][j] - mesh->x[i-1][j]));  // x_csi(i,j+1/2)
	        HP->yhj_csi[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[i-1][j+1]) + (mesh->y[i+1][j] - mesh->y[i-1][j]));  // y_csi(i,j+1/2)
	        HP->xhi_eta[i][j] = 0.25*((mesh->x[i+1][j+1] - mesh->x[i+1][j-1]) + (mesh->x[i][j+1] - mesh->x[i][j-1]));  // x_eta(i+1/2,j)
	        HP->yhi_eta[i][j] = 0.25*((mesh->y[i+1][j+1] - mesh->y[i+1][j-1]) + (mesh->y[i][j+1] - mesh->y[i][j-1]));  // y_eta(i+1/2,j)	
    	}

        
    }
    
    void HPjacobian(HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {
        IP->J[i][j]    = 1.0/(IP->x_csi[i][j]*IP->y_eta[i][j] - IP->x_eta[i][j]*IP->y_csi[i][j]);
        HP->J_hi[i][j] = 1.0/(HP->xhi_csi[i][j]*HP->yhi_eta[i][j] - HP->xhi_eta[i][j]*HP->yhi_csi[i][j]);
        HP->J_hj[i][j] = 1.0/(HP->xhj_csi[i][j]*HP->yhj_eta[i][j] - HP->xhj_eta[i][j]*HP->yhj_csi[i][j]);
    }
    
    void HPderivNM(HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {
        IP->csi_x[i][j] =  IP->J[i][j]*IP->y_eta[i][j];
        IP->csi_y[i][j] = -IP->J[i][j]*IP->x_eta[i][j];
        IP->eta_x[i][j] = -IP->J[i][j]*IP->y_csi[i][j];
        IP->eta_y[i][j] =  IP->J[i][j]*IP->x_csi[i][j];
        
        HP->csi_xhi[i][j] =  HP->J_hi[i][j]*HP->yhi_eta[i][j];
        HP->csi_yhi[i][j] = -HP->J_hi[i][j]*HP->xhi_eta[i][j];
        HP->eta_xhi[i][j] = -HP->J_hi[i][j]*HP->yhi_csi[i][j];
        HP->eta_yhi[i][j] =  HP->J_hi[i][j]*HP->xhi_csi[i][j];
        
        HP->csi_xhj[i][j] =  HP->J_hj[i][j]*HP->yhj_eta[i][j];
        HP->csi_yhj[i][j] = -HP->J_hj[i][j]*HP->xhj_eta[i][j];
        HP->eta_xhj[i][j] = -HP->J_hj[i][j]*HP->yhj_csi[i][j];
        HP->eta_yhj[i][j] =  HP->J_hj[i][j]*HP->xhj_csi[i][j];
    }
    
    void HPmetrics(HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {
        // Integer-Point Metric Parameters
        IP->A1[i][j] = IP->csi_x[i][j]*IP->csi_x[i][j] + IP->csi_y[i][j]*IP->csi_y[i][j];
        IP->A2[i][j] = IP->csi_x[i][j]*IP->eta_x[i][j] + IP->csi_y[i][j]*IP->eta_y[i][j];
        IP->A3[i][j] = IP->eta_x[i][j]*IP->eta_x[i][j] + IP->eta_y[i][j]*IP->eta_y[i][j];

        HP->A1_hi[i][j] = HP->csi_xhi[i][j]*HP->csi_xhi[i][j] + HP->csi_yhi[i][j]*HP->csi_yhi[i][j];
        HP->A2_hi[i][j] = HP->csi_xhi[i][j]*HP->eta_xhi[i][j] + HP->csi_yhi[i][j]*HP->eta_yhi[i][j];
        HP->A3_hi[i][j] = HP->eta_xhi[i][j]*HP->eta_xhi[i][j] + HP->eta_yhi[i][j]*HP->eta_yhi[i][j];

        HP->A1_hj[i][j] = HP->csi_xhj[i][j]*HP->csi_xhj[i][j] + HP->csi_yhj[i][j]*HP->csi_yhj[i][j];
        HP->A2_hj[i][j] = HP->csi_xhj[i][j]*HP->eta_xhj[i][j] + HP->csi_yhj[i][j]*HP->eta_yhj[i][j];
        HP->A3_hj[i][j] = HP->eta_xhj[i][j]*HP->eta_xhj[i][j] + HP->eta_yhj[i][j]*HP->eta_yhj[i][j]; 
    }
    
    void FreeVars(HalfPoints * HP, IntegerPoints * IP)
    {
        int i;
    
        for(i = 0; i < IMAX; i++) { 
            
            free(IP->x_csi[i]);
            free(IP->x_eta[i]);
            free(IP->y_csi[i]);
            free(IP->y_eta[i]);
            
            free(IP->csi_x[i]);
            free(IP->csi_y[i]);
            free(IP->eta_x[i]);
            free(IP->eta_y[i]);
            
            free(HP->xhi_csi[i]);
            free(HP->yhi_csi[i]);
            free(HP->xhj_eta[i]);
            free(HP->yhj_eta[i]);
                
            free(HP->xhj_csi[i]);
            free(HP->yhj_csi[i]);
            free(HP->xhi_eta[i]);
            free(HP->yhi_eta[i]);
             
            free(HP->csi_xhi[i]);
            free(HP->csi_yhi[i]);
            free(HP->eta_xhi[i]);
            free(HP->eta_yhi[i]);
             
            free(HP->csi_xhj[i]);
            free(HP->csi_yhj[i]);
            free(HP->eta_xhj[i]);
            free(HP->eta_yhj[i]);

        }
        
        // IP 
        free(IP->x_csi);   
        free(IP->x_eta);  
        free(IP->y_csi);   
        free(IP->y_eta);   
        
        free(IP->csi_x);   
        free(IP->csi_y);   
        free(IP->eta_x);   
        free(IP->eta_y);   
        
        // HP
        free(HP->xhi_csi); 
        free(HP->yhi_csi);
        free(HP->xhj_eta); 
        free(HP->yhj_eta); 
         
        free(HP->xhj_csi); 
        free(HP->yhj_csi); 
        free(HP->xhi_eta); 
        free(HP->yhi_eta); 
          
        free(HP->csi_xhi); 
        free(HP->csi_yhi); 
        free(HP->eta_xhi); 
        free(HP->eta_yhi); 
        
        free(HP->csi_xhj); 
        free(HP->csi_yhj); 
        free(HP->eta_xhj); 
        free(HP->eta_yhj); 
   
    }

    void ConstSetCal(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP)
    {
        // Setting mesh
        InverseMesh(mesh);
        
        int i=0, j=0;
        FILE * fw;
        fw = fopen("inverse_mesh.dat","w");
        for(i=0;i<=IMAX-1;i++)
        {
            for(j=0;j<=JMAX-1;j++)
            {
            	fprintf(fw,"%lf %lf\n",mesh->x[i][j],mesh->y[i][j]);
                // Calculate Integer-Point derivatives x_csi, x_eta, y_csi and y_eta
                // and Half-Point derivatives xhi_csi, xhi_eta, yhi_csi, yhi_eta
                //                            xhj_csi, xhj_eta, yhj_csi, yhj_eta
                HPderivPM(mesh,HP,IP,i,j); // Integer-Point/Half-Point derivatives on Physical Mesh
                
                // Calculate Half-Point and Integer-Point Jacobians
                //                      J, J_hi and J_hj
                HPjacobian(HP,IP,i,j);

                // Calculate Integer-Point derivatives csi_x, csi_y, eta_x and eta_y
                // and Half-Point derivatives csi_xhi, csi_yhi, eta_xhi, eta_yhi
                //                            csi_xhj, csi_yhj, eta_xhj, eta_yhj
                HPderivNM(HP,IP,i,j); // Integer-Point/Half-Point derivatives on Numerical Mesh
                
                // Calculate Half-Point and Integer-Point Metric Parameters
                //                      A1   , A2   , A3
                //                      A1_hi, A2_hi, A3_hi
                //                      A1_hj, A2_hj, A3_hj
                HPmetrics(HP,IP,i,j);
            }
        }
        fclose(fw);
        // It can be done a free() function for those unnecessary variables
        FreeVars(HP,IP);
        
        printf("\nAll pre-settings have been done successfully!\n");
        printf("We are ready to start it!\n");
        
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* ---------- Variables Settings and Calculations ------------- */
// ------------------------------------------------------------ //
    
    void HPfiDeriv(HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {

    	if(i==0 && j>0  && j<JMAX-1)
    	{
    		IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[IMAX-2][j]);
	        IP->fi_eta[i][j] = 0.5*(IP->fi[i][j+1] - IP->fi[i][j-1]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[i+1][j] - IP->fi[i][j]);
	        HP->fi_eta_hi[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i+1][j-1]) + (IP->fi[i][j+1] - IP->fi[i][j-1]));
	        
	        HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[IMAX-2][j+1]) + (IP->fi[i+1][j] - IP->fi[IMAX-2][j])); 
	        HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);			
    	}
    	else if(i>0 && j==0  && i<IMAX-1)
    	{
    		IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[i-1][j]);
	        IP->fi_eta[i][j] = 0.5*(-3.0*IP->fi[i][j] + 4.0*IP->fi[i][j+1] - IP->fi[i][j+2]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[i+1][j] - IP->fi[i][j]);
	        HP->fi_eta_hi[i][j] = 0.25*(-3.0*(IP->fi[i+1][j]+IP->fi[i][j]) +4.0*(IP->fi[i+1][j+1]+IP->fi[i][j+1]) -(IP->fi[i+1][j+2] + IP->fi[i][j+2]));
	        
	        HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[i+1][j] - IP->fi[i-1][j])); 
	        HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);	

    	}
    	else if(i==0 && j==0)
    	{
    		IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[IMAX-2][j]);
	        IP->fi_eta[i][j] = 0.5*(-3.0*IP->fi[i][j] + 4.0*IP->fi[i][j+1] - IP->fi[i][j+2]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[i+1][j] - IP->fi[i][j]);
	        HP->fi_eta_hi[i][j] = 0.25*(-3.0*(IP->fi[i+1][j]+IP->fi[i][j]) +4.0*(IP->fi[i+1][j+1]+IP->fi[i][j+1]) -(IP->fi[i+1][j+2] + IP->fi[i][j+2]));
	        
	        HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[IMAX-2][j+1]) + (IP->fi[i+1][j] - IP->fi[IMAX-2][j])); 
	        HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);		
    	}
    	else if(i==IMAX-1 && j <JMAX-1)
    	{
    		IP->fi_csi[i][j] = 0.5*(IP->fi[1][j] - IP->fi[i-1][j]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[1][j] - IP->fi[i][j]);
	        HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);			

	        if(j==0)
	        {
	        	IP->fi_eta[i][j] = 0.5*(-3.0*IP->fi[i][j] + 4.0*IP->fi[i][j+1] - IP->fi[i][j+2]);
	        	HP->fi_eta_hi[i][j] = 0.25*(-3.0*(IP->fi[1][j]+IP->fi[i][j]) +4.0*(IP->fi[1][j+1]+IP->fi[i][j+1]) -(IP->fi[1][j+2] + IP->fi[i][j+2]));
	        	HP->fi_csi_hj[i][j] = 0.25*((IP->fi[1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[1][j] - IP->fi[i-1][j])); 
	        }
	        else
	        {
	        	IP->fi_eta[i][j] = 0.5*(IP->fi[i][j+1] - IP->fi[i][j-1]);
		        HP->fi_eta_hi[i][j] = 0.25*((IP->fi[1][j+1] - IP->fi[1][j-1]) + (IP->fi[i][j+1] - IP->fi[i][j-1]));
		        HP->fi_csi_hj[i][j] = 0.25*((IP->fi[1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[1][j] - IP->fi[i-1][j])); 
		        
	        }
    	}
    	else if(i<IMAX-1 && j==JMAX-1)
    	{
	    	if(i==0)
	    	{
	    		IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[IMAX-2][j]);
	    	}
	    	else
	    	{
	    		IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[i-1][j]);
	    	}
	    	
	        IP->fi_eta[i][j] = 0.5*(3.0*IP->fi[i][j] -4.0*IP->fi[i][j-1] + IP->fi[i][j-2]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[i+1][j] - IP->fi[i][j]);
	        // HP->fi_eta_hi[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i+1][j-1]) + (IP->fi[i][j+1] - IP->fi[i][j-1]));
	        
	        // HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[i+1][j] - IP->fi[i-1][j])); 
	        // HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);	    
    	}
    	else if(i==IMAX-1 && j==JMAX-1)
    	{
    		IP->fi_csi[i][j] = 0.5*(IP->fi[1][j] - IP->fi[i-1][j]);
	        IP->fi_eta[i][j] = 0.5*(3.0*IP->fi[i][j] -4.0*IP->fi[i][j-1] + IP->fi[i][j-2]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[1][j] - IP->fi[i][j]);
	        // HP->fi_eta_hi[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i+1][j-1]) + (IP->fi[i][j+1] - IP->fi[i][j-1]));
	        
	        // HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[i+1][j] - IP->fi[i-1][j])); 
	        // HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);		
    	}
    	else
    	{
    		IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[i-1][j]);
	        IP->fi_eta[i][j] = 0.5*(IP->fi[i][j+1] - IP->fi[i][j-1]);
	        
	        HP->fi_csi_hi[i][j] = (IP->fi[i+1][j] - IP->fi[i][j]);
	        HP->fi_eta_hi[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i+1][j-1]) + (IP->fi[i][j+1] - IP->fi[i][j-1]));
	        
	        HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[i+1][j] - IP->fi[i-1][j])); 
	        HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);
    	}
        
    }

    void HPvel(HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {

        IP->U[i][j] = IP->A1[i][j]*IP->fi_csi[i][j] + IP->A2[i][j]*IP->fi_eta[i][j];
        IP->V[i][j] = IP->A2[i][j]*IP->fi_csi[i][j] + IP->A3[i][j]*IP->fi_eta[i][j];
        
        HP->U_hi[i][j] = HP->A1_hi[i][j]*HP->fi_csi_hi[i][j] + HP->A2_hi[i][j]*HP->fi_eta_hi[i][j];
        HP->U_hj[i][j] = HP->A1_hj[i][j]*HP->fi_csi_hj[i][j] + HP->A2_hj[i][j]*HP->fi_eta_hj[i][j];
        
        HP->V_hi[i][j] = HP->A2_hi[i][j]*HP->fi_csi_hi[i][j] + HP->A3_hi[i][j]*HP->fi_eta_hi[i][j];
        HP->V_hj[i][j] = HP->A2_hj[i][j]*HP->fi_csi_hj[i][j] + HP->A3_hj[i][j]*HP->fi_eta_hj[i][j];
        
    }
    
    void HPrho(HalfPoints * HP, IntegerPoints * IP, int i, int j)
    {
    	// if(j==0)
    	// {
    	// 	IP->rho[i][j]   = ;
    	// 	HP->rho_hi[i][j]= ;
    	// }
    	// else
    	// {
    		IP->rho[i][j]    = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(IP->U[i][j]*IP->fi_csi[i][j] + IP->V[i][j]*IP->fi_eta[i][j])),(1.0/(gamma-1.0)));
        
        	HP->rho_hi[i][j] = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(HP->U_hi[i][j]*HP->fi_csi_hi[i][j] + HP->V_hi[i][j]*HP->fi_eta_hi[i][j])),(1.0/(gamma-1.0)));	
    	// }
        
        printf("HP->U_hi[%d][%d]=%lf\n",i,j,HP->U_hi[i][j]);
        printf("HP->V_hi[%d][%d]=%lf\n",i,j,HP->V_hi[i][j]);
		printf("HP->fi_csi_hi[%d][%d]=%lf\n",i,j,HP->fi_csi_hi[i][j]);
        printf("HP->fi_eta_hi[%d][%d]=%lf\n",i,j,HP->fi_eta_hi[i][j]);
        printf("HP->rho_hi[%d][%d]=%lf\n",i,j,HP->rho_hi[i][j]);
        printf("IP->rho[%d][%d]=%lf\n",i,j,IP->rho[i][j]);
printf("calc=%lf\n",pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(HP->U_hi[i][j]*HP->fi_csi_hi[i][j] + HP->V_hi[i][j]*HP->fi_eta_hi[i][j])),(1.0/(gamma-1.0))));
printf("calc2=%lf\n",pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(22.91484091469)),(1.0/(gamma-1.0))));

        printf("\n");
        // HP->rho_hj[i][j] = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(IP->U_hj[i][j]*IP->fi_csi_hj[i][j] + IP->V_hj[i][j]*IP->fi_eta_hj[i][j])),(1.0/(gamma-1.0)));
    }
    
    void HPrhoArtf(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0, r=0 ,s=0;
        double C1,C2,C;
        
        C1 = pow(2.0/(gamma+1.0),(1.0/(gamma-1.0)));
        C2 = 4.9325;
        C  = 1.0;
        for(i=0;i<=IMAX-1;i++)
        {
            for(j=0;j<= JMAX-1;j++)
            {
                if (HP->U_hi[i][j] < 0.0)
                {
                    r = 1;
                    // HP->nu_hi[i][j] = ((IP->M[i+1][j]*IP->M[i+1][j]-1.0)*C > 0.0) ? (IP->M[i+1][j]*IP->M[i+1][j]-1.0)*C : 0.0;
                    if(i==IMAX-1)
                    {
                    	HP->nu_hi[i][j] = ((C1-IP->rho[1][j])*C2*C > 0.0) ? (C1-IP->rho[1][j])*C2*C : 0.0;
                    }
                    else
                    {
                    	HP->nu_hi[i][j] = ((C1-IP->rho[i+1][j])*C2*C > 0.0) ? (C1-IP->rho[i+1][j])*C2*C : 0.0;
                    } 
                    	
                }
                else
                {
                    r = -1;
                    // HP->nu_hi[i][j] = ((IP->M[i][j]*IP->M[i][j]-1.0)*C > 0.0) ? (IP->M[i][j]*IP->M[i][j]-1.0)*C : 0.0;
                    HP->nu_hi[i][j] = ((C1-IP->rho[i][j])*C2*C > 0.0) ? (C1-IP->rho[i][j])*C2*C : 0.0;
                }
                
                if(i==IMAX-1 && r == 1)
                {
                	HP->rho_til_hi[i][j] = (1.0-HP->nu_hi[i][j])*HP->rho_hi[i][j] + HP->nu_hi[i][j]*HP->rho_hi[1][j];
                }
                else if(i==0 && r == -1)
                {
                	HP->rho_til_hi[i][j] = (1.0-HP->nu_hi[i][j])*HP->rho_hi[i][j] + HP->nu_hi[i][j]*HP->rho_hi[IMAX-2][j];
                }
                else
                {
                	HP->rho_til_hi[i][j] = (1.0-HP->nu_hi[i][j])*HP->rho_hi[i][j] + HP->nu_hi[i][j]*HP->rho_hi[i+r][j];
                } 
                
            }
        }
        for(i=0;i<=IMAX-1;i++)
        {
            for(j=0;j<=JMAX-1;j++)
            {
                if (HP->V_hj[i][j] < 0.0)
                {
                    s = 1;
                    // HP->nu_hj[i][j] = ((IP->M[i][j+1]*IP->M[i][j+1]-1.0)*C > 0.0) ? (IP->M[i][j+1]*IP->M[i][j+1]-1.0)*C : 0.0;
                    
                    if(j==JMAX-1)
                    {
                    	// HP->nu_hj[i][j] = ((C1-IP->rho[i][j+1])*C2*C > 0.0) ? (C1-IP->rho[i][j+1])*C2*C : 0.0;
                    }
                    else
                    {
                    	HP->nu_hj[i][j] = ((C1-IP->rho[i][j+1])*C2*C > 0.0) ? (C1-IP->rho[i][j+1])*C2*C : 0.0;
                    } 
                    	
                }
                else
                {
                    s = -1;
                    // HP->nu_hj[i][j] = ((IP->M[i][j]*IP->M[i][j]-1.0)*C > 0.0) ? (IP->M[i][j]*IP->M[i][j]-1.0)*C : 0.0;
                    HP->nu_hj[i][j] = ((C1-IP->rho[i][j])*C2*C > 0.0) ? (C1-IP->rho[i][j])*C2*C : 0.0;
                }
                

                if(j==JMAX-1 && s == 1)
                {
                	// HP->rho_bar_hj[i][JMAX-1] = HP->rho_bar_hj[i][JMAX-2];
                	HP->rho_bar_hj[i][j] = (1.0-HP->nu_hj[i][j])*HP->rho_hj[i][j] + HP->nu_hj[i][j]*HP->rho_hj[i][JMAX-2];
                }
                else if(j==0 && s == -1)
                {
                	// HP->rho_bar_hj[i][j] = (1.0-HP->nu_hj[i][j])*HP->rho_hj[i][j] + HP->nu_hj[i][j]*HP->rho_hj[i][j+s];
                }
                else
                {
                	HP->rho_bar_hj[i][j] = (1.0-HP->nu_hj[i][j])*HP->rho_hj[i][j] + HP->nu_hj[i][j]*HP->rho_hj[i][j+s];
					// printf("HP->rho_hj[%d][%d] = %lf\n",i,j,HP->rho_hj[i][j]);
					// printf("HP->rho_hi[%d][%d] = %lf\n",i,j,HP->rho_hi[i][j]);
					// printf("HP->rho_hj[%d][%d+%d] = %lf\n",i,j,s,HP->rho_hj[i][j+s]);
					// printf("HP->nu_hj[%d][%d] = %lf\n",i,j,HP->nu_hj[i][j]);
     //            	printf("HP->rho_bar_hj[%d][%d] = %lf\n",i,j,HP->rho_bar_hj[i][j]);
                } 
            }
        }
    }
    void VarSetCal(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=0;i<=IMAX-1;i++)
        {
            for(j=0;j<=JMAX-1;j++)
            {
                // Calculate Half-Point Fi derivatives - fi_csi and fi_eta
                HPfiDeriv(HP,IP,i,j);
                
                // Calculate Half-Point Velocities
                HPvel(HP,IP,i,j);
                
                // Calculate Half-Point Densities
                HPrho(HP,IP,i,j);
            }
        }
        for(i=0;i<=IMAX-1;i++)
        {
            for(j=0;j<=JMAX-1;j++)
            { 
            	// Class Note Suggestion:
            	if(i==0 && j<JMAX-1)
                {
                	HP->rho_hj[i][j] = 0.25*(HP->rho_hi[i][j] + HP->rho_hi[i][j+1] + HP->rho_hi[IMAX-2][j+1] + HP->rho_hi[IMAX-2][j]);
                }
                else if(i>0 && j==JMAX-1)
                {
                	// HP->rho_hj[i][j] = 0.25*(HP->rho_hi[i][j] + HP->rho_hi[i][j+1] + HP->rho_hi[i-1][j+1] + HP->rho_hi[i-1][j]);
                } 
                else if(i==0 && j==JMAX-1)
                {
                	// HP->rho_hj[i][j] = 0.25*(HP->rho_hi[i][j] + HP->rho_hi[i][j+1] + HP->rho_hi[IMAX-2][j+1] + HP->rho_hi[IMAX-2][j]);
                } 
                else
                {
                	HP->rho_hj[i][j] = 0.25*(HP->rho_hi[i][j] + HP->rho_hi[i][j+1] + HP->rho_hi[i-1][j+1] + HP->rho_hi[i-1][j]);
                }
		        // Class Note Suggestion:
		        // HP->rho_hj[i][j] = 0.25*(HP->rho_hi[i][j] + HP->rho_hi[i][j+1] + HP->rho_hi[i-1][j+1] + HP->rho_hi[i-1][j]);
            }
        	
        }
        for(j=0;j<=JMAX-1;j++)
        { 
			HP->rho_hj[0][j] = HP->rho_hj[IMAX-1][j];
		}

        // Calculate Half-Point Artificial Densities
        HPrhoArtf(HP,IP);
    }
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Residue Calculation ------------------- */
// ------------------------------------------------------------ //
    void Residue(HalfPoints * HP, IntegerPoints * IP, int n)
    {

    	int i=0,j=0;
    	if(n>0)
    	{
    		IP->resMax = (double *) realloc(IP->resMax, n+1);
	    	IP->resL1  = (double *) realloc(IP->resL1 , n+1);
	    	IP->resL2  = (double *) realloc(IP->resL2 , n+1);
    	}
    	
        IP->resMax[n] = 0.0;
        IP->resL1[n]  = 0.0;
        IP->resL2[n]  = 0.0;

        for(i=1;i<=IMAX-1;i++)
        {
            for(j=1;j<=JMAX-1;j++)
            {
            	// printf("J[%d][%d]=%lf;  J_hi[%d][%d]=%lf;  J_hj[%d][%d]=%lf\n",i,j,IP->J[i][j],i,j,HP->J_hi[i][j],i,j,HP->J_hj[i][j]);
                if(i==0 && j<JMAX-1)
            	{
            		// printf("J_hi[%d][%d]=%lf;  J_hj[%d][%d]=%lf\n",i,j,HP->J_hi[i][j],i,j,HP->J_hj[i][j]);
            		// printf("J_hi[%d][%d]=%lf;  J_hj[%d][%d]=%lf\n",IMAX-2,j,HP->J_hi[IMAX-2][j],i,j-1,HP->J_hj[i][j-1]);
	             //    printf("HP->rho_til_hi[%d][%d]=%lf;  HP->rho_bar_hj[%d][%d]=%lf\n",i,j,HP->rho_til_hi[i][j],i,j,HP->rho_bar_hj[i][j]);
	             //    printf("HP->rho_til_hi[%d][%d]=%lf;  HP->rho_bar_hj[%d][%d]=%lf\n",IMAX-2,j,HP->rho_til_hi[IMAX-2][j],i,j-1,HP->rho_bar_hj[i][j-1]);
	             //    printf("HP->U_hi[%d][%d]=%lf;  HP->V_hj[%d][%d]=%lf\n",i,j,HP->U_hi[i][j],i,j,HP->V_hj[i][j]);
	             //    printf("HP->U_hi[%d][%d]=%lf;  HP->V_hj[%d][%d]=%lf\n",IMAX-2,j,HP->U_hi[IMAX-2][j],i,j-1,HP->V_hj[i][j-1]);

            		IP->residue[i][j] = ((HP->rho_til_hi[i][j]*HP->U_hi[i][j]/HP->J_hi[i][j]) - (HP->rho_til_hi[IMAX-2][j]*HP->U_hi[IMAX-2][j]/HP->J_hi[IMAX-2][j])) + ((HP->rho_bar_hj[i][j]*HP->V_hj[i][j]/HP->J_hj[i][j]) - (HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]));
            	}
            	else if(j==JMAX-1 && i>0)
            	{
            		IP->residue[i][j] = ((HP->rho_til_hi[i][j]*HP->U_hi[i][j]/HP->J_hi[i][j]) - (HP->rho_til_hi[i-1][j]*HP->U_hi[i-1][j]/HP->J_hi[i-1][j])) + ((-HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]) - (HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]));
            	}
            	else if(j==JMAX-1 && i==0)
            	{
            		IP->residue[i][j] = ((HP->rho_til_hi[i][j]*HP->U_hi[i][j]/HP->J_hi[i][j]) - (HP->rho_til_hi[IMAX-2][j]*HP->U_hi[IMAX-2][j]/HP->J_hi[IMAX-2][j])) + ((-HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]) - (HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]));
            	}
            	else
            	{
            		IP->residue[i][j] = ((HP->rho_til_hi[i][j]*HP->U_hi[i][j]/HP->J_hi[i][j]) - (HP->rho_til_hi[i-1][j]*HP->U_hi[i-1][j]/HP->J_hi[i-1][j])) + ((HP->rho_bar_hj[i][j]*HP->V_hj[i][j]/HP->J_hj[i][j]) - (HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]));
            	}
                IP->resMax[n] = (IP->resMax[n] < fabs(IP->residue[i][j]))? fabs(IP->residue[i][j]) : IP->resMax[n];
                IP->resL1[n]  = IP->resL1[n]  + fabs(IP->residue[i][j]);
                IP->resL2[n]  = IP->resL2[n]  + IP->residue[i][j]*IP->residue[i][j];
				
				printf("J_hi[%d][%d]=%lf;  J_hj[%d][%d]=%lf\n",i,j,HP->J_hi[i][j],i,j,HP->J_hj[i][j]);
        		// printf("J_hi[%d][%d]=%lf;  J_hj[%d][%d]=%lf\n",i-1,j,HP->J_hi[i-1][j],i,j-1,HP->J_hj[i][j-1]);
                printf("HP->rho_til_hi[%d][%d]=%lf;  HP->rho_bar_hj[%d][%d]=%lf\n",i,j,HP->rho_til_hi[i][j],i,j,HP->rho_bar_hj[i][j]);
                // printf("HP->rho_til_hi[%d][%d]=%lf;  HP->rho_bar_hj[%d][%d]=%lf\n",i-1,j,HP->rho_til_hi[i-1][j],i,j-1,HP->rho_bar_hj[i][j-1]);
                printf("HP->U_hi[%d][%d]=%lf;  HP->V_hj[%d][%d]=%lf\n",i,j,HP->U_hi[i][j],i,j,HP->V_hj[i][j]);
                // printf("HP->U_hi[%d][%d]=%lf;  HP->V_hj[%d][%d]=%lf\n",i-1,j,HP->U_hi[i-1][j],i,j-1,HP->V_hj[i][j-1]);                

                printf("residuo[%d][%d] = %lf\n",i,j,IP->residue[i][j]);
            }
        }
        // IP->resMax = log10(IP->resMax);
        IP->resL2[n] = sqrt(IP->resL2[n]);
        printf("residuoL2 = %lf; residuoMax = %lf\n",IP->resL2[n],IP->resMax[n]);
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Iteration Method ---------------------- */
// ------------------------------------------------------------ //
    
    // Thomas Algorithm
    void Thomas_Algorithm(double ** C , ThomStr * T, int k, int dir, int I, int J)
    {
        int i,j;

        if(dir==1)  // line on x
        {
            for(i=0;i<I;i++)
            {   
                if(i==0)
                {
                    T->c[i] = (double) T->c[i]/T->b[i];
                    T->RHS[i] = (double) T->RHS[i]/T->b[i];
                }
                else
                {
                    T->c[i] = (double) (T->c[i]/(T->b[i]-T->a[i]*T->c[i-1]));
                    T->RHS[i] = (double)((T->RHS[i]-T->a[i]*T->RHS[i-1])/(T->b[i]-T->a[i]*T->c[i-1]));
                }
            }
            for(i=I-1;i>-1;i--)
            {   
                if(i==I-1)
                {
                    C[i][k] = T->RHS[i];
                }
                else
                {
                    C[i][k] = T->RHS[i]-T->c[i]*C[i+1][k];
                }
            }
        }
        else if(dir==2) // line on y
        {
            for(j=1;j<J-1;j++)   
            {
                if(j==1)
                {
                    T->c[j] = (double)(T->c[j]/T->b[j]);
                    T->RHS[j] = (double)(T->RHS[j]/T->b[j]);
                }
                else
                {
                    T->c[j] = (double)(T->c[j]/(T->b[j]-T->a[j]*T->c[j-1]));
                    T->RHS[j] = (double)((T->RHS[j]-T->a[j]*T->RHS[j-1])/(T->b[j]-T->a[j]*T->c[j-1]));
                }
            }
            for(j=J-2;j>0;j--)
            {   
                if(j==J-2)
                {
                    C[k][j] = T->RHS[j];
                }
                else
                {
                    C[k][j] = T->RHS[j]-T->c[j]*C[k][j+1];
                }
            }
        }
        else if(dir==3)  // line on x
        {
            for(i=1;i<I-1;i++)
            {   
                if(i==1)
                {
                    T->c[i] = (double) T->c[i]/T->b[i];
                    T->RHS[i] = (double) T->RHS[i]/T->b[i];
                }
                else
                {
                    T->c[i] = (double) (T->c[i]/(T->b[i]-T->a[i]*T->c[i-1]));
                    T->RHS[i] = (double)((T->RHS[i]-T->a[i]*T->RHS[i-1])/(T->b[i]-T->a[i]*T->c[i-1]));
                }
            }
            for(i=I-2;i>0;i--)
            {   
                if(i==I-2)
                {
                    C[i][k] = T->RHS[i];
                }
                else
                {
                    C[i][k] = T->RHS[i]-T->c[i]*C[i+1][k];
                }
            }
        }
    }
    
    // AF2
    void IterMethod(IntegerPoints * IP, HalfPoints * HP, ThomStr * T, int iter)
    {
        int i, j;

        double Ai,Ai1;
        double Aj,Aj1;
        double alpha_n;
        double alp1H, alp2H, alp1L, alp2L;
        double Rt, betaH, betaL;
        double ** f;

        alp1H = 1.0;
        alp1L = 0.003;
        alp2H = 1.0;
        alp2L = 0.003;

        int M=10;
        f = calloc((IMAX+2), sizeof(double *));
        for(i=0;i<IMAX+2;i++)
        {
            f[i] = calloc((JMAX+1), sizeof(double)); 
        }
        
        // Step 1:
        // thomas in y-direction
        alpha_n = alp2H*pow((alp2L/alp2H),((iter-1)%M)/(M-1));
        for(i=0;i<IMAX;i++)
        {
            f[i][JMAX] = 0.0;
            for(j=JMAX-1;j>0;j--)
            {
                Aj  = HP->rho_bar_hj[i][j-1]*HP->A3_hj[i][j-1]/HP->J_hj[i][j-1];
                Aj1 = HP->rho_bar_hj[i][j]*HP->A3_hj[i][j]/HP->J_hj[i][j];
                
                // f[i][j] = (alpha_n*omega*IP->residue[i][j] - Aj1*f[i][j-1])/(alpha_n - Aj);
                f[i][j] = (alpha_n*omega*IP->residue[i][j] + Aj1*f[i][j+1])/(alpha_n - Aj);

                // 0        IMAX-2    IMAX-2
                // 1        0         0
                // 2        1         1
                // 3        2         2
                // ...      ...       ...
                // IMAX-1   IMAX-2    IMAX-2
                // IMAX     IMAX-1    0

                // Boundary Condition ?
            }
        }

        if(iter>M)
    	{
    		Rt = (IP->resL2[iter]/IP->resL2[iter-M]) + (IP->resMax[iter]/IP->resMax[iter-M]);
    	}

    	if(iter>0) IP->beta = (double *) realloc(IP->beta, iter+1);

        // Step 2:
        // thomas in x-direction
        for(j=1;j<=JMAX-1;j++)
        {
            for(i=0;i<IMAX;i++)
            {	
            	if(i==0)
            	{
            		Ai  = HP->rho_til_hi[IMAX-2][j]*HP->A1_hi[IMAX-2][j]/HP->J_hi[IMAX-2][j];
	                Ai1 = HP->rho_til_hi[i][j]*HP->A1_hi[i][j]/HP->J_hi[i][j];
	            
	                T->RHS[i] = f[i][j] + alpha_n*IP->corr[i][j-1];
            	}
            	else
            	{
            		Ai  = HP->rho_til_hi[i-1][j]*HP->A1_hi[i-1][j]/HP->J_hi[i-1][j];
	                Ai1 = HP->rho_til_hi[i][j]*HP->A1_hi[i][j]/HP->J_hi[i][j];
	            
	                T->RHS[i] = f[i][j] + alpha_n*IP->corr[i][j-1];
            	}

            	if(IP->rho[i][j] > pow((2.0/(gamma-1.0)),(1.0/(gamma-1.0))))
            	{
            	// Supersonic region
            		IP->beta[iter]=4.5;
            		betaH = IP->beta[0] +1.0;
    				betaL = IP->beta[0] -1.0;
            		
            		if(iter>0) 
            		{
            			if(Rt < 2.0)
			    		{
			    			IP->beta[iter] = 0.98*IP->beta[iter-1];
			    		}
			    		else if(Rt > 2.1)
			    		{
			    			IP->beta[iter] = 1.1*IP->beta[iter-1];
			    		}

			    		if(IP->beta[iter] > betaH)
			    		{
			    			IP->beta[iter] = betaH;
			    		}
			    		else if(IP->beta[iter] < betaL)
			    		{
			    			IP->beta[iter] = betaL;
			    		}
            		}
            		
            		T->a[i]= -Ai1;
	                T->b[i]= Ai + Ai1 + alpha_n + alpha_n*IP->beta[iter];
	                T->c[i]= -Ai - alpha_n*IP->beta[iter];
            	}
            	else
            	{
            	// Subsonic region
            		IP->beta[iter]=0.3;

            		T->a[i]= -Ai1 - alpha_n*IP->beta[iter];
	                T->b[i]= Ai + Ai1 + alpha_n + alpha_n*IP->beta[iter];
	                T->c[i]= -Ai;
            	}    
            }
            // Thomas_Algorithm(IP->corr,T,j,3,IMAX+1,JMAX);
            Thomas_Algorithm(IP->corr,T,j,3,IMAX+1,JMAX);
        }

        for(i=0;i<IMAX+2;i++)
        {
            free(f[i]);        
        }

        free(f);
            
    }
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Boundary Conditions ------------------- */
// ------------------------------------------------------------ //
    void BCs(IntegerPoints * IP, HalfPoints * HP, MeshGrid * mesh)
    {
        int i=0,j=0;
        double fi_inf, U_INF;
    
        U_INF = pow((gamma+1.0)/(gamma-1.0+2.0/(IP->M_INF*IP->M_INF)),0.5);
    
        for(i=0;i<IMAX;i++)
        {
            // Airfoil Boundary - Wall Condition
            HP->V_hj[i][JMAX-1] = -HP->V_hj[i][JMAX-2];
            HP->rho_bar_hj[i][JMAX-1] = HP->rho_bar_hj[i][JMAX-2];
            HP->J_hj[i][JMAX-1] = HP->J_hj[i][JMAX-2];
            
            HP->U_hj[i][JMAX-1] = HP->U_hj[i][JMAX-2];
            


            // HP->V_hi[i][JMAX-1] = -HP->V_hi[i][JMAX-2];
            // HP->rho_til_hi[i][JMAX-1] = HP->rho_til_hi[i][JMAX-2];
            // HP->J_hi[i][JMAX-1] = HP->J_hi[i][JMAX-2];
            
            // HP->U_hi[i][JMAX-1] = HP->U_hi[i][JMAX-2];


            IP->V[i][JMAX-1] = 0.0;
            IP->fi_eta[i][JMAX-1] = - (IP->A2[i][JMAX-1]*IP->fi_csi[i][JMAX-1]/IP->A3[i][JMAX-1]);        

            IP->rho[i][0] = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(U_INF*U_INF)),(1.0/(gamma-1.0)));
            HP->rho_hi[i][0] = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(U_INF*U_INF)),(1.0/(gamma-1.0)));


            // External Boundary Condition
            fi_inf = U_INF*mesh->x[i][0];
            IP->fi[i][0] = fi_inf;

            // IP->U[i][JMAX-1] = IP->U[i][JMAX-2];
            // IP->fi[i][JMAX-1] = IP->fi[i][JMAX-2];
            
        }   
        
        for(j=0;j<JMAX;j++)
        {
            // Weak Boundary Condition
            IP->fi[0][j] = IP->fi[IMAX-1][j];
        }
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Printing Functions -------------------- */
// ------------------------------------------------------------ //

    // Printing when converges
    void PrintConv(int n, IntegerPoints * IP)
    {
        printf("\nConverged at %d iteration!! Residue = %f\n",n,IP->resMax[n]);
    }
    
    // Printing when diverges
    void PrintDiverg(int n, IntegerPoints * IP)
    {
        printf("\nDiverged at %d iteration!! Residue = %f\n",n,IP->resMax[n]);
    }
    
    // Printing iteration
    void PrintIter(int n, IntegerPoints * IP, double timeIter)
    {
        printf("\nIteration = %d ..... Residue = %e .... Time = %lf \
                seconds\n",n,IP->resMax[n],timeIter);
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------------- Free-up ------------------------- */
// ------------------------------------------------------------ //
    
    
    void FreeUp(IntegerPoints * IP, HalfPoints * HP, MeshGrid * mesh, ThomStr * T)
    {
        int i;
    
        for(i = 0; i < IMAX; i++) { 
            free(mesh->x[i]);
            free(mesh->y[i]);
            free(IP->resMax);
            free(IP->resL1);
			free(IP->resL2);
			free(IP->beta);
            free(IP->residue[i]);
            free(IP->fi[i]);
            free(IP->corr[i]);
            free(IP->rho[i]);
            // free(IP->x_csi[i]);
            // free(IP->x_eta[i]);
            // free(IP->y_csi[i]);
            // free(IP->y_eta[i]);
            
            // free(IP->csi_x[i]);
            // free(IP->csi_y[i]);
            // free(IP->eta_x[i]);
            // free(IP->eta_y[i]);
            free(IP->J[i]);
            free(IP->A1[i]);
            free(IP->A2[i]);
            free(IP->A3[i]);
              
            free(IP->fi_csi[i]);
            free(IP->fi_eta[i]);
            free(IP->U[i]);
            free(IP->V[i]);
            
            // free(HP->xhi_csi[i]);
            // free(HP->yhi_csi[i]);
            // free(HP->xhj_eta[i]);
            // free(HP->yhj_eta[i]);
                
            // free(HP->xhj_csi[i]);
            // free(HP->yhj_csi[i]);
            // free(HP->xhi_eta[i]);
            // free(HP->yhi_eta[i]);
             
            // free(HP->csi_xhi[i]);
            // free(HP->csi_yhi[i]);
            // free(HP->eta_xhi[i]);
            // free(HP->eta_yhi[i]);
             
            // free(HP->csi_xhj[i]);
            // free(HP->csi_yhj[i]);
            // free(HP->eta_xhj[i]);
            // free(HP->eta_yhj[i]);
            
            free(HP->A1_hi[i]);
            free(HP->A2_hi[i]);
            free(HP->A3_hi[i]);
             
            free(HP->A1_hj[i]);
            free(HP->A2_hj[i]);
            free(HP->A3_hj[i]);
        
            free(HP->J_hi[i]);
            free(HP->J_hj[i]);
            
            free(HP->fi_csi_hi[i]);
            free(HP->fi_csi_hj[i]);
            free(HP->fi_eta_hi[i]);
            free(HP->fi_eta_hj[i]);
            
            free(HP->U_hi[i]);
            free(HP->U_hj[i]);
            free(HP->V_hi[i]);
            free(HP->V_hj[i]);
            
            free(HP->rho_hi[i]);
            free(HP->rho_hj[i]);
            
            free(HP->rho_til_hi[i]);
            free(HP->rho_bar_hj[i]);
            
            free(HP->nu_hi[i]);
            free(HP->nu_hj[i]);

        }
        // mesh
        free(mesh->x);     
        free(mesh->y);     
        // IP
        free(IP->residue); 
        free(IP->fi);      
        free(IP->corr);    
        free(IP->rho);     
        // free(IP->x_csi);   
        // free(IP->x_eta);  
        // free(IP->y_csi);   
        // free(IP->y_eta);   
        
        // free(IP->csi_x);   
        // free(IP->csi_y);   
        // free(IP->eta_x);   
        // free(IP->eta_y);   
        free(IP->J);       
        free(IP->A1);      
        free(IP->A2);      
        free(IP->A3);      
        
        free(IP->fi_csi);  
        free(IP->fi_eta);  
        free(IP->U);       
        free(IP->V);       
        
        // HP
        // free(HP->xhi_csi); 
        // free(HP->yhi_csi);
        // free(HP->xhj_eta); 
        // free(HP->yhj_eta); 
         
        // free(HP->xhj_csi); 
        // free(HP->yhj_csi); 
        // free(HP->xhi_eta); 
        // free(HP->yhi_eta); 
          
        // free(HP->csi_xhi); 
        // free(HP->csi_yhi); 
        // free(HP->eta_xhi); 
        // free(HP->eta_yhi); 
        
        // free(HP->csi_xhj); 
        // free(HP->csi_yhj); 
        // free(HP->eta_xhj); 
        // free(HP->eta_yhj); 
        
        free(HP->A1_hi);   
        free(HP->A2_hi);   
        free(HP->A3_hi);   
      
        free(HP->A1_hj);   
        free(HP->A2_hj);   
        free(HP->A3_hj);   
   
        free(HP->J_hi);    
        free(HP->J_hj);    
        
        free(HP->fi_csi_hi);  
        free(HP->fi_csi_hj);  
        free(HP->fi_eta_hi);  
        free(HP->fi_eta_hj);  
        
        free(HP->U_hi);       
        free(HP->U_hj);       
        free(HP->V_hi);       
        free(HP->V_hj);       
        
        free(HP->rho_hi);     
        free(HP->rho_hj);     
      
        free(HP->rho_til_hi); 
        free(HP->rho_bar_hj); 
        
        free(HP->nu_hi);      
        free(HP->nu_hj);      
        
        // T
        free(T->a);           
        free(T->b);           
        free(T->c);           
        free(T->RHS);         
        free(T->LHS); 
        
        
        free(mesh);
        free(IP);
        free(HP);
        free(T);
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------------- Output -------------------------- */
// ------------------------------------------------------------ //
    void Output(MeshGrid * mesh, IntegerPoints * IP)
    {
        int i,j;
        FILE * fw;
        fw = fopen("result.dat","w");
        
        for(i=0;i<IMAX;i++)
        {
            for(j=0;j<JMAX;j++)
            {
                fprintf(fw,"%lf   %lf   %lf   %lf   %lf\n",mesh->x[i][j],mesh->y[i][j],IP->rho[i][j],IP->U[i][j],IP->V[i][j]);
            }
        }
        fclose(fw);
    }
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //
    