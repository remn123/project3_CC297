#include <stdio.h>
#include "structs.h"

/* This file contains all functions that were used on project 3 */
/* ------------------------------------------------------------ */

/* ----------------------- Verifying Input -------------------- */
// ------------------------------------------------------------ //
    int VerInput( int argc, char ** argv)
    {
        // bla bla bla
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* --------------------- Allocating Memory -------------------- */
// ------------------------------------------------------------ //
    void AllocateAll(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP)
    {
        // bla bla bla
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* --------------------- Common Functions --------------------- */
// ------------------------------------------------------------ //

    // Centered second-order finite-differences
    double d_csi(double ** f,int i, int j)
    {
        double delta=0.0; 
        
        delta = 0.5*(f[i+1][j] - f[i-1][j]);
        
        return delta;
    }
    
    double d_eta(double ** f,int i, int j)
    {
        double delta=0.0; 
        
        delta = 0.5*(f[i][j+1] - f[i][j-1]);
        
        return delta;
    }
    
    // Backward finite-differences
    double b_csi(double ** f, int i, int j)
    {
        double delta=0.0;
        
        delta = f[i][j] - f[i-1][j];
        
        return delta;
    }
    
    double b_eta(double ** f, int i, int j)
    {
        double delta=0.0;
        
        delta = f[i][j] - f[i][j-1];
        
        return delta;
    }
    
    // Forward finite-differences
    double f_csi(double ** f, int i, int j)
    {
        double delta=0.0;
        
        delta = f[i+1][j] - f[i][j];
        
        return delta;
    }
    
    double f_eta(double ** f, int i, int j)
    {
        double delta=0.0;
        
        delta = f[i][j+1] - f[i][j];
        
        return delta;
    }
    
    // Averages
    double a_csi(double ** f, int i, int j)
    {
        double average=0.0;
        
        average = 0.5*(f[i+1][j] + f[i][j]);
        
        return average;
    }
    
    double a_eta(double ** f, int i, int j)
    {
        double average=0.0;
        
        average = 0.5*(f[i][j+1] + f[i][j]);
        
        return average;
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* --------------------- Read Mesh ---------------------------- */
// ------------------------------------------------------------ //
    void ReadMesh(MeshGrid * mesh, char * filename)
    {
        int nelm=0, i=0, j=0;
        double x=0.0, y=0.0;
        
        FILE * fr;
        
        fr = fopen(filename,"r");
        fscanf(fr,"%d %d",&IMAX,&JMAX); // reading the mesh size on global variables
        
        while(fscanf(fr,"%f %f",&x,&y)!=EoF)
        {
            mesh->x[i][j] = x;
            mesh->y[i][j] = y;
            
            if(i==IMAX-1)
            {
                j++;
                i=0;
            }
            else
            {
                i++;
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
            for(j=0;j<(JMAX+1)/2;j++)
            {
                auX = mesh->x[i][j];
                mesh->x[i][j] = mesh->x[IMAX-1+i][j];
                mesh->x[IMAX-1+i][j] = auX;
                
                auY = mesh->y[i][j];
                mesh->y[i][j] = mesh->y[IMAX-1+i][j];
                mesh->y[IMAX-1+i][j] = auY;
            }
        }
        
        // Inverting eta direction
        for(i=0;i<(IMAX+1)/2;i++)
        {
            for(j=0;j<(JMAX+1)/2;j++)
            {
                auX = mesh->x[i][j];
                mesh->x[i][j] = mesh->x[i][JMAX-1+j];
                mesh->x[i][JMAX-1+j] = auX;
                
                auY = mesh->y[i][j];
                mesh->y[i][j] = mesh->y[i][JMAX-1+j];
                mesh->y[i][JMAX-1+j] = auY;
            }
        }
    }
    
    void HPderivPM(MeshGrid * mesh, IntegerPoints * IP, HalfPoints * HP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
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
    }
    
    void HPjacobian(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
               IP->J[i][j] = 1.0/(IP->x_csi[i][j]*IP->y_eta[i][j] - IP->x_eta[i][j]*IP->y_csi[i][j]);
               HP->J_hi[i][j] = 1.0/(HP->xhi_csi[i][j]*HP->yhi_eta[i][j] - HP->xhi_eta[i][j]*HP->yhi_csi[i][j]);
               HP->J_hj[i][j] = 1.0/(HP->xhj_csi[i][j]*HP->yhj_eta[i][j] - HP->xhj_eta[i][j]*HP->yhj_csi[i][j]);
            }
        }
        
        
    }
    
    void HPderivNM(IntegerPoints * IP, HalfPoints * HP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
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
        }
    }
    
    void HPmetrics(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        // Integer-Point Metric Parameters
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
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
        }
        
    }
    
    void ConstSetCal(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP)
    {
        // Setting mesh
        InverseMesh(mesh);
        
        // Calculate Integer-Point derivatives x_csi, x_eta, y_csi and y_eta
        // and Half-Point derivatives xhi_csi, xhi_eta, yhi_csi, yhi_eta
        //                            xhj_csi, xhj_eta, yhj_csi, yhj_eta
        HPderivPM(mesh,IP,HP); // Integer-Point/Half-Point derivatives on Physical Mesh
        
        // Calculate Half-Point and Integer-Point Jacobians
        //                      J, J_hi and J_hj
        HPjacobian(HP,IP);
        
        // Calculate Integer-Point derivatives csi_x, csi_y, eta_x and eta_y
        // and Half-Point derivatives csi_xhi, csi_yhi, eta_xhi, eta_yhi
        //                            csi_xhj, csi_yhj, eta_xhj, eta_yhj
        HPderivNM(IP,HP); // Integer-Point/Half-Point derivatives on Numerical Mesh
        
        // Calculate Half-Point and Integer-Point Metric Parameters
        //                      A1   , A2   , A3
        //                      A1_hi, A2_hi, A3_hi
        //                      A1_hj, A2_hj, A3_hj
        HPmetrics(HP,IP);
        
        // It can be done a free() function for those unnecessary variables
        // freeVar(HP,IP);
        
        printf("\nAll pre-settings have been done successfully!\n");
        printf("We are ready to start it!\n");
        
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* ---------- Variables Settings and Calculations ------------- */
// ------------------------------------------------------------ //
    
    void HPfiDeriv(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                IP->fi_csi[i][j] = 0.5*(IP->fi[i+1][j] - IP->fi[i-1][j]);
                IP->fi_eta[i][j] = 0.5*(IP->fi[i][j+1] - IP->fi[i][j-1]);
                
                HP->fi_csi_hi[i][j] = (IP->fi[i+1][j] - IP->fi[i][j]);
                HP->fi_eta_hi[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i+1][j-1]) + (IP->fi[i][j+1] - IP->fi[i][j-1]));
                
                HP->fi_csi_hj[i][j] = 0.25*((IP->fi[i+1][j+1] - IP->fi[i-1][j+1]) + (IP->fi[i+1][j] - IP->fi[i-1][j])); 
                HP->fi_eta_hj[i][j] = (IP->fi[i][j+1] - IP->fi_csi[i][j]);
                
            }
        }
    }
    void HPvel(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                IP->U[i][j] = IP->A1[i][j]*IP->fi_csi[i][j] + IP->A2[i][j]*IP->fi_eta[i][j];
                IP->V[i][j] = IP->A2[i][j]*IP->fi_csi[i][j] + IP->A3[i][j]*IP->fi_eta[i][j];
                
                HP->U_hi[i][j] = HP->A1_hi[i][j]*HP->fi_csi_hi[i][j] + HP->A2_hi[i][j]*HP->fi_eta_hi[i][j];
                HP->V_hj[i][j] = HP->A2_hj[i][j]*HP->fi_csi_hj[i][j] + HP->A3_hj[i][j]*HP->fi_eta_hj[i][j];
                
            }
        }
    }
    void HPrho(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                IP->rho[i][j]    = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(IP->U[i][j]*IP->fi_csi[i][j] + IP->V[i][j]*IP->fi_eta[i][j])),(1.0/(gamma-1.0)));
                
                HP->rho_hi[i][j] = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(IP->U_hi[i][j]*IP->fi_csi_hi[i][j] + IP->V_hi[i][j]*IP->fi_eta_hi[i][j])),(1.0/(gamma-1.0)));
                HP->rho_hj[i][j] = pow((1.0 - ((gamma-1.0)/(gamma+1.0))*(IP->U_hj[i][j]*IP->fi_csi_hj[i][j] + IP->V_hj[i][j]*IP->fi_eta_hj[i][j])),(1.0/(gamma-1.0)));
                
            }
        }
    }
    void HPrhoArtf(HalfPoints * HP, IntegerPoints * IP)
    {
        int i=0, j=0, r=0 ,s=0;
        double nu = 0.0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                if (HP->U_hi[i][j] < 0)
                {
                    r = 1;
                    HP->nu_hi[i][j] = ((IP->M[i+1][j]*IP->M[i+1][j]-1.0)*C > 0.0) ? (IP->M[i+1][j]*IP->M[i+1][j]-1.0)*C : 0.0;
                }
                else
                {
                    r = -1;
                    HP->nu_hi[i][j] = ((IP->M[i][j]*IP->M[i][j]-1.0)*C > 0.0) ? (IP->M[i][j]*IP->M[i][j]-1.0)*C : 0.0;
                }
                if (HP->V_hj[i][j] < 0)
                {
                    s = 1;
                    HP->nu_hj[i][j] = ((IP->M[i][j+1]*IP->M[i][j+1]-1.0)*C > 0.0) ? (IP->M[i][j+1]*IP->M[i][j+1]-1.0)*C : 0.0;
                }
                else
                {
                    s = -1;
                    HP->nu_hj[i][j] = ((IP->M[i][j]*IP->M[i][j]-1.0)*C > 0.0) ? (IP->M[i][j]*IP->M[i][j]-1.0)*C : 0.0;
                }
                
                HP->rho_til_hi[i][j] = (1.0-HP->nu_hi[i][j])*HP->rho_hi[i][j]) + HP->nu_hi[i][j]*HP->rho_hi[i+r][j];
                HP->rho_bar_hj[i][j] = (1.0-HP->nu_hj[i][j])*HP->rho_hj[i][j]) + HP->nu_hj[i][j]*HP->rho_hj[i][j+s];
                
            }
        }
    }
    void VarSetCal(HalfPoints * HP, IntegerPoints * IP)
    {
        // Calculate Half-Point Fi derivatives - fi_csi and fi_eta
        HPfiDeriv(HP,IP);
        
        // Calculate Half-Point Velocities
        HPvel(HP,IP);
        
        // Calculate Half-Point Densities
        HPrho(HP,IP);
        
        // Calculate Half-Point Artificial Densities
        HPrhoArtf(HP,IP);
    }
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Residue Calculation ------------------- */
// ------------------------------------------------------------ //
    void Residue(HalfPoints * HP, IntegerPoints * IP)
    {
        IP->resMax = 0.0;
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                
                IP->residue[i][j] = ((HP->rho_til_hi[i][j]*HP->U_hi[i][j]/HP->J_hi[i][j]) - (HP->rho_til_hi[i-1][j]*HP->U_hi[i-1][j]/HP->J_hi[i-1][j])) + ((HP->rho_bar_hj[i][j]*HP->V_hj[i][j]/HP->J_hj[i][j]) - (HP->rho_bar_hj[i][j-1]*HP->V_hj[i][j-1]/HP->J_hj[i][j-1]));
                IP->resMax = IP->resMax + IP->residue[i][j]
            }
        }
        
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Iteration Method ---------------------- */
// ------------------------------------------------------------ //
    
    // Thomas Algorithm
    void Thomas(IntegerPoints * IP)
    //void Thomas_Algorithm(double ** C , double * a, double * b, double * c, double * RHS, int k, int dir)
    {
    	int i,j;
    	
    	if(dir==1)	// line on x
    	{
    		for(i=1;i<IMAX-1;i++)
    		{	
    			if(i==1)
    			{
    				c[i] = (double) c[i]/b[i];
    				RHS[i] = (double) RHS[i]/b[i];
    			}
    			else
    			{
    				c[i] = (double) (c[i]/(b[i]-a[i]*c[i-1]));
    				RHS[i] = (double)((RHS[i]-a[i]*RHS[i-1])/(b[i]-a[i]*c[i-1]));
    			}
    		}
    		for(i=IMAX-2;i>0;i--)
    		{	
    			if(i==IMAX-2)
    			{
    				C[i][k] = RHS[i];
    			}
    			else
    			{
    				C[i][k] = RHS[i]-c[i]*C[i+1][k];
    			}
    		}
    	}
    	else if(dir==2) // line on y
    	{
    		for(j=1;j<JMAX-1;j++)	
    		{
    			if(j==1)
    			{
    				c[j] = (double)(c[j]/b[j]);
    				RHS[j] = (double)(RHS[j]/b[j]);
    			}
    			else
    			{
    				c[j] = (double)(c[j]/(b[j]-a[j]*c[j-1]));
    				RHS[j] = (double)((RHS[j]-a[j]*RHS[j-1])/(b[j]-a[j]*c[j-1]));
    			}
    		}
    		for(j=JMAX-2;j>0;j--)
    		{	
    			if(j==JMAX-2)
    			{
    				C[k][j] = RHS[j];
    			}
    			else
    			{
    				C[k][j] = RHS[j]-c[j]*C[k][j+1];
    			}
    		}
    	}
    }
    
    // AF2
    void IterMethod(IntegerPoints * IP, )
    {
        // Step 1:
            
        // Step 2:
            
    }
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Boundary Conditions ------------------- */
// ------------------------------------------------------------ //
    void BCs(IntegerPoints * IP)
    {
        
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Printing Functions -------------------- */
// ------------------------------------------------------------ //

    // Printing when converges
    void PrintConv(int n, IntegerPoints * IP)
    {
        printf("\nConverged at %d iteration!! Residue = %f\n",n,IP->resMax);
    }
    
    // Printing when diverges
    void PrintDiverg(int n, IntegerPoints * IP)
    {
        printf("\nDiverged at %d iteration!! Residue = %f\n",n,IP->resMax);
    }
    
    // Printing iteration
    void PrintIter(int n, IntegerPoints * IP, double timeIter)
    {
        printf("\nIteration = %d ..... Residue = %f .... Time = %f \
                seconds\n",n,IP->resMax,timeIter);
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------------- Output -------------------------- */
// ------------------------------------------------------------ //
    void Output(IntegerPoints * IP)
    {
        
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //
    