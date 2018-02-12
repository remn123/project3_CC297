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
    
    void HPderivPM(MeshGrid * mesh, IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                IP->x_csi[i][j] = d_csi(mesh->x,i,j); 
                IP->y_csi[i][j] = d_csi(mesh->y,i,j);
                IP->x_eta[i][j] = d_eta(mesh->x,i,j); 
                IP->y_eta[i][j] = d_eta(mesh->y,i,j); 
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
            }
        }
        
        // Average in csi
        for(i=1;i<IMAX-2;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                HP->J_hi[i][j] = a_csi(IP->J,i,j);      // J(i+1/2,j)
            }
        }
        // Average in eta
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-2;j++)
            {
                HP->J_hj[i][j] = a_eta(IP->J,i,j);      // J(i,j+1/2)
            }
        }
        
    }
    
    void HPderivNM(IntegerPoints * IP)
    {
        int i=0, j=0;
        
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                IP->csi_x[i][j] = IP->J[i][j]*IP->y_eta[i][j];
                IP->csi_y[i][j] = -IP->J[i][j]*IP->x_eta[i][j];
                IP->eta_x[i][j] = -IP->J[i][j]*IP->y_csi[i][j];
                IP->eta_y[i][j] = IP->J[i][j]*IP->x_csi[i][j];
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
            }
        }
        
        // Half-Point Metric Parameters
        // Average in csi
        for(i=1;i<IMAX-2;i++)
        {
            for(j=1;j<JMAX-1;j++)
            {
                HP->A1_hi[i][j] = a_csi(IP->A1,i,j);    // A1(i+1/2,j)
                HP->A2_hi[i][j] = a_csi(IP->A2,i,j);    // A2(i+1/2,j)
                HP->A3_hi[i][j] = a_csi(IP->A3,i,j);    // A3(i+1/2,j)
            }
        }
        // Average in eta
        for(i=1;i<IMAX-1;i++)
        {
            for(j=1;j<JMAX-2;j++)
            {
                HP->A1_hj[i][j] = a_eta(IP->A1,i,j);    // A1(i,j+1/2)
                HP->A2_hj[i][j] = a_eta(IP->A2,i,j);    // A2(i,j+1/2)
                HP->A3_hj[i][j] = a_eta(IP->A3,i,j);    // A3(i,j+1/2)
            }
        }
        
    }
    
    void ConstSetCal(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP)
    {
        // Setting mesh
        InverseMesh(mesh);
        
        // Calculate Integer-Point derivatives x_csi, x_eta, y_csi and y_eta
        HPderivPM(mesh,IP); // Integer-Point derivatives on Physical Mesh
        
        // Calculate Half-Point and Integer-Point Jacobians
        HPjacobian(HP,IP);
        
        // Calculate Integer-Point derivatives csi_x, csi_y, eta_x and eta_y
        HPderivNM(IP); // Integer-Point derivatives on Numerical Mesh
        
        // Calculate Half-Point and Integer-Point Metric Parameters
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
        
    }
    void HPvel(HalfPoints * HP)
    {
        
    }
    void HPrho(HalfPoints * HP)
    {
        
    }
    
    void VarSetCal(HalfPoints * HP, IntegerPoints * IP)
    {
        // Calculate Half-Point Fi derivatives - fi_csi and fi_eta
        HPfiDeriv(HP,IP);
        
        // Calculate Half-Point Velocities
        HPvel(HP);
        
        // Calculate Half-Point Densities
        HPrho(HP);
    }
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Residue Calculation ------------------- */
// ------------------------------------------------------------ //
    void Residue(HalfPoints * HP, IntegerPoints * IP)
    {
        
    }
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Iteration Method ---------------------- */
// ------------------------------------------------------------ //
    
    // Thomas Algorithm
    void Thomas(IntegerPoints * IP)
    {
        
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
    