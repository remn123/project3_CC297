#ifndef STRUCTS_H
#define STRUCTS_H
    
    /* ------- STRUCT 1 --------
    Half-Point Struct (HPstruct):
         
         (1) - x derivatives:
            - x_csi:
               x_csi = d_csi(x)
            
            - x_eta:
               x_eta = d_eta(x)
         
         (2) - y_derivatives:
            - y_csi:
               y_csi = d_csi(y)
         
            - y_eta:
               y_eta = d_eta(y)
               
         (3) - metrics:
            - A1:
               A1 = csi_x*csi_x + csi_y*csi_y
            - A2:
               A2 = csi_x*eta_x + csi_y*eta_y
            - A3:
               A3 = eta_x*eta_x + eta_y*eta_y
         
         (4) - Jacobian:
            - J:
               J = 1.0/(x_csi*y_eta - x_eta*y_csi);
         
         (5) - Fi Derivatives :
            - fi_csi:
               fi_csi = 
            
            - fi_eta:
               fi_eta = 
            
         (6) - Velocities
            - U:
               U = A1*fi_csi + A2*fi_eta
            - V:
               V = A2*fi_csi + A3*fi_eta
         
         (7) - Densities
            - rho_hp:
               rho_hp = 
            
            - nu:
               nu =
            
            - rho_tilda:
               rho_tilda = 
               
            - rho_hat:
               rho_hat =
         
    ---------------------------- */
    typedef struct HalfPoints
    {
       // Metric Parameters
       double ** xhi_csi;
       double ** yhi_csi;
       double ** xhj_eta;
       double ** yhj_eta;
       
       double ** xhj_csi;
       double ** yhj_csi;       
       double ** xhi_eta;
       double ** yhi_eta;

       double ** csi_xhi;
       double ** csi_yhi;
       double ** eta_xhi;
       double ** eta_yhi;
       
       double ** csi_xhj;
       double ** csi_yhj;
       double ** eta_xhj;
       double ** eta_yhj;
       
       double ** A1_hi; // half i
       double ** A1_hj; // half j
       
       double ** A2_hi; // half i
       double ** A2_hj; // half j
       
       double ** A3_hi; // half i
       double ** A3_hj; // half j
       
       // Jacobian
       double ** J_hi;  // half i
       double ** J_hj;  // half j
       
       // Fi derivatives
       double ** fi_csi_hi;   // half i
       double ** fi_csi_hj;   // half j
       
       double ** fi_eta_hi;   // half i
       double ** fi_eta_hj;   // half j
       
       // Velocities
       double ** U_hi;  // half i
       double ** U_hj;  // half j
       double ** V_hi;  // half i
       double ** V_hj;  // half j
       
       // Half-point Density
       double ** rho_hi;   // half i
       double ** rho_hj;   // half j
       
       // Density limiter
       double ** nu_hi; // half i
       double ** nu_hj; // half j
       
       // Densities
       double ** rho_til_hi;   // half i;
       double ** rho_bar_hj;   // half j;
       
       
    } * HP;
    /* ------- STRUCT 2 --------
    Integer-Point Struct (IPstruct):
         
         (1) - residue
         
         (2) - fi
         
         (3) - correction
            
         (4) - density
            
    ---------------------------- */
    typedef struct IntegerPoints
    {
       double ** residue;
       double ** fi;
       double ** corr;
       double ** rho;
       double resMax;
       
       // Derivatives
       double ** x_csi;
       double ** x_eta;
       double ** y_csi;
       double ** y_eta;
       
       double ** csi_x;
       double ** csi_y;
       double ** eta_x;
       double ** eta_y;
       
       double ** J;
       
       // Metric Parameters
       double ** A1;
       double ** A2;
       double ** A3;
       
       // Fi derivatives
       double ** fi_csi;
       double ** fi_eta;
       
       // Velocities
       double ** U;
       double ** V;
       
       // Mach
       double M_INF;
       
    } * IP;
    
    // Mesh struct
    typedef struct MeshGrid
    {
       double ** x;
       double ** y;
       
    } * mesh;
    
    int IMAX;
    int JMAX;
    
    // Thomas Algorithm input
    typedef struct ThomStr
    {
       double * a;
       double * b;
       double * c;
       double * RHS;
       double * LHS;
       int k;
       
    } * T;
#endif