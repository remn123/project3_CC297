#ifndef FUNCTIONS_H
#define FUNCTIONS_H

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
    void AllocateAll(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP, ThomStr * T);
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Initial Condition --------------------- */
// ------------------------------------------------------------ //
    void InitialCond(MeshGrid * mesh, IntegerPoints * IP);
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* --------------------- Read Mesh ---------------------------- */
// ------------------------------------------------------------ //
    void ReadMesh(MeshGrid * mesh, char * filename, IntegerPoints * IP, HalfPoints * HP, ThomStr * T);
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* ------- Constant Initial Settings and Calculations --------- */
// ------------------------------------------------------------ //
    void InverseMesh(MeshGrid * mesh);
    void HPderivPM(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void HPjacobian(HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void HPderivNM(HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void HPmetrics(HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void FreeVars(HalfPoints * HP, IntegerPoints * IP);

    void ConstSetCal(MeshGrid * mesh, HalfPoints * HP, IntegerPoints * IP);
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Boundary Conditions ------------------- */
// ------------------------------------------------------------ //
    void BCs(IntegerPoints * IP, HalfPoints * HP, MeshGrid * mesh);
/* ---------- Variables Settings and Calculations ------------- */
// ------------------------------------------------------------ //
    
    void HPfiDeriv(HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void HPvel(HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void HPrho(HalfPoints * HP, IntegerPoints * IP, int i, int j);
    void HPrhoArtf(HalfPoints * HP, IntegerPoints * IP);
    void VarSetCal(HalfPoints * HP, IntegerPoints * IP, MeshGrid * mesh);
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Residue Calculation ------------------- */
// ------------------------------------------------------------ //
    void Residue(HalfPoints * HP, IntegerPoints * IP, int n);
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Iteration Method ---------------------- */
// ------------------------------------------------------------ //
    
    // Thomas Algorithm
    void Thomas_Algorithm(double ** C , ThomStr * T, int k, int dir, int I, int J);
    
    // AF2
    void IterMethod(IntegerPoints * IP, HalfPoints * HP, ThomStr * T, int iter);
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //


// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------- Printing Functions -------------------- */
// ------------------------------------------------------------ //

    // Printing when converges
    void PrintConv(int n, IntegerPoints * IP);

    // Printing when diverges
    void PrintDiverg(int n, IntegerPoints * IP);
    
    // Printing iteration
    void PrintIter(int n, IntegerPoints * IP, double timeIter);
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------------- Free-up ------------------------- */
// ------------------------------------------------------------ //
    
    void FreeUp(IntegerPoints * IP, HalfPoints * HP, MeshGrid * mesh, ThomStr * T);

// ------------------------------------------------------------ //
// ------------------------------------------------------------ //

/* -------------------------- Output -------------------------- */
// ------------------------------------------------------------ //
    void Output(MeshGrid * mesh, IntegerPoints * IP);
    
// ------------------------------------------------------------ //
// ------------------------------------------------------------ //
#endif