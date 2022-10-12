# GCRDr
 
README file for GCRDr (Greatest Common Right Divisor)

This directory contains Matlab files to compute the
 
Compact Greatest Common Right Divisor of a polynomial matrix
 
    P(s) = Q(s).G(s)

where P, Q and G are polynomial matrices

P(s) is mxn (m.ge.n) of normal rank r (r.le.n)
N(s) is mxr and is left invertible over the ring of polynomials
G(s) is rxn and has the Smith zeros and right minimal indices of P(s)

The polynomial matrices are stored in Matlab as 3D arrays.

All files are Matlab function files except for the GCRDTestDriver 
which is a Matlab script that produces test results for the codes.

Each Matlab function file is commented and has a built-in help.