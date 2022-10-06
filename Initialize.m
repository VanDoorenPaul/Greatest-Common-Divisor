function [C3,Q,Z,E,A,r3] = Initialize(P,tol)
%
%  function [C3,Q,Z,E,A,r3] = Initialize(P,tol)
%
%  compresses the columns of the compound matrix[Pd, ... , P1, P0]
%  to the right and constructs the pencil to be put in staircase form
%  This follows the notation of the paper 
%
%  A state-space approah for polynomial GCD extraction
%
%  The tolerance for the renk check is given via tol.
%  P(s) is given as a three dimensional array of size mxnx(d+1)
%
mnd=size(P);m=mnd(1);n=mnd(2);d=mnd(3)-1;
%  Construct C
C=P(:,:,1);
for i=1:d, C=[P(:,:,i+1) C]; end
[Z,C3,r3]=ColCompR(C,tol);
Q=eye(d*n,d*n);
E=Z(n+1:n+d*n,:);
A=Z(1:d*n,:);



