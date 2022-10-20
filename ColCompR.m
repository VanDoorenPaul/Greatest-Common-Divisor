function [Qright,Mr,r] = ColCompR(M,tol)
%
%   function [Qright,Mr,r]=ColCompR(M,tol)
%
%  compresses the columns of M to r linearly independent  
%  columns to the right side of the matrix Mr = M*Qright.
%  Here r=rank(M,tol) with tol a chosen tolerance.
%
mn = size(M); m = mn(1); n = mn(2);
[U,S,V] = svd(M); 
minmn=min([m,n]); s = diag(S(1:minmn,1:minmn));
r = sum(s > tol);
if r == 0, Qright = eye(n,n); Mr = zeros(size(M)); return, end
if r == n, Qright = eye(n,n); Mr = M; return, end
Qright = [V(:,r+1:n) V(:,1:r)];
Mr = [zeros(m,n-r) U*S(:,1:r)];
end

