function [Qleft,Mr,r] = RowCompT(M,tol)
%
%   function [Qleft,Mr,r]=RowCompT(M,tol)
%
%  compresses the rows M to r linearly independent rows 
%  at the top of the matrix Mr = Qleft'*M.
%  Here r=rank(M,tol) with tol a chosen tolerance.
%
mn = size(M); m = mn(1); n = mn(2); 
[U,S,V] = svd(M); 
minmn=min([m,n]); s = diag(S(1:minmn,1:minmn));
r = sum(s > tol);
if r == 0, Qleft = eye(m,m); Mr = zeros(size(M)); return, end
if r == m, Qleft = eye(m,m); Mr = M; return, end
Qleft = U;
Mr = [S(1:r,:)*V';zeros(m-r,n)];