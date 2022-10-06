function [Morth]=RowOrthCompl(M)
%
%   function [Morth]=RowOrthCompl(M)
%
%  produces an orthonormal matrix Morth whose rows 
%  are orthogonal to the row space of M.
%  M is assumed of full row rank (i.e. m.le.n)
%
mn = size(M); m = mn(1); n = mn(2); 
if m == 0, Morth = eye(n,n); return, end
if m == n, Morth = zeros(0,n); return, end
[Q,R] = qr(M'); Morth = Q(:,m+1:n)'; 

