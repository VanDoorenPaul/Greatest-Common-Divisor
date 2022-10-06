% This script generates 10 mxn polynomial matrices of size mxn 
% and with regular part in Smith form of size rxr and degree d 
% The left and right factors M and N of P=M.La.N have degree k
% 
d=2;m=4;n=3;r=2;k=2;imax=3;
Table=zeros(10,7);
for itest=1:10,
La=zeros(m,n,d+1);smithzero=[];
for i=1:r,
    poly=randn(d+1,1);
    La(i,i,:)=poly;
    smithzero=[smithzero;1./roots(poly)];
end
% Now construct random polynomial matrices M and N of
% degree k and form P= M.La.N with local Smith form La.
M=randn(m,m,k+1);N=randn(n,n,k+1);
% Its degree will be 2k+d
P=PxN(PxN(M,La),N);
normP=norm(P(:));
P=P/normP;
% We also set tol 
tol=10000*eps;
P=Trim(P,tol)
% now run the GCRDr algorithm
[N,G]=GCRDr(P,tol)
G=Trim(G,eps);
N=Trim(N,eps);
% residual errors for factorization and roots
Res=Trim(PxN(N,G),tol)-P;
dn=size(N,3)-1;dg=size(G,3)-1;
resFactor=norm(Res(:));
resZero=ResGzero(G,smithzero)
smithzero'
Table(itest,:)=[norm(N(:)),norm(G(:)),resFactor,resZero];
end
format short e
Table
