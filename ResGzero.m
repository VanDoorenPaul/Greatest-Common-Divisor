function [zeroval] = ResGzero(G,sm)
%
% [zeroval] = ResGzero(G,smithzero) 
%
%  computes the smallest singular value of a polynomial G(s)
%  at each of its presumed Smith zeros stored in sm
mnd=size(G);d=mnd(3)-1;r=min(mnd(1:2));k=length(sm);
for i=1:k, 
    M=G(:,:,1);
    for j=1:d
        M=M+G(:,:,j+1)*sm(i)^j;
    end
    s=svd(M);zeroval(i)=s(r)/s(1);
end