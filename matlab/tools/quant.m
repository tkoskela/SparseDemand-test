function [q,iq]=quant(x,p)
% Compute p'th quantile of x
[x,ix]=sort(x);
n=length(p);
q = zeros(n,1);
iq = zeros(n,1);
for i1=1:n
  p0=p(i1)*length(x);
  p1=floor(p0);
  p2=ceil(p0);
  lam=p0-p1;
  q(i1)=lam*x(p2)+(1-lam)*x(p1);
  iq(i1) = ix(p2);
end