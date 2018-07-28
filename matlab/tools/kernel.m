function f = kernel(x,xi,bw)

n   = size(x,1);
nf  = size(xi,1);
sig = sqrt(var(x));
if nargin<3
  bw =  1.06 * sig * n ^(1/5);
end

% sum( K( (x - xi)/bw) ) / (n*h)

z = (repmat(x',nf,1) - repmat(xi,1,n))/bw;

f = sum(exp(-0.5*z.*z)/sqrt(2*pi),2)/(n*bw);
