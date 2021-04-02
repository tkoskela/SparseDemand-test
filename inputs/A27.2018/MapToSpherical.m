function [r,phi]=MapToSpherical(x)
% Map vector x to spherical coordinates
n=length(x);
x2=x.*x;
r=sqrt(sum(x2));
phi=zeros(n-1,1);
for i=1:n-2 
  phi(i)=atan2(sqrt(sum(x2(i+1:n))),x(i));
end
phi(n-1)=atan2(x(n),x(n-1));
  