function [D,phi]=MatrixToSphere(C)
% Convert (upper triangular) cholesky decomposition of a matrix into 
% spherical coordinates
% C =   (K x K) cholesky decomposition of a matrix
% D =   (K x 1) radius of each column
% phi = (K*(K-1)/2 x 1) angles for each column

k=size(C,1);
D = zeros(k,1);
phi = zeros(k*(k-1)/2,1);
D(1) = C(1,1);

for i1=2:k
  [D(i1),phi0]=MapToSpherical(C(1:i1,i1));  
  phi((i1-1)*(i1-2)/2+(1:i1-1)') = phi0.';
end