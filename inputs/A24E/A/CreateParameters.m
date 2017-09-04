
% 1) C.raw           (nC x 1)  nC = K*(K-1)/2 + (J-K)*(K-1)   
% 2) D.raw           (J x 1)
% 3) MUE.raw         (K x 1)
% 4) InvCDiag.raw    (K x 1)
% 5) InvCOffDiag.raw (nInvC x 1)  K*(K-1)/2 
% 6A) BC_Z  (nBC x BC_z_dim)
% 6B) BC_beta        (BC_z_dim x 1)
% 7A  BD_Z           (J x BD_z_dim)
% 7B) BD_beta        (BD_z_dim x 1)
% 8) BC_CDiag        (nBC x 1)
% 9) BC_COffDiag     (nBC_C x 1)  nBC_C = BC_eta_dim*(BC_eta_dim-1)/2 + (nBC-BC_eta_dim)*(BC_eta_dim-1)
% 10) BD_CDiag       (J x 1) 
% 11) BD_COffDiag    (nBD_C x 1)  nBD_C = BD_eta_dim*(BD_eta_dim-1)/2 + (J-BD_eta_dim) *(BD_eta_dim-1)

J  = 24;
K  = 5;
nBC = K*(K-1)/2 + (J-K)*(K-1);
BC_z_dim = J-1;
BD_z_dim = J;
BC_eta_dim = 3;
BD_eta_dim = 3;

% 1) C.raw           (nC x 1)  nC = K*(K-1)/2 + (J-K)*(K-1)   
C = 0.5*pi * ones(nBC,1);
fid = fopen('C.raw','w');
fprintf(fid,'%25.16f \n',C);
fclose(fid);

% 2) D.raw           (J x 1)
D = 0.5 + rand(J,1);
fid = fopen('D.raw','w');
fprintf(fid,'%25.16f \n',D);
fclose(fid);

% 3) MUE.raw         (K x 1)
MUE = 20*ones(K,1);
fid = fopen('MUE.raw','w');
fprintf(fid,'%25.16f \n',MUE);
fclose(fid);

% 4) InvCDiag.raw    (K x 1)
% 5) InvCOffDiag.raw (nInvC x 1)  K*(K-1)/2 
sig = eye(K);
InvS = inv(sig);
InvC = chol(InvS,'lower');
[InvCDiag,InvCOffDiag] = MatrixToSphere(InvC);

fid = fopen('INVCDiag.raw','w');
fprintf(fid,'%25.16f \n',InvCDiag);
fclose(fid);
fid = fopen('INVCOffDiag.raw','w');
fprintf(fid,'%25.16f \n',InvCOffDiag);
fclose(fid);

% 6A) BC_Z  (nBC x BC_z_dim)
BC_Z = zeros(nBC,BC_z_dim);
for i1=2:J
  if i1<=K
    index = (i1-1)*(i1-2)/2 + (1:i1-1)';
  else 
    index = K*(K-1)/2 + (K-1)*(i1-K-1) + (1:K-1)';
  end
  BC_Z(index,i1) = 1;
end
fid = fopen('BC_Z.raw','w');
fprintf(fid,'%25.16f \n',BC_Z);
fclose(fid);

% 6B) BC_beta           (BC_z_dim x 1)
BC_beta = randn(BC_z_dim,1);
fid = fopen('BC_beta.raw','w');
fprintf(fid,'%25.16f \n',BC_beta);
fclose(fid);

% 7A) BD_Z        (J x BD_z_dim)
BD_Z = eye(J);
fid = fopen('BD_Z.raw','w');
fprintf(fid,'%25.16f \n',BD_Z);
fclose(fid);

% 7) BD_beta        (BD_z_dim x 1)
BD_beta = sort(0.25*rand(BD_z_dim,1)-0.25,1,'descend');
fid = fopen('BD_beta.raw','w');
fprintf(fid,'%25.16f \n',BD_beta);
fclose(fid);

% 8) BC_CDiag          (nBC x 1)
% 9) BC_COffDiag       (nBC_C x 1)  nBC_C = BC_eta_dim*(BC_eta_dim-1)/2 + (nBC-BC_eta_dim)*(BC_eta_dim-1)
BC_CDiag = 0.1*ones(nBC,1);
fid = fopen('BC_CDiag.raw','w');
fprintf(fid,'%25.16f \n',BC_CDiag);
fclose(fid);

nBC_C = BC_eta_dim*(BC_eta_dim-1)/2 + (nBC-BC_eta_dim)*(BC_eta_dim-1);
BC_COffDiag = pi*rand(nBC_C,1);
fid = fopen('BC_COffDiag.raw','w');
fprintf(fid,'%25.16f \n',BC_COffDiag);
fclose(fid);

% 10) BD_CDiag         (J x 1) 
BD_CDiag = 0.1*ones(J,1);
fid = fopen('BD_CDiag.raw','w');
fprintf(fid,'%25.16f \n',BD_CDiag);
fclose(fid);

% 11) BD_COffDiag      (nBD_C x 1)  nBD_C = BD_eta_dim*(BD_eta_dim-1)/2 + (J-BD_eta_dim) *(BD_eta_dim-1)
nBD_COffDiag = BD_eta_dim*(BD_eta_dim-1)/2 + (BD_eta_dim-1)*(J-BD_eta_dim);
BD_COffDiag = pi*rand(nBD_COffDiag,1);
fid = fopen('BD_COffDiag.raw','w');
fprintf(fid,'%25.16f \n',BD_COffDiag);
fclose(fid);

% create sigP
SigP = 2*rand(J,J);
SigP = SigP'*SigP;
fid = fopen('sigp.raw','w');
fprintf(fid,'%25.16f \n',SigP);
fclose(fid);


