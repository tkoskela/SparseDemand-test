J = 4;
K = 3;

D = 0.5+rand(J,1);
save('D.raw','D','-ascii','-double');

nphi = K*(K-1)/2+(J-K)*(K-1);
phi = (0.40+0.2*rand(nphi,1))*pi;
save('PHI.raw','phi','-ascii','-double');

mue = 0.5*ones(K,1);
save('MUE.raw','mue','-ascii','-double');

InvC_D = ones(K,1);
save('INVC_D.raw','InvC_D','-ascii','-double');

InvC_phi = 0.5*pi*ones(K*(K-1)/2,1);
save('INVC_PHI.raw','InvC_phi','-ascii','-double');

% (mu_d,sig_d,gamma_d)
% (mu_c,sig_c,gamma_c)



