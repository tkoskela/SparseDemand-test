FortranOutDir = '../output/A27_2018MAY/fortran_output';
%OutDir = '/SAN/economics/Nesheim-IO/FruitDemand/output/A27_2018MAY';
qdata_file     = fullfile(FortranOutDir,'qdata.csv');
qdata_hat_file = fullfile(FortranOutDir,'qdata_hat.csv');

addpath('tools');
addpath('ImportTools');
q = ImportRawData(qdata_file);
qhat = ImportRawData(qdata_hat_file);

q0    = sum(q==0,1)';
qhat0 = sum(qhat==0,1)';
qplus    = sum(q>0,1)';
qhatPlus = sum(qhat>0,1)';

J = size(q,2);
N = size(q,1);
nf = 20;
p = linspace(0.05,0.95,nf)';
plotflag=2;
for i1=1:J
  figure(i1)
  iq = (q(:,i1)>0);
  nq=sum(iq);
  iqhat = qhat(:,i1)>0;
  nqhat = sum(iqhat);

  if plotflag==1
    qi = quant(q(iq,i1),p);
    fq=kernel(q(iq,i1),qi);
    fqhat=kernel(qhat(iqhat,i1),qi);  
    plot(qi,[fq,fqhat])
  elseif plotflag==2
    hold off  
    histogram(q(iq,i1))
    hold on
    histogram(qhat(iqhat,i1))
  end    
  str1 = ['data: nq = ',int2str(nq)];
  str2 = ['model: nq = ',int2str(nqhat)];
  legend(str1,str2);

end




