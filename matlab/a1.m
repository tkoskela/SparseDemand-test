cols.ihh      = 1;
cols.jp       = 2;  % index of price that varies
cols.ix       = 3;
cols.eta      = (4:5)';
cols.n        = 6;
cols.iNonZero = (7:11)';
cols.q        = (12:16)';
cols.p        = (17:16+24)';

J = 24;
K = 5;
nhh=9;
np=30;
OutDir = '/home/uctpln0/FruitDemand/code/fortran/output/2017SEP/A24D';

i0 = IndividualDemand(:,cols.n)==0;
if any(i0)
  IndividualDemand(i0,cols.iNonZero) = 0;
end
for i1=1:K-1
  i0 = IndividualDemand(:,cols.n)==i1;
  if any(i0)
    IndividualDemand(i0,cols.iNonZero(i1+1:K))=0;
  end
end

FullQ         = zeros(np*J,J,nhh);
UniqueNonZero = cell(nhh,1);
iMaxQ         = zeros(nhh,1);
nprices=2;
fig=0;

close all

for ihh=1:nhh
  % Extract demand for current HH
  iTemp = IndividualDemand(:,1)==ihh;
  demand = IndividualDemand(iTemp,:);
  
  % Find unique product id's of products purchased by ihh
  iNonZero           = demand(:,cols.iNonZero);
  iNonZero           = iNonZero(:);
  UniqueNonZero{ihh} = unique(iNonZero(iNonZero~=0));
  
  % convert sparse demand matrix to full demand matrix
  for j1=1:J
  for k1=1:K
    iq = demand(:,cols.iNonZero(k1))==j1;
    if any(iq)
      FullQ(iq,j1,ihh) = demand(iq,cols.q(k1));
    end
  end 
  end
  
  nq = length(UniqueNonZero{ihh});
  for j1=1:nprices
    % plot demands vs p(j1)
    jp = UniqueNonZero{ihh}(j1);  % index of current price

    fig=fig+1;
    figure(fig);
    FigFile = fullfile(OutDir,['hh',int2str(ihh),'_p',int2str(jp),'.eps']);
    
    for j2=1:nq
      subplot(ceil(nq/2),2,j2)
      jq = UniqueNonZero{ihh}(j2);
      i0 = demand(:,cols.jp)==jp;
      plot(demand(i0,cols.p(jp)),FullQ(i0,jq))
      xlabel(['p(',int2str(jp),')'])
      ylabel(['q(',int2str(jq),')'])
      print(FigFile,'-depsc2')
    end      
  end  
end        
