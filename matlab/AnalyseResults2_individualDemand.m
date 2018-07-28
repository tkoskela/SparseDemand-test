
file1 = fullfile(FortranOutDir,'IndividualDemand.txt');
IndividualDemand = ImportIndividualDemand(file1);

% Column indexes of variables
cols.ihh      = 1;  % HH ID
cols.jp       = 2;  % index of price that varies
cols.ix       = 3;           
cols.eta      = (4:5)';      % columns with draw of random coefficient
cols.n        = 6;           % number of goods purchased
cols.iNonZero = (7:11)';     % indexes of goods with nonzero demand
cols.q        = (12:16)';    % columns with quantities
cols.p        = (17:16+24)'; % columns with prices

nhh=9;  % number of distinct households in data

% Set all demands to 0 if n==0 (q might be >0 but less than crit)  
i0 = IndividualDemand(:,cols.n)==0;
if any(i0)
  IndividualDemand(i0,cols.iNonZero) = 0;
end
% Set elements of quantity to zero that should be zero
%  (some values are >0 but less than crit. these must be set to zero)
for i1=1:K-1
  i0 = IndividualDemand(:,cols.n)==i1;
  if any(i0)
    IndividualDemand(i0,cols.iNonZero(i1+1:K))=0;
  end
end

FullQ         = zeros(np*J,J,nhh);
UniqueNonZero = cell(nhh,1);
iMaxQ         = zeros(nhh,1);

fig=0;

close all

eta = zeros(nhh,2);
nprices = zeros(nhh,1);
for ihh=1:nhh
  % Extract demand for current HH
  iTemp = IndividualDemand(:,cols.ihh)==ihh;
  demand = IndividualDemand(iTemp,:);
  eta(ihh,:) = demand(1,cols.eta);
  
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
  for j1=1:nq
    % plot demands vs p(j1)
    jp = UniqueNonZero{ihh}(j1);  % index of current price
    i0 = demand(:,cols.jp)==jp;
    iNonZero = demand(i0,cols.iNonZero);
    iNonZero = iNonZero(:);
    SubsetUnique = unique(iNonZero(iNonZero~=0));
    nq1 = length(SubsetUnique);
    
    fig=fig+1;
    figure(fig);
    FigFile = fullfile(OutDir,['hh',int2str(ihh),'_p',int2str(jp),'.eps']);
    
    for j2=1:nq1
      subplot(ceil(nq1/2),2,j2)
      jq = SubsetUnique(j2);
      plot(demand(i0,cols.p(jp)),FullQ(i0,jq,ihh))
      xlabel([FruitLabels{jp},' price (GBP)']);
      ylabel(FruitLabels{jq});
      %xlabel(['p(',int2str(jp),')'])
      %ylabel(['q(',int2str(jq),')'])
      print(FigFile,'-depsc2')
    end      
  end  
end 

% Plot density of nitems for each household type
fig =fig+1;
figure(fig)
FigFile = fullfile(OutDir,'nitems.eps');
for ihh=1:nhh
  subplot(3,3,ihh)
  itemp = IndividualDemand(:,cols.ihh)==ihh;
  nitems = IndividualDemand(itemp,cols.n);
  histogram(nitems)
end
print(FigFile,'-depsc2');

% Create table listing names of products for each consumer
CreateTableProductNames(eta,UniqueNonZero,FruitLabels, ...
                        fullfile(OutDir,'baskets.tex'));
                    
close all