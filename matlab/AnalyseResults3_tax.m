addpath('tools')

taxfile1       = fullfile(WorkDir,'taxresults_aggregate.csv');
taxfile2       = fullfile(WorkDir,'taxresults_hh.csv');
taxfile_tables = fullfile(OutDir,'taxresults.tex');

% Import aggregate tax results
% J x ntax
% q0 q1 --- q5,p0,p1 --- p5

taxresultsaggregate = ImportTaxAggregate(taxfile1);
taxresultshh        = ImportTaxHH(taxfile2);

ntax = (size(taxresultsaggregate,2)-2)/2;
J    = size(taxresultsaggregate,1);
N    = size(taxresultshh,1);


% quantities and prices
q0   = taxresultsaggregate(:,1);
qtax = taxresultsaggregate(:,1+(1:ntax)');
p0   = taxresultsaggregate(:,ntax+2);
ptax = taxresultsaggregate(:,(ntax+2)+(1:ntax)');

% expenditure and welfare
e0   = taxresultshh(:,1);
etax = taxresultshh(:,1+(1:ntax)');
u0   = taxresultshh(:,ntax+2);
utax = taxresultshh(:,(ntax+2)+(1:ntax)');

tau   = [0.1 0.25 0.5 0.75 0.9]';
e0_summary = [quant(e0,tau); mean(e0)];
u0_summary = [quant(u0,tau); mean(u0)];

etax_summary = zeros(size(e0_summary,1),ntax);
utax_summary = etax_summary;
for i1=1:ntax
  etax_summary(:,i1) = [quant(etax(:,i1),tau);mean(etax(:,i1))];    
  utax_summary(:,i1) = [quant(utax(:,i1),tau);mean(utax(:,i1))];
end

CreateTaxTables(taxfile_tables,FruitLabels,tau,N, ...
                q0,qtax,p0,ptax,            ...
                e0_summary,etax_summary,    ...
                u0_summary,utax_summary);



