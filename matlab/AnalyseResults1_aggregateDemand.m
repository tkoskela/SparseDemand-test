% 1) load demand data from current results
% 2) plot predicted demand vs actual
% 3) plot demand curve for each category
% 4) Create table of elasticities
addpath('ImportTools')

%FortranOutDir = '/SAN/economics/Nesheim-IO/FruitDemand/output/A27_20171116';
%FortranOutDir = '/SAN/economics/Nesheim-IO/FruitDemand/output/A27_20180101';
FortranOutDir = '../output/A27_2018MAY/fortran_output';
OutDir        = '../output/A27_2018MAY/figures';
%FortranOutDir = '/SAN/economics/Nesheim-IO/FruitDemand/output/A27_2018MAY';
%OutDir  = '/home/uctpln0/FruitDemand/code/fortran/output/A27_2018MAY';
N       = 26514;
J       = 27;
K       = 5;
np      = 30;

if exist(FortranOutDir,'dir')~=7
  display(['FortranOutDir = ',FortranOutDir]);    
  error('Directory does not exist.')
end
if exist(OutDir,'dir')~=7
  e=mkdir(OutDir); 
  if (e~=1) 
    error(['Failed to create OutDir = ',OutDir]);
  end
end

[qdata,qhat,qaverage] = ...
    ImportDemandData(fullfile(FortranOutDir,'demand_data.csv'), 1, J);

fig1=1;
FontSize = 16;

% Plot qdata vs qhat
% Plot qhat vs qaverage
[~,i1]=sort(qdata);
figure(fig1)
subplot(2,1,1)
hold off
plot(log(qdata(i1)+1),log(qhat(i1)+1),'rx')
hold on
plot(log(qdata(i1)+1),log(qdata(i1)+1),'k')
xlabel('Actual log consumption','FontSize',FontSize);
ylabel('Predicted log demand','FontSize',FontSize);

[~,i1]=sort(qhat);
subplot(2,1,2)
hold off
plot(log(qhat(i1)+1),log(qaverage(i1)+1),'rx')
hold on
plot(log(qhat(i1)+1),log(qhat(i1)+1),'k')
xlabel('Predicted log demand: actual prices','FontSize',FontSize);
ylabel('Predicted log demand: average prices','FontSize',FontSize);
print(fig1,fullfile(OutDir,'fig1_model_fit.eps'),'-depsc');

% for each product, plot demand
price = zeros(np,J);
quantity = zeros(np,K+1,J);
qtemp    = zeros(np,K+1);
FruitLabels = {'Apricots'; ...
               'Avocados'; ...
               'Bananas'; ...
               'Berries'; ...
               'Cherries'; ...
               'Dates'; ...
               'Apples'; ...
               'Easy Peelers'; ...
               'Grapes'; ...
               'Grapefruits'; ...
               'Kiwis'; ...
               'Lemons'; ...
               'Limes';       ...
               'Lychees'; ...
               'Mangos'; ...
               'Melons'; ...
               'Nectarines'; ...
               'Oranges'; ...
               'Passion fruits'; ...
               'Paw-paws'; ...
               'Peaches';  ...
               'Pears'; ...
               'Pineapples'; ...
               'Plums'; ...
               'Pomegranates'; ...
               'Rhubarb'; ...
               'Sharon fruits'};
for j1=1:J
  if (j1<10)
    str1 = ['demand0',int2str(j1),'.csv'];
  else
    str1 = ['demand',int2str(j1),'.csv'];
  end
  str1 = fullfile(FortranOutDir,str1);
  str2 = sprintf('%c%c%c%c',FruitLabels{j1});
  [p,q1,q2,q3,q4,q5,q6]=ImportDemandData2(str1, 1, np);
  [price(:,j1),      ...
   quantity(:,1,j1), ...
   quantity(:,2,j1), ...
   quantity(:,3,j1), ...
   quantity(:,4,j1), ...
   quantity(:,5,j1), ...
   quantity(:,6,j1)] = ImportDemandData2(str1, 1, np);  
   fig1=fig1+1;
   figure(fig1)
   
   NewPlotFlag=1;
   if NewPlotFlag==0   
     subplot(2,1,1)
     plot(price(:,j1),quantity(:,1,j1)/N)
     title(['Demand for ',FruitLabels{j1}],'FontSize',FontSize);
     xlabel('price (GBP)','FontSize',FontSize)
     ylabel('quantity (kilograms)','FontSize',FontSize)
     subplot(2,1,2)
     plot(price(:,j1),reshape(quantity(:,2:K+1,j1),[np K]));
     legend('One item','Two items','Three items','Four items','Five items')
     print(fig1, ...
         fullfile(OutDir,['fig',int2str(j1),'_',str2,'.eps']), ...
         '-depsc');
   elseif NewPlotFlag==1
     baseprice = (2/3) * min(price(:,j1)) + (1/3)*max(price(:,j1));  
     area(price(:,j1),quantity(:,2:K+1,j1))
     title(['Demand for ',FruitLabels{j1}],'FontSize',FontSize);
     xlabel(['price (GBP: average market price = ', ...
             num2str(baseprice,2),' GBP)'],'FontSize',FontSize)
     ylabel('quantity (kilograms)','FontSize',FontSize)
     legend('One item','Two items','Three items','Four items','Five items')
     print(fig1, ...
         fullfile(OutDir,['fig',int2str(j1),'_',str2,'.eps']), ...
         '-depsc');
   end
end

% Load elasticities
elas = ImportElasticity(fullfile(FortranOutDir,'elas.csv'),J, 1, J);
%elas = ImportElasticity2(fullfile(FortranOutDir,'elas.csv'), 1, J);

CreateElasTable(fullfile(OutDir,'elas1.tex'),elas,FruitLabels);
close all;
