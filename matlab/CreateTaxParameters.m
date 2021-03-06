addpath('ImportTools')

% import files containing frequencies of country and store
InputDir     = '/home/uctpln0/FruitDemand/code/fortran/inputs/A27/rawparms';
TaxParmsFile = fullfile(InputDir,'taxparms.csv');
countryfile  = fullfile(InputDir,'origin_uk_eu.csv');
storefile    = fullfile(InputDir,'fascia_asda_sain.csv');

% Import data on frequencies of country of origin by fruit
[UK,EU,OTHER,fruit] = ImportCountries(countryfile);
[OTHER1,ASDA,fruit1] = ImportStores(storefile);

J          = 27;
ntax       = 5;
total      = UK+EU+OTHER;
shares     = [UK./total, EU./total, OTHER./total];
storeshare = [ASDA./(ASDA+OTHER1) OTHER1./(ASDA+OTHER1)];

taxid = (1:ntax)';
taxlabel = {'20 percent increase in price of all fruit'; ...
            '10 percent increase in price of European fruit'; ...
            '10 percent increase in price UK fruit'; ...
            '10 percent reduction in price of all fruit';    ...
            '5 percent increase in price Asda and Sainsbury'};
taxtype = {'ad valorem'; ...
           'ad valorem'; ...
           'ad valorem'; ...
           'ad valorem'; ...
           'ad valorem'};
tax = [1.2*ones(J,1), ...
       1 + 0.1*shares(:,2), ...
       1 + 0.1*shares(:,1), ...
       0.9*ones(J,1),       ...
       1 + 0.05*storeshare(:,1)];
   
%Save tax parms file
fid=fopen(TaxParmsFile,'w');
for i1=1:ntax
  fprintf(fid,'%4i',taxid(i1));  
  if (i1<ntax)
    fprintf(fid,',');
  end
end
fprintf(fid,'\n');
for i1=1:ntax
  fprintf(fid,'%s',taxlabel{i1});    
  if (i1<ntax)
    fprintf(fid,',');
  end
end
fprintf(fid,'\n');

for i1=1:ntax
  fprintf(fid,'%s',taxtype{i1});    
  if (i1<ntax)
    fprintf(fid,',');
  end
end
fprintf(fid,'\n');

for j1=1:J
  for i1=1:ntax
    fprintf(fid,'%12.4f',tax(j1,i1));
    if (i1<ntax)
      fprintf(fid,',');
    end
  end
  if j1<J
    fprintf(fid,'\n');
  end
end
fclose(fid);
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
