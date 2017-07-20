% plot likelihood data
OutDir = '/SAN/economics/Nesheim-IO/FruitDemand/output/A24_20170629';
LikeDataFile = fullfile(OutDir,'LikeData.txt');
ResultsFile  = fullfile(OutDir,'results.txt');
n=24; % 285 = full
%data=importdata(LikeDataFile);
fid=fopen(LikeDataFile);
data = textscan(fid,'%d %s %d %f %f','Delimiter',',');
fclose(fid);
[VarName,x0,x,gradient]=ImportResults(ResultsFile,2,n+1);


% i2 = variable name
% i1 = row of variable
% x  = value of varaible
% L  = likelihood
% (i2,i1,x(i1,i2),L(i1))

i2Max = max(data{1});
for i2=1:i2Max
  figure(double(i2))
  hold off
  index = (data{1}==i2);
  Lmax = 100*min(data{5}(index));
  plot(data{4}(index),min(data{5}(index),Lmax))
  hold on;
  plot(x(i2),min(data{5}(index)),'ro')
  TempLabel = data{2}(index);
  xlabel(TempLabel{1})
  ylabel('Likelihood')
  title(['L vs. ',TempLabel{1}]);
end
