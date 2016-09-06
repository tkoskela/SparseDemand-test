% plot likelihood data
LikeDataFile = '../../temp/LikeData.txt';
%data=importdata(LikeDataFile);
fid=fopen(LikeDataFile);
data = textscan(fid,'%d %s %d %f %f','Delimiter',',');
fclose(fid);

% i2 = variable name
% i1 = row of variable
% x  = value of varaible
% L  = likelihood
% (i2,i1,x(i1,i2),L(i1))

i2Max = max(data{1});
for i2=1:i2Max
  index = data{1}==i2;
  figure(i2)
  Lmax = 100*min(data{5}(index));
  plot(data{4}(index),min(data{5}(index),Lmax))
  TempLabel = data{2}(index);
  xlabel(TempLabel{1})
  ylabel('Likelihood')
  title(['L vs. ',TempLabel{1}]);
end
