% Matlab file to create regressor index

%RegressorIndex = (9:15)';
% Full vector of regressors: (9:214)'
RegressorIndex = (9:214)';
fid=fopen('inputs/RegressorIndex.txt','w');
n=size(RegressorIndex,1);
for i1=1:n
  ecode = fprintf(fid,'%4d \n',RegressorIndex(i1));
end
fclose(fid);
