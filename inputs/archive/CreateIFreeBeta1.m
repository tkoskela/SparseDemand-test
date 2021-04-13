iFreeBeta1 = (1:206)';
fid=fopen('inputs/iFreeBeta1.txt','w');
n=size(iFreeBeta1,1);
for i1=1:n
  ecode = fprintf(fid,'%-4d \n',iFreeBeta1(i1));
end
fclose(fid);
