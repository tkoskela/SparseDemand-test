iFreeBeta2 = (1:206)';
fid=fopen('inputs/iFreeBeta2.txt','w');
n=size(iFreeBeta2,1);
for i1=1:n
  ecode = fprintf(fid,'%-4d \n',iFreeBeta2(i1));
end
fclose(fid);
