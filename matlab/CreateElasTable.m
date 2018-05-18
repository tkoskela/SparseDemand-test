function CreateElasTable(filename,elas,q_labels)
system(['rm ',filename]);
J     = size(elas,1);
ncol0 = 7;
prec  = 2;  % precision

diary(filename);

disp('\documentclass[11pt]{article}');
disp('\usepackage{geometry}');
disp('\geometry{a4paper}');
disp('\usepackage{graphicx}');
disp('\usepackage{amssymb}');
disp('\usepackage{epstopdf}');
disp('\DeclareGraphicsRule{.tif}{png}{.png}{`convert #1 `dirname #1`/`basename #1 .tif`.png}');
disp('');
disp('\title{Elasticities}');
disp('\author{Lars Nesheim}');
disp('\date{}');
disp('');
disp('\begin{document}');
disp('\maketitle');


for j1=1:ceil(J/ncol0)

c1 = ncol0*(j1-1)+1;
c2 = min(ncol0*j1,J);
ncol = (c2-c1)+1;

disp('');
disp('\begin{table}[h]');
disp(['\caption{Elasticities (',int2str(j1),')}']);
disp(['\label{Table: elasticities ',int2str(j1),'}']);
disp('\begin{center}');

% Number of columns = ncol+1
str1 = '\begin{tabular}{l';
for i1=1:ncol
  str1 = strcat(str1,'c');
end
str1 = strcat(str1,'} \hline \hline');
disp(str1);

str1 = 'Price';
for i1=1:ncol
    str1 = strcat(str1,' & ',q_labels{ncol*(j1-1)+i1});
end
str1 = strcat(str1,' \\ \hline');
disp(str1);

for i1=1:J
  str1 = ['p$_{',q_labels{i1},'}$ & '];
  for i2=1:ncol
    str1=strcat(str1,' & ',num2str(elas(i1,ncol*(j1-1)+i2),prec));
  end
  str1 = strcat(str1,' \\');
  disp(str1);
end
disp('\end{tabular}');
disp('\end{center}');
disp('\end{table}');
disp('');
disp('');
disp('');

end % for j1=1:

disp('\end{document}');  
diary off;