function CreateTaxTables(filename,productlabels,tau,N, ...
                         q0,qtax,p0,ptax, ...
                         e0,etax,u0,utax)
% Create tables of results from tax simulations
% (q0,qtax,p0,ptax)   quantity and price
% (e0,etax,u0,utax)   expenditure and welfare
% tau = vector of percentiles at which quantiles are computed
% N   = population

system(['rm ',filename]);
J    = size(q0,1);
ntax = size(qtax,2);
prec = 3;  % precision

diary(filename);

disp('\documentclass[11pt]{article}');
disp('\usepackage{geometry}');
disp('\geometry{a4paper}');
disp('\usepackage{graphicx}');
disp('\usepackage{amssymb}');
disp('\usepackage{epstopdf}');
disp('\DeclareGraphicsRule{.tif}{png}{.png}{`convert #1 `dirname #1`/`basename #1 .tif`.png}');
disp('');
disp('\title{Aggregate tax results}');
disp('\author{Lars Nesheim}');
disp('\date{}');
disp('');
disp('\begin{document}');
disp('\maketitle');

%% Table 1: results on quantities
disp('');
disp('\begin{table}[h]');
disp('\caption{Tax impact on demand}');
disp('\label{table:tax impact 1}');
disp('\begin{center}');
% Number of columns = 2 + ntax
str1 = '\begin{tabular}{lc';
for t1=1:ntax
  str1 = strcat(str1,'c');
end
str1 = [str1,'} \hline \hline'];
disp(str1);

% column labels
str1 = 'Fruit & Baseline';
for i1=1:ntax
  str1 = strcat(str1,' & Case ',int2str(i1));
end
str1 = [str1,' \\'];
disp(str1);

irow = 0;
for j1=1:J
    
  str1 = strcat(productlabels{j1},' & ',num2str(q0(j1),prec));
  for i1=1:ntax
    str1 = strcat(str1,' & ',num2str(qtax(j1,i1),prec));
  end
  str1 = strcat(str1,' \\');
  disp(str1);
  irow=irow+1;
  if irow==5
    disp('\hline')
    irow=0;
  end
end
disp(' \hline \hline')
disp('\end{tabular}');
disp('\end{center}');
disp('\end{table}');
disp('');

%% Table 2: impacts on prices
disp('\begin{table}[h]');
disp('\caption{Tax impact on prices}');
disp('\label{table:tax impact 2}');
disp('\begin{center}');
% Number of columns = 2 + ntax
str1 = '\begin{tabular}{lc';
for t1=1:ntax
  str1 = strcat(str1,'c');
end
str1 = [str1,'} \hline \hline'];
disp(str1);

% column labels
str1 = 'Fruit & Baseline';
for i1=1:ntax
  str1 = strcat(str1,' & Case ',int2str(i1));
end
str1 = [str1,' \\'];
disp(str1);

irow = 0;
for j1=1:J
    
  str1 = strcat(productlabels{j1},' & ',num2str(p0(j1),prec));
  for i1=1:ntax
    str1 = strcat(str1,' & ',num2str(ptax(j1,i1),prec));
  end
  str1 = strcat(str1,' \\');
  disp(str1);
  irow=irow+1;
  if irow==5
    disp('\hline')
    irow=0;
  end
end
disp(' \hline \hline')
disp('\end{tabular}');
disp('\end{center}');
disp('\end{table}');
disp('');

%% Table 3:  Expenditure and welfare
disp('\begin{table}[h]');
disp('\caption{Tax impact on expenditure and welfare}');
disp('\label{table:tax impact welfare}');
disp('\begin{center}');
% Number of columns = 2 + ntax
str1 = '\begin{tabular}{lc';
for t1=1:ntax
  str1 = strcat(str1,'c');
end
str1 = [str1,'} \hline \hline'];
disp(str1);

% column labels
str1 = ' & Baseline';
for i1=1:ntax
  str1 = strcat(str1,' & Case ',int2str(i1));
end
str1 = [str1,' \\'];
disp(str1);

str1 = 'Consumer expenditure & ';
for i1=1:ntax
  str1 = strcat(str1,' & ');
end
str1 = strcat(str1,' \\');
disp(str1);

ntau = length(tau);
for i1=1:ntau
  str1 = [num2str(100*tau(i1),2),'th percentile &',num2str(e0(i1),prec) ];
  for i2=1:ntax
    str1 = strcat(str1,' & ',num2str(etax(i1,i2),prec));
  end
  str1 = strcat(str1,' \\');
  disp(str1);
end
disp('\hline ')

str1 = 'Consumer welfare & ';
for i1=1:ntax
  str1 = strcat(str1,' & ');
end
str1 = strcat(str1,' \\');
disp(str1);

ntau = length(tau);
for i1=1:ntau
  str1 = [num2str(100*tau(i1)),'th percentile &',num2str(u0(i1),prec) ];
  for i2=1:ntax
    str1 = strcat(str1,' & ',num2str(utax(i1,i2),prec));
  end
  str1 = strcat(str1,' \\');
  disp(str1);
end
disp('\hline ')

str1 = 'Per capita effects & ';
for i1=1:ntax
  str1 = strcat(str1,' & ');
end
str1 = strcat(str1,' \\');
disp(str1);

str1 = ['Welfare & ',num2str(u0(ntau+1),prec)];
for i1=1:ntax
  str1 = strcat(str1,' 7 ',num2str(utax(ntau+1,i1),prec));    
end
str1 = strcat(str1,' \\');
disp(str1)

str1 = 'Tax revenue & 0.0 ';
for i1=1:ntax
  str1 = strcat(str1,' & ',num2str(((ptax(:,i1)-p0).'*qtax(:,i1))/N,prec));
end
str1 = strcat(str1,' \\');
disp(str1);

str1 = ['Firm Revenue & ',num2str((p0'*q0)/N,prec)];
for i1=1:ntax
  str1 = strcat(str1,' & ',num2str((p0'*qtax(:,i1))/N,prec));    
end
str1 = strcat(str1,' \\');
disp(str1);

disp(' \hline \hline')
disp('\end{tabular}');
disp('\end{center}');
disp('\end{table}');
disp('');

disp('\end{document}');  
diary off;