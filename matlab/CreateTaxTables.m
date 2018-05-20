function CreateTaxTables(filename,productlabels,tau,N,TaxLabel, ...
                         q0,qtax,p0,ptax, ...
                         e0,etax,u0,utax)
% Create tables of results from tax simulations
% (q0,qtax,p0,ptax)   quantity and price
% (e0,etax,u0,utax)   expenditure and welfare
% tau = vector of percentiles at which quantiles are computed
% N   = population
% TaxLabel = labels for each tax scenario

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
disp('');
disp('\begin{document}');

%% Table 1: results on quantities
disp('');
disp('\begin{table}[h]');
disp('\caption{Percentage change in demand due to tax/price change}');
disp('\label{table:tax impact 1}');
disp('\begin{center}');
% Number of columns = 2 + ntax
disp('\resizebox{1 \textwidth}{!}{')
str1 = '\begin{tabular}{lc';
for t1=1:ntax
  str1 = strcat(str1,'c');
end
str1 = [str1,'} \hline \hline'];
disp(str1);

% column labels
str1 = ' & ';
str2 = 'Fruit & Baseline (kg)';
for i1=1:ntax
  str1 = [str1, ' & Scenario ',int2str(i1)];
  str2 = strcat(str2,' & ',TaxLabel{i1});
end
str1 = [str1,' \\ '];
str2 = [str2,' \\ \hline'];
disp(str1);
disp(str2);

irow = 0;
for j1=1:J
    
  str1 = strcat(productlabels{j1},' & ',num2str(q0(j1)/N,prec));
  for i1=1:ntax
    % percent change in demand
    x1 = 100*(qtax(j1,i1)/q0(j1)-1)  ;
    str1 = strcat(str1,' & ',num2str(x1,prec),'\%');
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
str1 = ['\multicolumn{',int2str(2+ntax),'}{p{1.0 \textwidth}}{',...
        'Note: The first column shows baseline demand for each fruit ', ...
        '(kilograms per household per shopping trip). The remaining ',...
        'columns show the percentage change in demand resulting from the ',...
        'change in tax or prices.}'];
disp(str1)
disp('\end{tabular}}');
disp('\end{center}');
disp('\end{table}');
disp('');

%% Table 2: impacts on prices
disp('\begin{table}[h]');
disp('\caption{Percentage change in price due to tax/price change}');
disp('\label{table:tax impact 2}');
disp('\begin{center}');
disp('\resizebox{1 \textwidth}{!}{');
% Number of columns = 2 + ntax
str1 = '\begin{tabular}{lc';
for t1=1:ntax
  str1 = strcat(str1,'c');
end
str1 = [str1,'} \hline \hline'];
disp(str1);

% column labels
str1 = ' & ';
str2 = 'Fruit & Baseline';
for i1=1:ntax
  % use cell array with strcat to keep trailing white space  
  str1 = [str1,' & Scenario ',int2str(i1)];
  str2 = strcat(str2,' & ',TaxLabel{i1});
end
str1 = [str1,' \\ '];
str2 = [str2,' \\ \hline'];
disp(str1);
disp(str2);

irow = 0;
for j1=1:J
  str1 = strcat(productlabels{j1},' & ',num2str(p0(j1),prec));
  for i1=1:ntax
    x1 = 100*(ptax(j1,i1)/p0(j1)-1);
    str1 = strcat(str1,' & ',num2str(x1,prec),'\%');
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
str1 = ['\multicolumn{',int2str(2+ntax),'}{p{1.0 \textwidth}}{',...
        'Note: The first column shows the baseline price for each fruit ', ...
        '(GBP per kilogram). The remaining ',...
        'columns show the percentage impact of the change in tax or prices.}'];
disp(str1);     
disp('\end{tabular}}');
disp('\end{center}');
disp('\end{table}');
disp('');

%% Table 3:  Expenditure and welfare
disp('\begin{table}[h]');
disp('\caption{Tax impact on expenditure and welfare}');
disp('\label{table:tax impact welfare}');
disp('\begin{center}');
disp('\resizebox{1 \textwidth}{!}{');
% Number of columns = 2 + ntax
str1 = '\begin{tabular}{lc';
for t1=1:ntax
  str1 = strcat(str1,'c');
end
str1 = [str1,'} \hline \hline'];
disp(str1);

% column labels
str1 = ' & ';
str2 = ' & Baseline';
for i1=1:ntax
  str1 = [str1,' & Scenario ',int2str(i1)];
  str2 = strcat(str1,' & ',TaxLabel{i1});
end
str1 = [str1,' \\'];
str2 = [str2,' \\ \hline'];
disp(str1);
disp(str2);

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
    x1 = 100*(etax(i1,i2)/e0(i1)-1);  
    str1 = strcat(str1,' & ',num2str(x1,prec),'\%');
  end
  str1 = strcat(str1,' \\');
  disp(str1);
end
disp('\hline ')

str1 = 'Change in consumer surplus (GBP) & ';
for i1=1:ntax
  str1 = strcat(str1,' & ');
end
str1 = strcat(str1,' \\');
disp(str1);

ntau = length(tau);
for i1=1:ntau
  str1 = [num2str(100*tau(i1)),'th percentile &',num2str(u0(i1),prec) ];
  for i2=1:ntax
    str1 = strcat(str1,' & ',num2str(utax(i1,i2)-u0(i1),prec));
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

str1 = ['Consumer surplus (GBP) & ',num2str(u0(ntau+1),prec)];
for i1=1:ntax
  str1 = strcat(str1,' & ',num2str(utax(ntau+1,i1)-u0(ntau+1),prec));    
end
str1 = strcat(str1,' \\');
disp(str1)

str1 = 'Tax revenue (GBP) & 0.0 ';
for i1=1:ntax
  str1 = strcat(str1,' & ',num2str(((ptax(:,i1)-p0).'*qtax(:,i1))/N,prec));
end
str1 = strcat(str1,' \\');
disp(str1);

str1 = ['Firm Revenue & ',num2str((p0'*q0)/N,prec)];
for i1=1:ntax
  x1 = 100*( (p0'*qtax(:,i1))/(p0'*q0)-1);
  str1 = strcat(str1,' & ',num2str(x1,prec),'\%');    
end
str1 = strcat(str1,' \\');
disp(str1);

disp(' \hline \hline');
str1 = ['\multicolumn{',int2str(2+ntax),'}{p{1.0 \textwidth}}{',...
        'Note: The first column shows the baseline values for ',...
        'expenditure, consumer surplus, firm revenue and ',...
        'tax revenue. All amounts are measured in pounds per ',...
        'household per shopping trip. Columns 2 - 7 show the ', ...
        'percentage change in expenditure, the absolute change in ',...
        'consumer surplus, the percentage change in firm revenue ', ...
        'and the absolute change in tax revenue arising in each scenario. ', ...
        'Because of quasilinear utility ',...
        'the change in consumer surplus equals compensating variation.}']; 
disp(str1);
disp('\end{tabular}}');
disp('\end{center}');
disp('\end{table}');
disp('');

disp('\end{document}');  
diary off;