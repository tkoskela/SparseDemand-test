function CreateElasTable(filename,elas,q_labels)

J = size(elas,1);

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


for j1=1:6

disp('');
disp('\begin{table}[h]');
disp(['\caption{Elasticities (',int2str(j1),')}']);
disp(['\label{Table: elasticities ',int2str(j1),'}']);
disp('\begin{center}');
disp('\begin{tabular}{ccccc}');

disp(['Price & ',q_labels{4*(j1-1)+1},' & ', ...
                 q_labels{4*(j1-1)+2},' & ', ...
                 q_labels{4*(j1-1)+3},' & ', ...
                 q_labels{4*(j1-1)+4},' \\ \hline']);
for i1=1:J
  disp(['p$_{',q_labels{i1},'}$ & ', ...
        num2str(elas(i1,4*(j1-1)+1)),' & ', ...
        num2str(elas(i1,4*(j1-1)+2)),' & ', ...
        num2str(elas(i1,4*(j1-1)+3)),' & ', ...
        num2str(elas(i1,4*(j1-1)+4)), ...
        ' \\ ']);
end
disp('\end{tabular}');
disp('\end{center}');
disp('\end{table}');
disp('');
disp('');
disp('');

end % for j1=1:6

disp('\end{document}');  
diary off;