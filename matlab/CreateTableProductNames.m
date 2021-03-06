function CreateTableProductNames(eta,UniqueNonZero,FruitLabels,filename)
system(['rm ',filename]);
nhh=size(eta,1);

diary(filename);

disp('\documentclass[11pt]{article}');
disp('\usepackage{geometry}');
disp('\geometry{a4paper}');
disp('\usepackage{graphicx}');
disp('\usepackage{amssymb}');
disp('\usepackage{epstopdf}');
disp('');
disp('\begin{document}');

disp('');
disp('\begin{table}[h]');
disp('\caption{Fruit baskets: by household type}');
disp('\label{table:fruit baskets}');
disp('\begin{center}');
disp('\tiny');
disp('\begin{tabular}{ll}');
disp('\hline \hline');
disp('HH ID & Fruit basket \\');
for ihh=1:nhh
  basket = FruitLabels(UniqueNonZero{ihh});
  nb = length(UniqueNonZero{ihh});
  for j1=1:nb
      if j1==1
        str1 = basket{j1};
      elseif (j1>1)
        str1 = strcat(str1,', ',basket{j1});
      end
  disp([int2str(ihh),' & ',str1,' \\']);
  end
end
disp('\hline \hline');
disp('\end{tabular}');
disp('\end{center}');
disp('\end{table}');
disp('');
disp('');
disp('');

disp('\end{document}');  
diary off;