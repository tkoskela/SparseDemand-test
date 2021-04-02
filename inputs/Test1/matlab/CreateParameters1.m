% Create initial parameter values to simulate multiple discrete choice model
% based on Lewbel and Nesheim (2019)
%   https://www.cemmap.ac.uk/publication/sparse-demand-systems-corners-and-complements/?highlight=nesheim
%
%  J      = number of products
%  K      = maximum number purchased
%  OutDir = subdirectory to store parameters 
J = 3;
K = 2;
OutDir = '../rawparms';
if exist(OutDir,'dir')==0
  system(['mkdir ',OutDir]);
end
CreateParameters2(J,K,OutDir);
