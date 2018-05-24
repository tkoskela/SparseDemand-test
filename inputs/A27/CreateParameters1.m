J = 27;
K = 5;
OutDir = 'rawparms';
if exist(OutDir,'dir')==0
  system(['mkdir ',OutDir])
end
CreateParameters2(J,K,OutDir);
