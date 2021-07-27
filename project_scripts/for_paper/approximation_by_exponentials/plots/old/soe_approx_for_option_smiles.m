%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%% Compute sum-of-exponentials approximation with 10,000 steps:
T = 0.1;
n = 10000; % Will have to sample n + 1 points in [delta,T] 
delta = T/n; % (check ApproximateKernel function to see this).

H = 0.05;
[c_H05,gamm_H05,K,K_hat,exitflag] = ApproximateKernel('BM2005','H',H,'T',T,'delta',delta,'n',n);           

H = 0.12;
[c_H12,gamm_H12,K,K_hat,exitflag] = ApproximateKernel('BM2005','H',H,'T',T,'delta',delta,'n',n);           

save([fileparts(matlab.desktop.editor.getActiveFilename),'\K_approximations_10000_steps.mat'],...
      'c_H05','gamm_H05','c_H12','gamm_H12');
  
%%