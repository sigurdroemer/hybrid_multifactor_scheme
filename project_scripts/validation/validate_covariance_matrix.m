%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%%
H = 0.1;
alpha = H - 0.5;
kappa = 10;
n = 50;

K = KernelFunctionClass(1,alpha,@(obj,t)(1),10^(-12));

SIGMA = GetVolterraCovarianceMatrix({K},kappa,1/n);
covM = CovMatrixHybrid(n,kappa,alpha,'double');

max(max(abs(SIGMA-covM)./SIGMA))

figure;
surf(SIGMA(3:end,3:end))


