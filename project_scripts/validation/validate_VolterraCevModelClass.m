%% Initialize Script
clear;
project_folder = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(project_folder));

%% Setup
model = VolterraCEVModelClass('xi',0.15^2,'eta',0.1,'beta',0.5,'rho',-0.67);
model.K = KernelFunctionClass(1,-0.38,{@(obj,t)(1./gamma(1+obj.alpha))},10^(-12));

model.pricerSettings.precision = 'single';
model.pricerSettings.N = 10000;
model.pricerSettings.n(:) = 500;
model.pricerSettings.numRepeat = 100;
model.pricerSettings.simulation.kernel_approximation.error_measure = 'l2_normalised';
model.pricerSettings.price_estimation.S.control_variate = 'none';
model.pricerSettings.price_estimation.VIX.control_variate = 'none';

%% Validation:
N = 10000;
n = 2000;
ttV = (0:1/n:1);
ttVIX = (0:0.1:1);
paths = model.Simulate(N,n,{'V','VIX2','S'},{ttV,ttVIX,ttV},'antithetic',false);

figure;
mu = mean(paths.V.values)';
se = std(paths.V.values)'./sqrt(N);
plot(paths.V.t,[mu-1.96*se,mu,mu+1.96*se],'-');hold on;

figure;
mu = mean(paths.VIX2.values)';
se = std(paths.VIX2.values)'./sqrt(N);
plot(paths.VIX2.t,[mu-1.96*se,mu,mu+1.96*se],'-');hold on;

figure;plot(paths.S.t,mean(paths.S.values));
ylim([0.99,1.01])

% Note that E(V) is significantly biased when model is very rough and/or 
% vol-vol is high. This is as expected for the particular scheme we are using.
