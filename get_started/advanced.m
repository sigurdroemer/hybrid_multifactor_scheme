%% Initialize Script
clear;
project_folder = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(project_folder));

%% Description:
% In this script we illustrate some more advanced features. See the function descriptions for 
% more details and options. We here focus on the rough Heston model defined as
%
%   V(t) = V(0) + int_0^t K(t-s) lambda*(theta - V(s))ds + int_0^t K(t-s) eta*sqrt(V(s)) dW(s)
%
%   where V(0),lambda,theta,eta > 0 and K(t) = t^(H-0.5)/gamma(H+0.5) for H in (0,1/2).
%


%% Set-up:
% General:
T = 1;
n = 1000;
delta = 1/n;
t = (0:1:n)'./n;

% Model:
v0 = 0.15^2;
lambda = 2;
theta = 0.15^2;
eta = 0.1;
H = 0.25;
K = @(t)(t.^(H - 0.5));

% Approximate kernel:
[c,gamm] = ApproximateKernel('BM2005','K',K,'epsilon',0.01,'n',n,...
                                     'delta',delta,'T',T,'approx_roots',true,'Nz',10^4);

N = 1;
kappa = 1;

%% Advanced settings:
% Use single precision for speed:
precision = 'single'; 

% Pre-simulate the underlying Gaussians (speeds up repeated evaluation):
M = floor2(n*T) + 1; % (using floor2 avoids round off errors)
Z = randn((M-1)*N,kappa+1,precision);

% Enforce positivity of the solution:
positive = true;

% Request U-factors as well:
returnU = true;

% Pre-compute (various) numerical integrals (can speed up repeated evaluation).
% Note: Can be replaced by analytical expressions where known.
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
SIGMA = GetVolterraCovarianceMatrix({K},kappa,1/n);
w = SIGMA(2:end,1);

% Request only a subset of the time points. Important remark: If the time 
% points do not fit into the simulation grid the points will be truncated
% down to the nearest grid point:
tX = t';
tU = t'; 
% Try e.g. changing the above to tU = (0.1:0.1:T)';

% Define drift and diffusion functions:
b = @(t,x)(lambda.*(theta-x));
sigma = @(t,x)(eta.*sqrt(x));


%% Simulate model:
rng(123)
[Xpaths,Upaths] = HybridMultifactorScheme(N,n,[],v0,b,sigma,gamm,c,kappa,'K',K,...
                                 'precision',precision,'positive',positive,...
                                 'returnU',true,'tX',tX,'tU',tU,'Z',Z,...
                                 'w',w,'SIGMA',SIGMA);

%close all;
figure;
plot(Xpaths.t,sqrt(Xpaths.values),'-');
title('Volatility');
xlabel('t');
ylabel('Values');

figure;
plot(Upaths.t,squeeze(Upaths.values),'-');
title('U-factors');
xlabel('t');
ylabel('Values');
