%% Initialize Script
clear;
serverRun = false;
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
t = (delta:delta:T)';

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
kappa = 4;

%% Advanced settings:
% Use single precision for speed:
precision = 'single'; 

% Pre-simulate the underlying Gaussians (speeds up repeated evaluation):
M = floor2(n*T)+1; % ~ floor(n*T) + 1 (avoids round off errors)
Z = randn((M-1)*N,kappa+1,precision);

% Enforce positivity of the solution:
positive = true;

% Request U-factors as well:
returnU = true;

% Pre-compute (various) numerical integrals (speeds up repeated evaluation).
% Note: Can be replaced by analytical expressions where known.
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
SIGMA = GetVolterraCovarianceMatrix({K},kappa,1/n);
w = SIGMA(2:end,1);

% Request only a subset of the time points. Important remark: If the time 
% points do not fit into the simulation grid the points will be truncated
% down to the nearest grid point:
tX = t';
tU = t'; 
% E.g. try changing the above to say tU = (0.1:0.1:T)';

% Define drift and diffusion functions:
b = @(t,x)(lambda.*(theta-x));
sigma = @(t,x)(eta.*sqrt(x));


%% Simulate model:
rng(123)
N = 50000;
tX = (0:delta:T-delta);
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,[],v0,b,sigma,gamm,c,kappa,'K',K,...
                                 'precision',precision,'positive',positive,...
                                 'returnU',true,'tX',tX,'tU',tU,'Z',[],...
                                 'w',w,'SIGMA',SIGMA,'returndW',true);

% Compute asset prices:
rho = -0.7;
s0 = 100;
logS = NaN(N,1);
V = Xpaths.values;
dWperp = sqrt(delta)*randn(N,size(tX,2));
dW2 = rho.*dW + sqrt(1-rho^2).*dWperp;
dlog_S = sqrt(V).*dW2 - 0.5*V*delta;
S = zeros(N,size(V,2));
S = exp(cumsum(dlog_S,2));

strikes = (80:5:110)';

[iv,se,S,optPrice,seOfPrice,idxCall] = GetOptionPrices(s0,0,0,strikes,T*ones(size(strikes)),...
                                                        N,n,[],[],S(:,end),T);

close all;
figure;
plot(strikes,iv,'--');hold on;


% Compute prices via rHeston implementation:
[price, iv2] = NumericalIntegrationRoughHeston(s0,v0,H+0.5,lambda,...
                                        theta,eta,rho,true,strikes,T,'N',252,...
                                        'disp_iter',true);

plot(strikes,iv2,'-x');




