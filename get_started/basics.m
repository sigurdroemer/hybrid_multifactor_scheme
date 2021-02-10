%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(project_folder));

%% Description:
% In this script we illustrate the approximation of (completely monotone) kernel functions
% and show how to simulate using the hybrid multifactor scheme (by example on the rough  
% Bergomi model).

%% Set-up:
H = 0.1;
T = 1;
n = 500;
delta = 1/n;
t = (delta:delta:T)';
K = @(t)(t.^(H-0.5));


%% Approximate kernel function:
[c,gamm,K,K_hat] = ApproximateKernel('BM2005','K',K,'epsilon',0.01,'n',n,...
                                     'delta',delta,'T',T,'approx_roots',true,'Nz',10^4);

% Inspect solution:
[c,gamm]

% Show approximation:
figure;
plot(t,K(t),'DisplayName','Truth','color','blue','LineWidth',1.5);hold on;
plot(t,K_hat(t),'--','DisplayName','Approximation','color','red','LineWidth',1.5);hold on;
xlim([delta-0.05,T+0.05]);
xlabel('t');ylabel('Value');
legend();

%% Simulate:
N = 10;
x0 = 0;
kappa = 1;
b = 0;
sigma = 1;

% Simulate Gaussian Volterra process:
Xpaths = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K);

% Inspect solution:
Xpaths

% Plot:
figure;
plot(Xpaths.t,Xpaths.values);
xlabel('t');
ylabel('Values');

% Construct and plot rough Bergomi paths:
xi = 0.15.^2;
eta = 2.1;
V = xi.*exp(eta*sqrt(2*H)*Xpaths.values - 0.5.*Xpaths.t.^(2*H)*eta^2);

figure;
plot(Xpaths.t,sqrt(V));
xlabel('t');
ylabel('Volatility');


%% Compare pricing against the hybrid scheme for TBSS processes:
% Additional inputs:
s0 = 100;
y = 0;
q = 0;
rho = -0.9;
strikes = (60:1:130)';
ttm = 1*ones(size(strikes));
N = 50000;
n = 1000;
kappa = 1;

% Pricing:
tic;
[iv_1,se_1] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,'HybridMultifactor',N,n,...
                                    kappa,'c',c,'gamm',gamm,'turbo',true);
toc;
tic;
[iv_2,se_2] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,'HybridTBSS',N,n,...
                                    kappa,'turbo',true);
toc;

% Plot:
figure;
plot(log(strikes./s0),iv_1,'LineWidth',1.5,'Color','blue');hold on;
plot(log(strikes./s0),iv_2,'--','LineWidth',1.5,'Color','red');
xlabel('Log-moneyness');
ylabel('Implied volatility');
legend('Hybrid multifactor scheme','Hybrid-TBSS scheme');

figure;
plot(log(strikes./s0),iv_1,'LineWidth',1.0,'Color','blue');hold on;
plot(log(strikes./s0),iv_1-se_1*1.96,'--','LineWidth',0.5,'Color','red');
plot(log(strikes./s0),iv_1+se_1*1.96,'--','LineWidth',0.5,'Color','red');
legend('Hybrid multifactor scheme','Confidence interval (95%)');
xlabel('Log-moneyness');
ylabel('Implied volatility');




