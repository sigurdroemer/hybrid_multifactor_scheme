%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));


%% Define model:
% rHeston parameters:
H = 0.05;
eta = 0.41;
rho = -0.67;
xi = 0.15^2;

model = QuadraticRoughHestonClass('zeta',0.05^2);
model.a = 0.384;
model.b = 0.095;
model.c = 0.0025;
model.eta = 1;
model.rho = 1;

alpha = H - 0.5;
z0 = 0.1;

model.K = KernelFunctionClass(1,alpha,{@(obj,t)(1./gamma(1+obj.alpha))},10^(-12));

s0 = 100;
model.s0 = s0;

model.pricerSettings.precision = 'single';
model.pricerSettings.N = 10000;
model.pricerSettings.n(:) = 500;
model.pricerSettings.numRepeat = 1;
model.pricerSettings.seed = [];
model.pricerSettings.price_estimation.S.control_variate = 'none';
model.pricerSettings.simulation.kernel_approximation.error_measure = 'l2_normalised';

% Load rHeston prices:
load(['C:\Users\sigur\GoogleDrive\university\phd\research\',...
      'rough_volatility_source_code\05_GitHub\hybrid_multifactor_scheme\',...
      'project_scripts\for_paper\rbergomi\rHeston_prices.mat']);

% Define contracts:
strikes = (80:2.5:110)';
ttm = 0.1;  
  
S_contracts = struct;
S_contracts.k = log(strikes./model.s0);
S_contracts.ttm = 0.1*ones(size(strikes));
S_contracts.mid = iv_rH(:,end);
S_contracts.bid = S_contracts.mid;
S_contracts.ask = S_contracts.mid;


%% Calibrate model:
diffChg = 0.001;
maxIter = 500;
options = optimoptions('lsqnonlin',...
                       'Display','iter','DiffMinChange',diffChg,'MaxIter',maxIter);                 
                   

parNames = {'c','a','eta','b','zeta'};  
par0 = {model.c;model.a;model.eta;model.b;model.zeta.values}';
lb = {0;0;0;-10;-10}';
ub = {1;30;30;10;10}';

f = @(obj,s)(1./gamma(obj.alpha+1));                         
fTmp = @(model,alpha)(KernelFunctionClass(1,alpha,f,10^(-12)));   

% Calibrate to SPX options:
model.pricerSettings.n(:) = 5000;
model.pricerSettings.N = 10000;
model.pricerSettings.price_estimation.antithetic = false;
rng(123);
res = model.Calibrate(S_contracts,[],parNames,par0,lb,ub,...
                                              'options',options,...
                                              'dispIter',true,...
                                              'plotFit',true,...
                                              'error_measure','dist_to_mid'); 

% Save model:
save([fileparts(matlab.desktop.editor.getActiveFilename),'/qrh_calibrated.mat'],'model');
% a = 0.32, b = 0.52, c = 0.003, eta = 0.95, rho = 1, alpha = -0.45, zeta = 0.42                                          
load([fileparts(matlab.desktop.editor.getActiveFilename),'/qrh_calibrated.mat']);


%% Compute convergence (S-options):
rng(123);
n_test = 10*[16,32,64,128,256,512]';

% Set parameters:
model.a = 0.32;
model.b = 0.52;
model.c = 0.003;
model.eta = 0.95;
model.rho = 1;
model.zeta.values = 0.42;
model.K.alpha = -0.45;

S_contracts = struct;
S_contracts.k = (-0.2:0.01:0.05)';
S_contracts.ttm = 0.1*ones(size(S_contracts.k));

strikes = model.s0*exp(S_contracts.k);

model.pricerSettings.N = 10000;
model.pricerSettings.numRepeat = 100;
model.pricerSettings.seed = [];

figure;

% Compute small n's:
cm = parula(size(n_test,1));
iv_mat_hmf = NaN(size(strikes,1),size(n_test,1));
for i=1:size(n_test,1)
    i
    model.pricerSettings.n(:) = n_test(i);
    pS = model.GetPrices(S_contracts,[]);
    iv_mat_hmf(:,i) = pS.p;
    plot(pS.k,pS.p,'-','color',cm(i,:),'LineWidth',1.25);hold on;
    drawnow;
end

% Compute with 10,000 steps:
% Load kernel approximation:
model.pricerSettings.n(:) = 100000;
model.pricerSettings.N = 1000;
model.pricerSettings.numRepeat = 1000;
pS_ref = model.GetPrices(S_contracts,[]);  
plot(pS_ref.k,pS_ref.p,'--','color','red','LineWidth',1.25); 
legend('\lfloornT\rfloor = 16','\lfloornT\rfloor = 32','\lfloornT\rfloor = 64',...
       '\lfloornT\rfloor = 128','\lfloornT\rfloor = 256','\lfloornT\rfloor = 512',...
       '\lfloornT\rfloor = 10000');
xlabel('Log-moneyness','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Quadratic rough Heston','(options on underlying asset)'},'interpreter','latex');
xlim([-0.225,0.075]);
ylim([0.07,0.35]);


%% Compute convergence (VIX-options):
rng(1234);
n_test = 10*[16,32,64,128,256,512]';

V_contracts = struct;
V_contracts.K = (0.12:0.005:0.35)';
V_contracts.ttm = 0.1*ones(size(V_contracts.K));

model.pricerSettings.simulation.kernel_approximation.n_cap = 20000;
model.pricerSettings.N = 10000;
model.pricerSettings.numRepeat = 100;
model.pricerSettings.n_VIX2_integrate(:) = 250;

figure;

% Compute small n's:
cm = parula(size(n_test,1));
iv_mat_hmf = NaN(size(V_contracts.K,1),size(n_test,1));
for i=1:size(n_test,1)
    i
    model.pricerSettings.n(:) = n_test(i);
    [~,pV] = model.GetPrices([],V_contracts);
    iv_mat_hmf(:,i) = pV.p;
    plot(pV.K,pV.p,'-','color',cm(i,:),'LineWidth',1.25);hold on;
    drawnow;
end

% Compute with 10,000 steps and 2500 integration points:
% Load kernel approximation:
model.pricerSettings.n_VIX2_integrate(:) = 2500; 
model.pricerSettings.n(:) = 100000;
model.pricerSettings.N = 1000;
model.pricerSettings.numRepeat = 1000;
[~,pV_ref] = model.GetPrices([],V_contracts);  
plot(pV_ref.K,pV_ref.p,'--','color','red','LineWidth',1.25); 
legend('\lfloornT\rfloor = 16','\lfloornT\rfloor = 32','\lfloornT\rfloor = 64',...
       '\lfloornT\rfloor = 128','\lfloornT\rfloor = 256','\lfloornT\rfloor = 512',...
       '\lfloornT\rfloor = 10000');
xlabel('Strike','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Quadratic rough Heston','(VIX options)'},'interpreter','latex');
%xlim([0.1,0.375]);
%ylim([0.4,1.8]);


%% Convergence in number of integration points:

nI_test = [10,25,50,100,250,2500]';

model.pricerSettings.n(:) = 5000;

figure;
% Compute small n's:
cm = parula(size(nI_test,1));
iv_mat_hmf = NaN(size(V_contracts.K,1),size(nI_test,1));
for i=1:size(nI_test,1)
    i
    rng(123);
    model.pricerSettings.n_VIX2_integrate(:) = nI_test(i);
    [~,pV] = model.GetPrices([],V_contracts);
    iv_mat_hmf(:,i) = pV.p;
    plot(pV.K,pV.p,'-','color',cm(i,:),'LineWidth',1.1);hold on;
    drawnow;
end

xlim([0.1,0.33]);
ylim([0.5,2.2]);

% Conclusion: 250 steps seems sufficient.





