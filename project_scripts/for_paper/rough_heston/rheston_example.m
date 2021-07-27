%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));


%% Define model:
model = VolterraCEVModelClass('xi',0.15^2);
model.eta = 0.29;
model.beta = 0.5;
H = 0.12;alpha = H - 0.5;
model.K = KernelFunctionClass(1,alpha,{@(obj,t)(1./gamma(1+obj.alpha))},10^(-12));
model.rho = -0.67;

s0 = 100;
model.s0 = s0;

model.pricerSettings.precision = 'single';
model.pricerSettings.N = 10000;
model.pricerSettings.n(:) = 500;
model.pricerSettings.numRepeat = 100;
model.pricerSettings.simulation.kernel_approximation.error_measure = 'l2_normalised';
model.pricerSettings.price_estimation.S.control_variate = 'none';
model.pricerSettings.price_estimation.VIX.control_variate = 'none';

% Define contracts:
ttm = 0.1;

S_contracts = struct;
S_contracts.k = (-0.2:0.01:0.05)';
S_contracts.ttm = ttm*ones(size(S_contracts.k));
strikes = model.s0*exp(S_contracts.k);


%% Parameter set A:
% * We use the default parameters.
n_test = [16,32,64,128,256,512,10000]*10;

% Compute "truth":
[price, iv_truth_par1] = NumericalIntegrationRoughHeston(s0,model.xi.values,model.K.alpha + 1,0,...
                                                model.xi.values,model.eta,model.rho,...
                                                true,strikes,ttm,'N',1500,...
                                                'disp_iter',true);

figure;
cm = parula(size(n_test,2));
for i=1:size(n_test,2)
    i
    model.pricerSettings.n(:) = n_test(i);
    pS = model.GetPrices(S_contracts,[]);
    plot(pS.k,pS.p,'-','color',cm(i,:),'LineWidth',1.25);hold on;
    drawnow;
end

plot(log(strikes./model.s0),iv_truth_par1,'--','LineWidth',1.25,'color','red');hold on;

xlim([-0.22,0.075]);
ylim([0.085,0.3]);
xlabel('Log-moneyness','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Rough Heston','(case A)'},'interpreter','latex');
legend('\lfloornT\rfloor = 16','\lfloornT\rfloor = 32','\lfloornT\rfloor = 64',...
       '\lfloornT\rfloor = 128','\lfloornT\rfloor = 256','\lfloornT\rfloor = 512',...
       '\lfloornT\rfloor = 10000','Fourier pricing');


%% Parameter set B:
model.K.alpha = 0.05 - 0.5;
model.eta = 0.41;
model.rho = -0.67;

% Compute "truth":
[price, iv_truth_par2] = NumericalIntegrationRoughHeston(s0,model.xi.values,model.K.alpha + 1,0,...
                                                model.xi.values,model.eta,model.rho,...
                                                true,strikes,ttm,'N',500,...
                                                'disp_iter',true);

n_test = [16,32,64,128,256,512,10000]*10;


figure;
cm = parula(size(n_test,2));
for i=1:size(n_test,2)
    i
    model.pricerSettings.n(:) = n_test(i);
    pS = model.GetPrices(S_contracts,[]);
    plot(pS.k,pS.p,'-','color',cm(i,:),'LineWidth',1.25);hold on;
    drawnow;
end

plot(log(strikes./model.s0),iv_truth_par2,'--','LineWidth',1.25,'color','red');hold on;


xlim([-0.22,0.075]);
ylim([0.05,0.375]);
xlabel('Log-moneyness','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Rough Heston','(case B)'},'interpreter','latex');
legend('\lfloornT\rfloor = 16','\lfloornT\rfloor = 32','\lfloornT\rfloor = 64',...
       '\lfloornT\rfloor = 128','\lfloornT\rfloor = 256','\lfloornT\rfloor = 512',...
       '\lfloornT\rfloor = 10000','Fourier pricing');


%% Compute mean value term structure:
% To argue where the bias comes from:
n = 500;
N = 10000;
tt = (1/n:1/n:1)';
paths = model.Simulate(N,n,{'V'},{tt'});

figure;
X = paths.V.values;
mu = mean(X)';
se = std(X)'./sqrt(N);
plot(tt,[mu-2*se,mu,mu+2*se]);

hold on;
plot(tt,ones(size(tt))*model.xi.values,'--','LineWidth',1.5,'color','red')
ylim([0.02,0.03])


%% Check convergence of Fourier solution (case A).
% Case B is checked in the relevant rBergomi script.

N_test = [1000,1500,2000];
[price,iv] = deal(NaN(size(strikes,1),size(N_test,2)));
for i=1:size(N_test,2)
[price(:,i), iv(:,i)] = NumericalIntegrationRoughHeston(s0,model.xi.values,model.K.alpha + 1,0,...
                                                model.xi.values,model.eta,model.rho,...
                                                true,strikes,ttm,'N',N_test(i),...
                                                'disp_iter',true);
end
                                            

%% VIX (case A)
model.eta = 0.29;
model.beta = 0.5;
H = 0.12;alpha = H - 0.5;
model.K = KernelFunctionClass(1,alpha,{@(obj,t)(1./gamma(1+obj.alpha))},10^(-12));
model.rho = -0.67;

rng(1234);
n_test = 10*[128,256,512,1024,2048,4096]';

V_contracts = struct;
%V_contracts.q = (0.01:0.01:0.99)';
V_contracts.K = (0.07:0.005:0.4)';
V_contracts.ttm = 0.1;

model.pricerSettings.simulation.kernel_approximation.n_cap = 20000;
model.pricerSettings.N = 10000;
model.pricerSettings.numRepeat = 100;

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

% Compute with 10,000 steps
model.pricerSettings.n(:) = 100000;
model.pricerSettings.N = 1000;
model.pricerSettings.numRepeat = 1000;
[~,pV_ref] = model.GetPrices([],V_contracts);  
plot(pV_ref.K,pV_ref.p,'--','color','red','LineWidth',1.25); 
legend('\lfloornT\rfloor = 128','\lfloornT\rfloor = 256',...
       '\lfloornT\rfloor = 512','\lfloornT\rfloor = 1024','\lfloornT\rfloor = 2048',...
       '\lfloornT\rfloor = 4096','\lfloornT\rfloor = 10000');
xlabel('Strike','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Rough Heston, case A','(options on VIX)'},'interpreter','latex');
%xlim([0.03,0.42]);
%ylim([0.4,1.8]);



%% VIX (case B)
model.eta = 0.41;
model.beta = 0.5;
H = 0.05;alpha = H - 0.5;
model.K = KernelFunctionClass(1,alpha,{@(obj,t)(1./gamma(1+obj.alpha))},10^(-12));
model.rho = -0.67;

rng(1234);
n_test = 10*[128,256,512,1024,2048,4096]';

V_contracts = struct;
%V_contracts.q = (0.01:0.01:0.99)';
V_contracts.K = (0.05:0.005:0.4)';
V_contracts.ttm = 0.1;

model.pricerSettings.simulation.kernel_approximation.n_cap = 20000;
model.pricerSettings.N = 10000;
model.pricerSettings.numRepeat = 100;

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

% Compute with 10,000 steps:
model.pricerSettings.n(:) = 100000;
model.pricerSettings.N = 1000;
model.pricerSettings.numRepeat = 1000;
[~,pV_ref] = model.GetPrices([],V_contracts);  
plot(pV_ref.K,pV_ref.p,'--','color','red','LineWidth',1.25); 
legend('\lfloornT\rfloor = 128','\lfloornT\rfloor = 256',...
       '\lfloornT\rfloor = 512','\lfloornT\rfloor = 1024','\lfloornT\rfloor = 2048',...
       '\lfloornT\rfloor = 4096','\lfloornT\rfloor = 10000');
xlabel('Strike','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Rough Heston, case B','(options on VIX)'},'interpreter','latex');
xlim([0.03,0.42]);
%ylim([0.4,1.8]);





