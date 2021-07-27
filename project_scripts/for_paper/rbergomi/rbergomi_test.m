%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));

%% Define model:
% rHeston parameters:
H = 0.05;
eta = 0.41;
rho = -0.67;
xi = 0.15^2;

% Define rBergomi model:
model = MixedrBergomiClass('zeta',xi,'theta',1,'c',0);
model.eta_1 = 2.1;
model.eta_2 = 0;
model.rho_12 = -0.8;
model.rho_13 = 0;
model.rho_23 = 0;


alpha = H - 0.5;

model.K_1 = KernelFunctionClass(1,alpha,{@(obj,t)( sqrt(2*obj.alpha+1) )},10^(-12));
model.K_2 = KernelFunctionClass(1,0,{@(obj,t)(1)},10^(-12));

s0 = 100;
model.s0 = s0;

model.pricerSettings.precision = 'single';
model.pricerSettings.N = 10000;
model.pricerSettings.n(:) = 500;
model.pricerSettings.numRepeat = 1;
model.pricerSettings.price_estimation.S.control_variate = 'none';
model.pricerSettings.simulation.kernel_approximation.error_measure = 'l2_normalised';

%% Compute rHeston smile:
% Define contracts:
strikes = (80:2.5:110)';
ttm = 0.1;

N_test = [1500,2000,2500]';

% Compute "truth":
iv_rH = NaN(size(strikes,1),size(N_test,1));
for i=1:size(N_test,1)
    [~, iv_rH(:,i)] = NumericalIntegrationRoughHeston(s0,xi,H + 0.5,0,...
                                               xi,eta,rho,...
                                               true,strikes,ttm,'N',N_test(i),...
                                               'disp_iter',true);
end

% Save result:
save([fileparts(matlab.desktop.editor.getActiveFilename),'/rHeston_prices.mat'],'strikes','iv_rH');
    
% Check convergence:
iv_rH

% Plot:
figure;
plot(log(strikes./model.s0),iv_rH(:,end),'--','DisplayName','Fourier pricing',...
     'LineWidth',1.5,'color','red');hold on;

S_contracts = struct;
S_contracts.k = log(strikes/model.s0);
S_contracts.ttm = 0.1*ones(size(strikes));
S_contracts.mid = iv_rH(:,end);
S_contracts.bid = S_contracts.mid;
S_contracts.ask = S_contracts.mid;


%% Calibrate:
diffChg = 0.001;
maxIter = 500;
options = optimoptions('lsqnonlin',...
                       'Display','off','DiffMinChange',diffChg,'MaxIter',maxIter,...
                       'TolX',10^(-10),'TolFun',10^(-10));                 
                   
parNames = {'eta_1','rho_12','zeta'};  
par0 = {model.eta_1;model.rho_12;model.zeta.values}';
lb = {0;-1;0}';
ub = {30;1;1}';

% Calibrate to SPX options:
model.pricerSettings.n(:) = 5000;
model.pricerSettings.N = 10000;
model.pricerSettings.price_estimation.antithetic = false;
model.pricerSettings.seed = 123;
result = model.Calibrate(S_contracts,[],parNames,par0,lb,ub,...
                          'options',options,...
                          'dispIter',true,...
                          'error_measure','dist_to_mid',...
                          'plotFit',true); 
model.pricerSettings.seed = [];

% Save model:
save([fileparts(matlab.desktop.editor.getActiveFilename),'/rBergomi_calibrated.mat'],'model');
% eta = 2.93, rho = -1.00, xi = 0.15^2


%% Compute convergence (HMF)
n_test = 10*[4,8,16,32,64,128,256,512]'; 

% Set parameters:
model.eta_1 = 2.93;
model.rho_12 = -1.00;
model.zeta.values = 0.15^2;
model.K_1.alpha = -0.45;

S_contracts = struct;
S_contracts.k = (-0.2:0.01:0.05)';
S_contracts.ttm = 0.1*ones(size(S_contracts.k));

model.pricerSettings.N = 10000;
model.pricerSettings.numRepeat = 100;

figure;

% Compute small n's:
cm = parula(size(n_test,1)-2);
iv_mat_hmf = NaN(size(S_contracts.k,1),size(n_test,1));
for i=1:size(n_test,1)
    i
    model.pricerSettings.n(:) = n_test(i);
    pS = model.GetPrices(S_contracts,[]);
    iv_mat_hmf(:,i) = pS.p;
    if i > 2
        hold on;plot(pS.k,pS.p,'-','color',cm(i-2,:),'LineWidth',1.5);
        drawnow;
    end
end

% Compute with 10,000 steps:
% Load kernel approximation:
model.pricerSettings.n(:) = 100000;
model.pricerSettings.N = 1000;
model.pricerSettings.numRepeat = 1000;
pS_ref = model.GetPrices(S_contracts,[]);  
hold on;plot(pS_ref.k(idx),pS_ref.p(idx),'--','color','red','LineWidth',1.5); 
legend('\lfloornT\rfloor = 16','\lfloornT\rfloor = 32','\lfloornT\rfloor = 64',...
       '\lfloornT\rfloor = 128','\lfloornT\rfloor = 256',...
       '\lfloornT\rfloor = 512','\lfloornT\rfloor = 10000');
xlabel('Log-moneyness','interpreter','latex');
ylabel('Implied volatility','interpreter','latex');
title({'Rough Bergomi','(Hybrid multifactor)'},'interpreter','latex');
xlim([-0.35,0.15]);
ylim([0.09,0.38]);


%% Compute convergence (Hybrid TBSS)
n_test = 10*[4,8,16,32,64,128,256,512]'; 
T = 0.1;
H = 0.05;

S_contracts = struct;
S_contracts.k = (-0.3:0.01:0.1)';
S_contracts.ttm = 0.1*ones(size(S_contracts.k));
strikes = model.s0*exp(S_contracts.k);

K = model.K_1.DeepCopy();
K.f{1} = @(obj,t)(1);

figure;
cm = parula(size(n_test,1));
kappa = 1;
N_per = 10000;
N_runs = 100;
iv_mat_htbss = NaN(size(S_contracts.k,1),size(n_test,1));
for i=1:size(n_test,1)
    i
    dt = 1./floor2(n_test(i)*T);
    p = zeros(size(S_contracts.k));
    for j=1:N_runs
        [iv,se] = GetPricesRoughBergomi(model.s0,0,0,model.zeta,model.eta_1,model.rho_12,...
                                        model.K_1.alpha+0.5,strikes,...
                                        S_contracts.ttm,'HybridTBSS',N_per,n_test(i),...
                                        kappa,'turbo',false);
        p_new = bscall(s0,strikes,0,S_contracts.ttm,iv.^2,0);
        p = p + p_new;
    end
    p = (p./N_runs);
    iv_mat_htbss(:,i) = blsimpv(model.s0,model.s0*exp(S_contracts.k),0,S_contracts.ttm,p);
    hold on;plot(S_contracts.k,iv_mat_htbss(:,i),'-','color',cm(i,:),'LineWidth',1.5);
    drawnow;
    
end


% % Add reference solution:
% hold on;plot(pS_ref.k,pS_ref.p,'--','color','red','LineWidth',1.5); 
% 
% legend('\lfloornT\rfloor = 16','\lfloornT\rfloor = 32','\lfloornT\rfloor = 64',...
%         '\lfloornT\rfloor = 128','\lfloornT\rfloor = 256','\lfloornT\rfloor = 528',...
%         '\lfloornT\rfloor = 10000');
% xlabel('Log-moneyness','interpreter','latex');
% ylabel('Implied volatility','interpreter','latex');
% title({'Rough Bergomi','(Hybrid-TBSS)'},'interpreter','latex');
% xlim([-0.35,0.15]);
% ylim([0.09,0.38]);


%% Bias:
idx = [26;31;36];
bias_htbss = abs(pS_ref.p(idx) - iv_mat_htbss(idx,:));
bias_hmf = abs(pS_ref.p(idx) - iv_mat_hmf(idx,:));

figure;
cm = hsv(3);
%leg_names1 = {'Hybrid multifactor ($k$ = -0.05)','Hybrid multifactor, ATM','Hybrid multifactor, OTM call'};
%leg_names2 = {'Hybrid TBSS, OTM put','Hybrid TBSS, ATM','Hybrid TBSS, OTM call'};
for i=1:size(idx,1)
    plot(n_test/10',bias_htbss(i,:),'-','color',cm(i,:),'MarkerFaceColor',cm(i,:),...
     'LineWidth',1.25,'MarkerSize',5);hold on; 
end
for i=1:size(idx,1)
    plot(n_test/10',bias_htbss(i,:),'o-','color',cm(i,:),'MarkerFaceColor',cm(i,:),...
     'LineWidth',1.25,'MarkerSize',5);hold on; 
end
for i=1:size(idx,1)
    plot(n_test/10',bias_hmf(i,:),'--x','color',cm(i,:),'MarkerFaceColor',cm(i,:),...
         'LineWidth',1.25,'MarkerSize',7);hold on;        
end

set(gca, 'XScale', 'log');
xlabel('\lfloornT\rfloor');ylabel('Bias (in implied volatility)','interpreter','latex');
ylim([-0.001,0.03]);xlim([0,2*10^3]);
hleg = legend('-0.05','0.00','0.05','interpreter','latex');
legend();
title(hleg,'Log-moneyness')
ax = gca; 
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax = gca;
ax.XAxis.Exponent = 0;
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
title({'Rough Bergomi','(Hybrid multifactor vs. hybrid TBSS)'},'Interpreter','Latex');



