%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));

%% rBergomi smiles:
rho = -1;
H = 0.05;
eta = 2.93;
xi = 0.15^2;
s0 = 100;
y = 0;
q = 0;

% Set strike and expiries:
strikes1 = (50:1:109)';
ttm1 = ones(size(strikes1));
strikes2 = (82:1:103)';
ttm2 = 0.1*ones(size(strikes2));


%% Plot maturity nr. 1:
rng(123);

% Find kernel:
n = 500;
[c,gamm] = ApproximateKernel('BM2005','K',@(t)(t.^(H-0.5)),'T',max(ttm1),'epsilon',0.001,...
                             'delta',1/n,'n',500,'error_measure','l2_normalised');
[c,gamm]

% Line colors:
cm = lines(5);

% Exact simulation:
% To avoid running out of memory we have to rerun it many times.
N_per = 1000;
N_runs = 1000;
n = 500;
kappa = n;
p = zeros(size(strikes1));
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
SIGMA = GetVolterraCovarianceMatrix({K},kappa,1/n);
scheme = 'HybridTBSS';
for i=1:N_runs
    i
    [iv_exact,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes1,ttm1,scheme,N_per,n,...
                                    kappa,'turbo',false,'conv_method','fft','SIGMA',SIGMA);
    p = p + bscall(s0,strikes1,0,ttm1,iv_exact.^2,0);
end
p = p./N_runs;
iv_exact = blsimpv(s0,strikes1,0,ttm1,p);

figure;
plot(log(strikes1./s0),iv_exact,'-','Color',[0.6 0.6 0.6],'LineWidth',4.5);hold on;

% Hybrid-TBSS:
N_runs = 100;
N_per = 10000;
scheme = 'HybridTBSS';
kappa = 1;
p = zeros(size(strikes1));
for i=1:N_runs
    i
    [iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes1,ttm1,scheme,...
                                    N_per,n,kappa,'turbo',false);
    p = p + bscall(s0,strikes1,0,ttm1,iv.^2,0);
end
p = p./N_runs;
iv = blsimpv(s0,strikes1,0,ttm1,p);
plot(log(strikes1./s0),iv,'-','Color',cm(4,:),'LineWidth',2);hold on

% Hybrid multifactor:
N_runs = 100;
N_per = 10000;
scheme = 'HybridMultifactor';
kappa_test = [0,1]; % Add kappa = 2,3 to test these too.
for i=1:size(kappa_test,2)
    
    p = zeros(size(strikes1));
    for j=1:N_runs
        j
        [iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes1,ttm1,scheme,N_per,n,...
                                        kappa_test(i),'c',c,'gamm',gamm,'turbo',false,...
                                        'explicit',false);
        p = p + bscall(s0,strikes1,0,ttm1,iv.^2,0);
    end
    p = p./N_runs;
    iv = blsimpv(s0,strikes1,0,ttm1,p);
    
    if i == 2
        plot(log(strikes1./s0),iv,'--','Color',cm(i+1,:),'LineWidth',2);hold on;
    else
        plot(log(strikes1./s0),iv,'--','Color',cm(i,:),'LineWidth',2);hold on;
    end
    
end

legend('Exact','Hybrid-TBSS ($\kappa = 1$)','$\kappa = 0$','$\kappa = 1$','$\kappa = 2$',...
       'Interpreter','Latex');
xlabel('Log-moneyness','Interpreter','Latex');
ylabel('Implied volatility','Interpreter','Latex');
title(['T = ', num2str(ttm1(1))],'Interpreter','Latex','FontSize',11);
legend boxoff;


%% Plot maturity nr. 2:
rng(123);

% Find kernel:
n = 5000;
[c,gamm] = ApproximateKernel('BM2005','K',@(t)(t.^(H-0.5)),'T',max(ttm2),'epsilon',0.001,...
                             'delta',1/n,'n',5000,'error_measure','l2_normalised');
[c,gamm]

% Line colors:
cm = lines(5);

% Exact simulation:
% To avoid running out of memory we have to rerun it many times.
scheme = 'HybridTBSS';
N_per = 1000;
N_runs = 1000;
kappa = floor(n*ttm2(1));
p = zeros(size(strikes2));
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
SIGMA = GetVolterraCovarianceMatrix({K},kappa,1/n);
for i=1:N_runs
    i
    [iv_exact,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes2,ttm2,scheme,N_per,n,...
                                    kappa,'turbo',false,'SIGMA',SIGMA,'conv_method','fft');
    p = p + bscall(s0,strikes2,0,ttm2,iv_exact.^2,0);
end
p = p./N_runs;
iv_exact = blsimpv(s0,strikes2,0,ttm2,p);

figure;
plot(log(strikes2./s0),iv_exact,'-','Color',[0.6 0.6 0.6],'LineWidth',4.5);hold on;


% Hybrid-TBSS:
scheme = 'HybridTBSS';
kappa = 1;
p = zeros(size(strikes2));
for i=1:N_runs
    i
    [iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes2,ttm2,scheme,...
                                    N_per,n,kappa,'turbo',false);
    p = p + bscall(s0,strikes2,0,ttm2,iv.^2,0);
end
p = p./N_runs;
iv = blsimpv(s0,strikes2,0,ttm2,p);
plot(log(strikes2./s0),iv,'-','Color',cm(4,:),'LineWidth',2);hold on


% Hybrid multifactor:
scheme = 'HybridMultifactor';
kappa_test = [0,1]; % Add kappa = 2,3 to test these too.
for i=1:size(kappa_test,2)
    
    p = zeros(size(strikes2));
    for j=1:N_runs
        j
        [iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes2,ttm2,scheme,N_per,n,...
                                        kappa_test(i),'c',c,'gamm',gamm,'turbo',false,...
                                        'explicit',false);
        p = p + bscall(s0,strikes2,0,ttm2,iv.^2,0);
    end
    p = p./N_runs;
    iv = blsimpv(s0,strikes2,0,ttm2,p);
    
    if i == 2
        plot(log(strikes2./s0),iv,'--','Color',cm(i+1,:),'LineWidth',2);hold on;
    else
        plot(log(strikes2./s0),iv,'--','Color',cm(i,:),'LineWidth',2);hold on;
    end
    
end

legend('Exact','Hybrid-TBSS ($\kappa = 1$)','$\kappa = 0$','$\kappa = 1$','$\kappa = 2$',...
       'Interpreter','Latex');
xlabel('Log-moneyness','Interpreter','Latex');
ylabel('Implied volatility','Interpreter','Latex');
title(['T = ', num2str(ttm2(1))],'Interpreter','Latex');%
legend boxoff;


