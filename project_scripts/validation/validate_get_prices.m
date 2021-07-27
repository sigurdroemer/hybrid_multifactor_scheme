%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%%
N = 10000;
n = 500;
delta = 1/n;
kappa = 1;
H = 0.1;
eta = 2.1;
xi = 0.15^2;
rho = -0.9;
y = 0;
q = 0;
s0 = 100;
T = 1;


%% Compare two schemes:
scheme = 'HybridMultifactor';
K = @(t)(t.^(H-0.5));
[c,gamm] = ApproximateKernel('BM2005','K',K,'delta',delta,'T',T,'n',n);
strikes = (70:5:130)';
ttm = ones(size(strikes))*T;
[iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,scheme,N,n,kappa,'c',c,'gamm',gamm);

figure;
plot(strikes,iv,'-','color','red');hold on;
plot(strikes,iv-se*1.96,'--','color','red');hold on;
plot(strikes,iv+se*1.96,'--','color','red');hold on;

scheme = 'HybridTBSS';
[iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,scheme,N,n,kappa);

plot(strikes,iv,'-','color','blue');hold on;
plot(strikes,iv-se*1.96,'--','color','blue');hold on;
plot(strikes,iv+se*1.96,'--','color','blue');hold on;


%% Compare turbo:
N = 50000;
seed = 123;
rng(seed);
[iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,scheme,N,n,kappa,'c',c,'gamm',gamm);
rng(seed);
[iv_2,se_2] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,scheme,N,n,kappa,'turbo',false,...
                                'c',c,'gamm',gamm);
                            
                            
figure;
plot(strikes,iv,'-','color','red');hold on;
plot(strikes,iv-se*1.96,'--','color','red');hold on;
plot(strikes,iv+se*1.96,'--','color','red');hold on;  

plot(strikes,iv_2,'-','color','blue');hold on;
plot(strikes,iv_2-se_2*1.96,'--','color','blue');hold on;
plot(strikes,iv_2+se_2*1.96,'--','color','blue');hold on;

                            
                            
%%     
y = CurveClass('gridpoints',0,'values',0)
[iv,se] = GetPricesRoughBergomi(s0,y,q,xi,eta,rho,H,strikes,ttm,scheme,N,n,kappa,'c',c,'gamm',gamm);
iv






                            