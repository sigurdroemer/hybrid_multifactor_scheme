%% Initialize Script
clear;
project_folder = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(project_folder));

%% Description:
% In this script we verify that our error fix of the problem with kappa > 1 has worked.

%% Check deterministic equation:
T = 1;
n = 1000;
delta = 1/n;
t = (delta:delta:T)';
N = 1;
kappa_test = (1:1:20)';

% Model:
v0 = 0.15^2;
lambda = 2;
theta = v0;
H = 0.25;
b = @(t,x)(lambda.*(theta-x));
sigma = 0;

for i=1:size(kappa_test,1)
    %i
    % Approximate kernel:
    kappa = kappa_test(i);
    K = @(t)(t.^(H - 0.5));
    [c,gamm] = ApproximateKernel('BM2005','K',K,'epsilon',0.01,'n',n,...
                                 'delta',kappa*delta,'T',T,'approx_roots',true,'Nz',10^4);
    K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
    Xpaths = HybridMultifactorScheme(N,n,T,v0,b,sigma,gamm,c,kappa,'K',K,...
                                     'precision','double','positive',true,...
                                     'returnU',false);      
                                 
    if i > 1
        max(abs(Xpaths.values-Xprev))
    end
	Xprev = Xpaths.values;
end

%% Check stochastic equation:
% We verify both mean of X and that of U:

N = 10000;
eta = 0.1;
sigma = @(t,x)(eta.*sqrt(x));
kappa_test = (1:10)';

for i=1:size(kappa_test,1)
    i
    % Approximate kernel:
    kappa = kappa_test(i);
    K = @(t)(t.^(H - 0.5));
    [c,gamm] = ApproximateKernel('BM2005','K',K,'epsilon',0.01,'n',n,...
                                     'delta',kappa*delta,'T',T,'approx_roots',true,'Nz',10^4);
    K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
    Xpaths = HybridMultifactorScheme(N,n,T,v0,b,sigma,gamm,c,kappa,'K',K,...
                                     'precision','double','positive',true,...
                                     'returnU',false);      
                                 
    plot(Xpaths.t,mean(Xpaths.values));hold on;
end

figure;
for i=1:size(kappa_test,1)
    i
    % Approximate kernel:
    kappa = kappa_test(i);
    K = @(t)(t.^(H - 0.5));
    [c,gamm,K,K_hat] = ApproximateKernel('BM2005','K',K,'epsilon',0.01,'n',n,...
                                     'delta',kappa*delta,'T',T,'approx_roots',true,'Nz',10^4);
    K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
    [Xpaths,Upaths] = HybridMultifactorScheme(N,n,T,v0,b,sigma,gamm,c,kappa,'K',K,...
                                     'precision','double','positive',true,...
                                     'returnU',true);      
                                 
    plot(Upaths.t,squeeze(mean(Upaths.values)));hold on;
    
    max(max(squeeze(std(Upaths.values))./sqrt(N)))
end

figure;
plot(squeeze(Upaths.values(22,:,:))');

%% Check pathwise:
N = 1;

for i=1:size(kappa_test,1)
    i
    % Approximate kernel:
    kappa = kappa_test(i);
    K = @(t)(t.^(H - 0.5));
    [c,gamm] = ApproximateKernel('BM2005','K',K,'epsilon',0.01,'n',n,...
                                     'delta',kappa*delta,'T',T,'approx_roots',true,'Nz',10^4);
    rng(312);
    K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
    Xpaths = HybridMultifactorScheme(N,n,T,v0,b,sigma,gamm,c,kappa,'K',K,...
                                     'precision','double','positive',true,...
                                     'returnU',false);      
                                 
    plot(Xpaths.t,Xpaths.values);hold on;
end





