%% Initialize Script
clear;
project_folder = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(project_folder));

%% Description:
% We verify the fix solving the problem where round off errors would cause some time points
% to be wrongly truncated down.

%% Basic setup:
N = 25;
T = 2;
x0 = 0;
b = 0;
sigma = 1;
H = 0.1;
kappa = 1;

                                         
%% Check when actual time points are inputted:
n_test = [50,100,500,1000,5000]';

for i=1:size(n_test,1)
    n = n_test(i);
    delta = 1/n;
    K = @(t)(t.^(H-0.5));
    [c,gamm,K] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);

    tX = (0:delta:T);
    tU = tX;
    K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
    [Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,[],x0,b,sigma,gamm,c,kappa,'K',K,...
                                                 'precision','double','tX',tX,'tU',tU,...
                                                 'returnU',true);

    abs(max(Xpaths.t - Xpaths.t_truncated))
    abs(max(Xpaths.t - tX))
    
    abs(max(Upaths.t - Upaths.t_truncated))
    abs(max(Upaths.t - tU))    

end

%% Check when actual time points are inputted (part 2)
n_test = [50,100,500,1000,5000]';

for i=1:size(n_test,1)
    n = n_test(i);
    delta = 1/n;
    K = @(t)(t.^(H-0.5));
    [c,gamm,K] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);

    tX = (0:1:2*n)/n;
    tU = tX;
    K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
    [Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,[],x0,b,sigma,gamm,c,kappa,'K',K,...
                                                 'precision','double','tX',tX,'tU',tU,...
                                                 'returnU',true);

    abs(max(Xpaths.t - Xpaths.t_truncated))
    abs(max(Xpaths.t - tX))
    
    abs(max(Upaths.t - Upaths.t_truncated))
    abs(max(Upaths.t - tU))    

end


%% Mixed inputs:
n = 500;
delta = 1/n;
K = @(t)(t.^(H-0.5));
[c,gamm,K] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);

tX = [1/n,2/n,0.123123,0.2,0.52123];
tU = tX;
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,[],x0,b,sigma,gamm,c,kappa,'K',K,...
                                             'precision','single','tX',tX,'tU',tU,...
                                             'returnU',true);

[Xpaths.t;Xpaths.t_truncated]
[Upaths.t;Upaths.t_truncated]

