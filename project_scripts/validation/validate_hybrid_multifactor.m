%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%% Description:
% Here we validate that indeed the implicit and explicit methods get the (approximately) same
% trajectories.
% The below code can be used to check all branches against itself and to compare tracjectories
% against the non-explicit method. All looks good.

%% Check all branches:
% b,sigma is (1) matrix, (2) function, (3) constant
b_fixed = 0.0;

seed = 123445;

N = 1;
n = 500;
delta = 1/n;
T = 1;
x0 = 0;
H = 0.1;
alpha = H - 0.5;
kappa = 1;
K = @(t)(t.^(H-0.5));
[c,gamm] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);

K = KernelFunctionClass(1,alpha,@(obj,t)(1),10^(-12));

eta = 3.1;
sigma_fixed = sqrt(2*H)*eta;

M = floor(n*T+eps(n*T))+1;

b_mat = repmat(b_fixed,N,M-1);
b_fun = @(t,x)(b_fixed*ones(size(x)));

sig_mat = repmat(sigma_fixed,N,M-1);
sig_fun = @(t,x)(sigma_fixed*ones(size(x)));

b_test = {b_fixed,b_mat,b_fun};
sig_test = {sigma_fixed,sig_mat,sig_fun};

rng(seed);
[Xpaths_0,Upaths_0,dW] = HybridMultifactorScheme(N,n,T,x0,b_test{1},sig_test{1},...
                                                gamm,c,kappa,'K',K,'returnU',true,...
                                                'explicit',true); 
figure;
plot(Xpaths_0.values.');hold on;
for i=1:size(b_test,2)
    for j=1:size(sig_test,2)
        j
        rng(seed);
        [Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b_test{i},sig_test{j},...
                                                     gamm,c,kappa,'K',K,'returnU',true,...
                                                     'explicit',false); 
                                                 
        plot(Xpaths.values.');hold on;
        %max(max(abs(Xpaths.values - Xpaths_0.values)))
        %max(max(max(abs(Upaths.values - Upaths_0.values))))     
                
    end
end









%% Check basics:
N = 25000;
n = 500;
delta = 1/n;
T = 2;
x0 = 0;
b = 0;
sigma = 1;
kappa = 1;
[c,gamm] = ApproximateKernel('BM2005','K',@(t)(K.Eval(t)),'T',T,'delta',delta,'n',n);
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K,...
                                             'precision','single');

%plot(Xpaths.t,Xpaths.values)

% Check mean and variance:
mu = mean(Xpaths.values);
sig = std(Xpaths.values);
figure;
plot(Xpaths.t,[mu-1.96*sig./sqrt(N);mu;mu+1.96*sig./sqrt(N)]);

figure;
plot(Xpaths.t,sig);hold on;
plot(Xpaths.t,sqrt(Xpaths.t.^(2*H)./(2*H)))

% Check U-factors in kappa = 0 case:
[Xpaths,Upaths,dW] = HybridMultifactorScheme(1,n,T,x0,b,sigma,gamm,c,0,'K',K,...
                                            'precision','single','returnU',true);

plot(sum(squeeze(Upaths.values.*c'),1),'-');hold on;
plot(Xpaths.values,'-');

Xpaths.values - sum(squeeze(Upaths.values.*c'),1)


% Check ad hoc time points:
rng(123);
tX = [0;0.1;0.10000001;0.2;0.6].';
tU = 1;
[Xpaths,Upaths,dW] = HybridMultifactorScheme(1,n,[],x0,b,sigma,gamm,c,kappa,'K',K,...
                                            'precision','single','returnU',true,'tX',tX,'tU',tU);

Xpaths.t-Xpaths.t_truncated

rng(123);
[Xpaths_2,Upaths_2,dW] = HybridMultifactorScheme(1,n,1,x0,b,sigma,gamm,c,kappa,'K',K,...
                                                 'precision','single','returnU',true);

idx_tX = find(ismember(Xpaths_2.t,tX));
idx_tU = find(ismember(Upaths_2.t,tU));

Xpaths_2.values(idx_tX)
Xpaths.values

Upaths.values
Upaths_2.values(:,:,idx_tU)

% Check dW output:
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K,...
                                             'precision','single','returndW',true);
                                        
figure;
mu = mean(dW).';
sig = std(dW).';
plot([mu-1.96*sig./sqrt(N),mu,mu+1.96*sig./sqrt(N)]);

figure;
plot(std(dW).^2)
1/n

% Z input:
N = 10;
kappa = 1;
M = floor2(n*T)+1;
Z = randn((M-1)*N,kappa+1);
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K,...
                                             'precision','single','Z',Z);
figure;
plot(Xpaths.values.');

% Positive paths:
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K,...
                                             'precision','single','positive',true);
plot(Xpaths.values.');


% Generate SIGMA from outside:
N = 10;
kappa = 1;
SIGMA = GetVolterraCovarianceMatrix({K},kappa,1/n);
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K,...
                                             'precision','single','SIGMA',SIGMA);                                   
plot(Xpaths.values.')

% Check U-factors correct even if kappa > 0
N = 10;
n = 500;
delta = 1/n;
T = 1;
x0 = 0;
b = 0;
sigma = 1;
H = 0.1;
K = @(t)(t.^(H-0.5));
[c,gamm] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);
K = KernelFunctionClass(1,alpha,@(obj,t)(1),10^(-12));


rng(123);
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,0,'K',K,...
                                             'precision','double','returnU',true);
rng(123);
[Xpaths_2,Upaths_2,dW] = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,4,'K',K,...
                                             'precision','double','returnU',true);                                         

max(max(max(abs(Upaths.values - Upaths_2.values))))

%% Check all branches:
% b,sigma is (1) matrix, (2) function, (3) constant
b_fixed = 0.1;
sigma_fixed = 1;
seed = 123;

N = 10;
n = 500;
delta = 1/n;
T = 1;
x0 = 0;
H = 0.1;
kappa = 1;
K = @(t)(t.^(H-0.5));
[c,gamm] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
M = sum((0:1:floor(n*T)+1)/n <= T);

b_mat = repmat(b_fixed,N,M-1);
b_fun = @(t,x)(b_fixed*ones(size(x)));

sig_mat = repmat(sigma_fixed,N,M-1);
sig_fun = @(t,x)(sigma_fixed*ones(size(x)));

b_test = {b_fixed,b_mat,b_fun};
sig_test = {sigma_fixed,sig_mat,sig_fun};

rng(seed);
[Xpaths_0,Upaths_0,dW] = HybridMultifactorScheme(N,n,T,x0,b_test{1},sig_test{1},...
                                             gamm,c,kappa,'K',K,'returnU',true); 

disp('Starting test');
for i=1:size(b_test,2)
    for j=1:size(sig_test,2)
        rng(seed);
        [Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,x0,b_test{i},sig_test{j},...
                                                     gamm,c,kappa,'K',K,'returnU',true); 
        max(max(abs(Xpaths.values - Xpaths_0.values)))
        max(max(max(abs(Upaths.values - Upaths_0.values))))
                
    end
end


plot(Xpaths.values.')


%% Test against rough Heston
% Will also verify the drift.
v_0 = 0.30.^2;v_bar=0.30.^2;alpha = 0.9;lambda = 2;xi = .5;rho = -0.7;
H = alpha - 0.5;
s0 = 100;strikes = (70:5:130)';T = 1;call = true;
%%

% Compute true prices:
[price, iv] = NumericalIntegrationRoughHeston(s0,v_0,alpha,lambda,...
                                              v_bar,xi,rho,call,strikes,T,...
                                              'disp_iter',true);

%% Monte Carlo:
N = 50000;
n = 500;
delta = 1/n;
dt = delta;
kappa = 5;
K = @(t)(t.^(alpha-1)./gamma(alpha));
b = @(t,x)(lambda.*(v_bar - x));
sigma = @(t,x)(xi.*sqrt(x));
[c,gamm] = ApproximateKernel('BM2005','K',K,'T',T,'delta',delta,'n',n);
K = KernelFunctionClass(1,alpha-1,@(obj,t)(1./gamma(alpha)),10^(-12));
[Xpaths,Upaths,dW] = HybridMultifactorScheme(N,n,T,v_0,b,sigma,gamm,c,kappa,'K',K,...
                                             'positive',true,'returndW',true); 

close all;
figure;
plot(sqrt(Xpaths.values(1:5,:))');

% Subset time points:
v = Xpaths.values;

% Simulate correlated Brownian motion:
dW_perp = normrnd(0,1,size(dW)).*sqrt(dt);
dW_1 = rho.*dW + sqrt(1-rho^2).*dW_perp;

corr(dW_1(:,1),dW(:,1))

% Compute asset price:
logs = zeros(size(Xpaths.values));
logs(:,1) = log(s0);
for i=2:size(logs,2)
    logs(:,i) = logs(:,i-1) + sqrt(v(:,i-1)).*dW_1(:,i-1) - 0.5.*v(:,i-1).*dt;
end

s = exp(logs);

figure;
plot(s(1:10,:).');

figure;
histogram(s(:,end),100)

figure;
plot(mean(s))

% Compute option prices:
[mc_price,mc_se_price,mc_se_iv,mc_iv] = deal(NaN(size(strikes)));

for i=1:size(strikes,1)
    X = max(s(:,end)-strikes(i),0);
    mc_price(i) = mean(X);
    mc_se_price(i) = std(X)./sqrt(N);
    mc_iv(i) = blsimpv(s0,strikes(i),0,T,mc_price(i));
    bs_vega = BSGreek('vega',[],s0,strikes(i),0,T,mc_iv(i),0);
    mc_se_iv(i) = mc_se_price(i)./bs_vega;
    
end


%%
close all;
figure;
plot(strikes,mc_iv,'-','color','red');hold on;
plot(strikes,mc_iv-1.96*mc_se_iv,'--','color','red');hold on;
plot(strikes,mc_iv+1.96*mc_se_iv,'--','color','red');hold on;
plot(strikes,iv,'-','color','blue');



%% Test differential equation:
% X'(t) = aX(t)dt => X(0)exp(a*t)

a = 2;
x0 = 1.2;
solution = @(t)(x0.*exp(a.*t));
K = KernelFunctionClass(1,0,@(obj,t)(1),10^(-12));
c = 1;
gamm = 0;
b = @(t,x)(a.*x);
sigma = 0;
kappa = 4;
n = 500;
T = 1;
N = 1;

Xpaths = HybridMultifactorScheme(N,n,T,x0,b,sigma,gamm,c,kappa,'K',K); 

plot(Xpaths.t,Xpaths.values)
hold on;
plot(Xpaths.t,solution(Xpaths.t))











                                        