%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%%
N = 10000;
n = 500;
T = 1;
sigma = 1;
H = 0.1;
alpha = H - 0.5;
kappa = 1;
[Y,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,[],[],[],[],[]);


% Check mean and variance:
mu = mean(Y);
sig = std(Y);
figure;
plot(t,[mu-1.96*sig./sqrt(N);mu;mu+1.96*sig./sqrt(N)]);

figure;
plot(t,sig);hold on;
plot(t,sqrt(t.^(2*H)./(2*H)))

% Z input:
M = floor(n*T+eps(n*T))+1;
Z = randn((M-1)*N,kappa+1);
[Y,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,Z,[],[],[],[]);

% Check and compare conv_methods:
N = 10;
n = 200;
seed = 123;
rng(seed);[Y1,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,[],'fft',[],[],[]);
rng(seed);[Y3,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,[],'conv_loop',[],[],[]);

max(max(abs(Y1 - Y3)))

%% Determine cut-off between FFT and conv2: (OLD - not needed anymore)
% Conclusion: 1200
N = 1000;
n_test = (50:200:2000)';
[tt_fft,tt_conv2] = deal(NaN(size(n_test)));
n_run = 10;
for i=1:size(n_test,1)
    i
    tic;
    for j=1:n_run
        HybridTBSSScheme(N,n_test(i),T,sigma,alpha,kappa,[],'fft',[],[],[]);
    end
    tt_fft(i) = toc/n_run;
    
    tic;
    for j=1:n_run
        HybridTBSSScheme(N,n_test(i),T,sigma,alpha,kappa,[],'conv2',[],[],[]);
    end
    tt_conv2(i) = toc/n_run;    
    
end

figure;
plot(n_test,tt_fft,'o-');
hold on;
plot(n_test,tt_conv2,'o-');
legend('fft','conv2');


%% Check and compare cov_methods:
rng(seed);[Y1,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,[],[],[],'exact',[]);
rng(seed);[Y2,t,dW] = HybridTBSSScheme(N,n,T,sigma,alpha,kappa,[],[],[],'numerical_integration',[]);
max(max(abs(Y1 - Y2)))

