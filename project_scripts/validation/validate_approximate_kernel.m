%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%% Test relation between l2 and uniform error:
H = 0.1; % Try different values here...
n = 500;
T = 1;
delta = T/n;
dt = T/(n+1);
tt = (delta:dt:T)';
K = @(t)(t.^(H-0.5));
epsilon_test = 10.^(-(1:10)');
[err_l2_norm,err_unif,m_out] = deal(NaN(size(epsilon_test)));
for i=1:size(epsilon_test,1)
    [c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon_test(i),...
                                                  'n',n,'delta',delta,'T',T,...
                                                  'disable_warn',false,'error_measure',...
                                                  'l2_normalised');
    err_l2_norm(i) = sqrt(sum((K(tt)-K_hat(tt)).^2)./sum(K(tt).^2));
    err_unif(i) = max(abs(K(tt)-K_hat(tt))./K(tt));
    m_out(i) = size(gamm,1);
end

[epsilon_test,err_l2_norm,err_unif,m_out]


%% Test uniform error and # m:
H = 0.1; % Try different values here...
n = 500;
T = 1;
delta = T/n;
dt = T/(n+1);
tt = (delta:dt:T)';
K = @(t)(t.^(H-0.5));
epsilon_test = 10.^(-(1:10)');
[err_unif,m_out] = deal(NaN(size(epsilon_test)));
for i=1:size(epsilon_test,1)
    [c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon_test(i),...
                                                  'n',n,'delta',delta,'T',T,...
                                                  'disable_warn',false,'error_measure',...
                                                  'uniform_relative');
    err_unif(i) = max(abs(K(tt)-K_hat(tt))./K(tt));
    m_out(i) = size(gamm,1);
end

[epsilon_test,err_unif,m_out]

%% Test BM2005 under various H choices:
% Try both uniform and l2.
H_test = (0.01:0.01:0.5)';
epsilon = 0.0001;
T = 1;
n = 500;
delta = T/n;
dt = T/(n+1);
tt = (delta:dt:T)';

[err_l2_norm,err_unif] = deal(NaN(size(H_test)));
for i=1:size(H_test,1)
    i
    K = @(t)(t.^(H_test(i)-0.5));
    [c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon,...
                                                  'n',n,'delta',delta,'T',T,...
                                                  'error_measure','l2_normalised');
    err_l2_norm(i) = sqrt(sum((K(tt)-K_hat(tt)).^2)./sum(K(tt).^2));
    err_unif(i) = max(abs(K(tt)-K_hat(tt))./K(tt));
end

close all;
figure;
plot(H_test,err_l2_norm,'.')

%figure;
%plot(H_test,err_unif,'.')


%% Test standard cases: (constant and sum of exponentials)
T = 1;
n = 500;
delta = 1/500;
epsilon = 0.01;
K = @(t)(23*ones(size(t)));
[c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon,...
                                              'n',n,'delta',delta,'T',T);
[c,gamm]

K = @(t)(2*exp(-2.*t) + 4*exp(-0.1.*t));
[c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon,...
                                              'n',n,'delta',delta,'T',T);
[c,gamm]


%% Test special cases and settings:
T = 1;
n = 500;
delta = T/n;
dt = T/(n+1);
tt = (delta:dt:T)';
H = 0.1;

% Do 'm' input:
K = @(t)(t.^(H-0.5));
[c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'n',n,'delta',delta,'T',T,'m',20);

[c,gamm]

close all;
figure;
plot(tt,K(tt),'-','color','blue','LineWidth',1.5);hold on;
plot(tt,K_hat(tt),'--','color','red','LineWidth',1.5);

%% Test default choice:
T = 1;
n = 500;
delta = T/n;
dt = T/(n+1);
tt = (delta:dt:T)';
H = 0.1;

% Do 'm' input:
K = @(t)(t.^(H-0.5));
[c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'n',n,'delta',delta,'T',T);

[c,gamm]

close all;
figure;
plot(tt,K(tt),'-','color','blue','LineWidth',1.5);hold on;
plot(tt,K_hat(tt),'--','color','red','LineWidth',1.5);
xlim([-0.05,1.05])
%ylim([0.90,1.5])











