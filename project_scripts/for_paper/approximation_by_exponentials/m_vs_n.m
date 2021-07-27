%% Initialize Script
clear;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));


%% Test delta vs m: 
H = 0.1; 
T = 1;
K = @(t)(t.^(H-0.5)./gamma(H+0.5));
epsilon = 0.001;
n_test = [2^4,2^6,2^8,2^10,2^12,2^14]';
[err_l2_norm,err_unif,m_out] = deal(NaN(size(n_test)));
for i=1:size(n_test,1)
    i
    n = n_test(i);
    delta = T/n;
    dt = T/(n+1);
    tt = (delta:dt:T)';
    [c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon,...
                                                  'n',n,'delta',delta,'T',T,...
                                                  'error_measure','l2_normalised'); 
    err_l2_norm(i) = sqrt(sum((K(tt)-K_hat(tt)).^2)./sum(K(tt).^2));
    err_unif(i) = max(abs(K(tt)-K_hat(tt))./K(tt));
    m_out(i) = size(gamm,1);
    plot(tt,K(tt));hold on;plot(tt,K_hat(tt),'--');
end

[n_test,err_l2_norm,err_unif,m_out]


%% Plots (fractional):
H_test = [0.01,0.1,0.2,0.4];

cm = hsv(size(H_test,2));

dot_style = {'x','o','*','+'};

g = @(x)(log2(x));

figure;
for j=1:size(H_test,2)
    H = H_test(j);
    T = 1;
    K = @(t)(t.^(H-0.5)./gamma(H+0.5));
    epsilon = 0.001;
    n_test = [2^4,2^6,2^8,2^10,2^12,2^14]';
    [err_l2_norm,err_unif,m_out] = deal(NaN(size(n_test)));
    for i=1:size(n_test,1)
        i
        n = n_test(i);
        delta = T/n;
        dt = T/(n+1);
        tt = (delta:dt:T)';
        [c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon,...
                                                      'n',n,'delta',delta,'T',T,...
                                                      'error_measure','l2_normalised'); 
        err_l2_norm(i) = sqrt(sum((K(tt)-K_hat(tt)).^2)./sum(K(tt).^2));
        err_unif(i) = max(abs(K(tt)-K_hat(tt))./K(tt));
        m_out(i) = size(gamm,1);
    end
    
    max(err_l2_norm)

    plot(g(n_test),m_out,dot_style{j},'DisplayName',['\alpha = ', num2str(H-0.5)],'color',cm(j,:));hold on;
    b = regress(m_out,g(n_test));
    plot(g(n_test),b*g(n_test),'color',cm(j,:),'DisplayName',[num2str(b),'\cdot log_2(n)']);hold on;
end

xlabel('$log_2(n)$','interpreter','latex');ylabel('$m(n)$','interpreter','latex');
legend();

set(get(gca,'YLabel'),'Rotation',0)
xlim([3,15]);
ylim([0,20]);

%% Plots (power law):
beta_test = [5,10,25,50];

cm = hsv(size(beta_test,2));

dot_style = {'x--','o--','*--','+--'};

g = @(x)(log2(x));

figure;
for j=1:size(beta_test,2)
    T = 1;
    K = @(t)((1+t).^(-beta_test(j)));
    epsilon = 0.001;
    n_test = [2^4,2^6,2^8,2^10,2^12,2^14]';
    [err_l2_norm,err_unif,m_out] = deal(NaN(size(n_test)));
    for i=1:size(n_test,1)
        i
        n = n_test(i);
        delta = T/n;
        dt = T/(n+1);
        tt = (delta:dt:T)';
        [c,gamm,K,K_hat,exitflag] = ApproximateKernel('BM2005','K',K,'epsilon',epsilon,...
                                                      'n',n,'delta',delta,'T',T,...
                                                      'error_measure','l2_normalised'); 
        err_l2_norm(i) = sqrt(sum((K(tt)-K_hat(tt)).^2)./sum(K(tt).^2));
        err_unif(i) = max(abs(K(tt)-K_hat(tt))./K(tt));
        m_out(i) = size(gamm,1);
    end
    
    max(err_l2_norm)

    plot(g(n_test),m_out,dot_style{j},'DisplayName',['\beta = ', num2str(beta_test(j))],...
         'color',cm(j,:));hold on;
end

xlabel('$n$','interpreter','latex');ylabel('$m(n)$','interpreter','latex');
legend();

set(get(gca,'YLabel'),'Rotation',0)
ylim([0,4])
xlim([3,15])


%% Plot power-laws:
tt = (0.001:0.001:1)';
figure;
for j=1:size(beta_test,2)
    K = @(t)((1+t).^(-beta_test(j)));
    plot(tt,K(tt));hold on;
end

