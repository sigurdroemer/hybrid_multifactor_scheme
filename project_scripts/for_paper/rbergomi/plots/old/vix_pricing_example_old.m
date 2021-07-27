%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));


%% Basic setup:
alpha = -0.40;
beta = -0.40;
K1 = KernelFunctionClass(1,alpha,{@(obj,t)( sqrt(2*obj.alpha+1) )},10^(-12));
K2 = KernelFunctionClass(1,beta,{@(obj,t)( sqrt(2*obj.alpha+1) )},10^(-12));
rho_12 = -0.6;
rho_23 = 1;
rho_13 = -0.6;
xi0 = 0.15^2;
eta = 2.5;
nu = 3.5;
theta = 0.2;
model = MixedrBergomiClass('c',0,'eta_1',eta,'eta_2',nu,'theta',theta,'K_1',K1,'K_2',K2,...
                           'rho_12',rho_12,'rho_23',rho_23,'rho_13',rho_13,'zeta',xi0);
model.pricerSettings.simulation.kernel_approximation.error_measure = 'l2_normalised';
model.pricerSettings.price_estimation.S.control_variate = 'none';


%% Can we estimate mean well enough?
N = 10000;
paths = model.Simulate(N,5000,{'V'},{0.1});    
mu = mean(paths.V.values);
se = std(paths.V.values)./sqrt(N);
[mu-1.96*se,mu,mu+1.96*se]


%% Compute full smile:
% For inspection and selection of strikes:
VIX_contracts = struct;
VIX_contracts.ttm = 0.1;
VIX_contracts.q = (0.01:0.01:0.99)';
model.pricerSettings.N = 10000;
model.pricerSettings.numRepeat = 10;
[~,pV] = model.GetPrices([],VIX_contracts);

figure;plot(pV.K,[pV.p-pV.se*1.96,pV.p,pV.p+pV.se*1.96]);



%% Setup for numerical examples:
t = 0.1;
nv_test = [2,4,8,16,32,64,128,256];
strikes = [0.10,0.20,0.30];
DELTA = 1/12;
N = 10000;
Npts = 101;


%% Verify covariance matrix computations (i.e. how many points?):
% Do for both methods!


%% Exact sampling:
rng(1234);
n_run = 1;

[p1_exact,p2_exact,p3_exact,tt_exact] = deal(NaN(n_run,size(nv_test,2)));
for j=1:n_run
    for i=1:size(nv_test,2)
        nv = nv_test(i);
        disp(['(i,j,nv) = ', '(',num2str(i),',',num2str(j),',',num2str(nv),')']);
        
        % Sample N(0,1)'s:
        Z = randn(2*nv,N);
        
        tic;
        
        % VIX integration points:
        du = DELTA/(nv-1);
        u = (t:du:t+DELTA)';          

        % Compute covariance matrix:
        SIGMA11 = GetCovMatrix(K1,K1,u,t,Npts,true);
        SIGMA12 = rho_23*GetCovMatrix(K1,K2,u,t,Npts,false);
        SIGMA22 = GetCovMatrix(K2,K2,u,t,Npts,true);
        
        SIGMA = [[SIGMA11,SIGMA12];...
                 [SIGMA12',SIGMA22']];
             
        % Factorize:
        A = chol_mod(SIGMA,10.^(-(16:-1:6)));
        
        % Sample:
        Y = A*Z;
        
        % Compute forward variances:
        drift1 = 0.5*eta^2*(u.^(2*alpha+1) - (u-t).^(2*alpha+1));
        drift2 = 0.5*nu^2*(u.^(2*beta+1) - (u-t).^(2*beta+1));
        xi = xi0*(theta*exp(eta*Y(1:nv,:)-drift1) + (1-theta)*exp(nu*Y(nv+1:end,:)-drift2));
        
        % Run to validate mean:        
        %mu = mean(xi')';
        %se = std(xi')'./sqrt(N);
        %figure;plot([mu-1.96*se,mu,mu+1.96*se]);
        %drawnow;
        
        % Compute VIX:
        VIX = sqrt((1/DELTA)*(0.5*(xi(1,:) + xi(end,:)) + sum(xi(2:end-1,:),1))*du);
        Vfut = mean(VIX); 
        
        % Strike 1:
        payoffs = max(VIX - strikes(1),0);
        p1_exact(j,i) = blsimpv(Vfut,strikes(1),0,t,mean(payoffs));
        
        % Strike 2:
        payoffs = max(VIX - strikes(2),0);
        p2_exact(j,i) = blsimpv(Vfut,strikes(2),0,t,mean(payoffs));
        
        % Strike 3:
        payoffs = max(VIX - strikes(3),0);
        p3_exact(j,i) = blsimpv(Vfut,strikes(3),0,t,mean(payoffs));        
        
        tt_exact(j,i) = toc;
        
    end
end


%% Hybrid multifactor approach (exact sampling):
rng(123);
n_run = 1;

[p1_hmf_exact,p2_hmf_exact,p3_hmf_exact,tt_hmf_exact] = deal(NaN(n_run,size(nv_test,2)));
for j=1:n_run
    % Pre-simulate random numbers:
    Z = randn(15,N);    
    for i=1:size(nv_test,2)
        nv = nv_test(i);
        disp(['(i,j,nv) = ', '(',num2str(i),',',num2str(j),',',num2str(nv),')']);
        
        tic;
        
        % VIX integration points:
        du = DELTA/(nv-1);
        u = (t:du:t+DELTA)';             
        
        % Fit sum-of-exponentials approximations:
        Npts_Kapprox = ceil((t+DELTA)*500);
        [c1,gamm1] = ApproximateKernel('BM2005','K',@(t)(K1.Eval(t)),'T',t+DELTA,...
                                       'delta',du,'N',Npts_Kapprox,'epsilon',0.001,...
                                       'error_measure','l2_normalised');  
        [c2,gamm2] = ApproximateKernel('BM2005','K',@(t)(K2.Eval(t)),'T',t+DELTA,...
                                       'delta',du,'N',Npts_Kapprox,'epsilon',0.001,...
                                       'error_measure','l2_normalised');                                  
        
        m1 = size(c1,1);
        m2 = size(c2,1);
        
        % Compute covariance matrix:
        gamm = [gamm1;gamm2];
        gamm_mat = repmat(gamm,1,size(gamm,1));
        
        % SIGMA_U:
        SIGMA_U = (1./(gamm_mat+gamm_mat')).*(1-exp(-(gamm_mat+gamm_mat')*t));
        SIGMA_U(m1+1:end,1:m1) = rho_23*SIGMA_U(m1+1:end,1:m1);
        SIGMA_U(1:m1,m1+1:end) = rho_23*SIGMA_U(1:m1,m1+1:end);

        % SIGMA_Y:
        SIGMA_Y = GetVolterraCovarianceMatrix({K1,K2},[1,1],t);
        SIGMA_Y = SIGMA_Y(2:3,2:3);
        SIGMA_Y(1,2) = rho_23*SIGMA_Y(1,2);
        SIGMA_Y(2,1) = rho_23*SIGMA_Y(2,1);
        
        % SIGMA_Y1_U1: 
        SIGMA_Y1_U1 = zeros(size(c1,1),1);
        for k=1:size(c1,1)
            Ktemp = KernelFunctionClass(1,0,@(obj,s)(exp(-gamm1(k)*s)),10^(-12));
            SIGMA_Y1_U1(k) = GetVolterraCovariance(K1,0,t,t,Ktemp,0,t,t);
        end
        
        % SIGMA_Y1_U2:
        SIGMA_Y1_U2 = zeros(size(c2,1),1);
        for k=1:size(c2,1)
            Ktemp = KernelFunctionClass(1,0,@(obj,s)(exp(-gamm2(k)*s)),10^(-12));
            SIGMA_Y1_U2(k) = rho_23*GetVolterraCovariance(K1,0,t,t,Ktemp,0,t,t);
        end        
        
        % SIGMA_Y2_U1:
        SIGMA_Y2_U1 = zeros(size(c1,1),1);
        for k=1:size(c1,1)
            Ktemp = KernelFunctionClass(1,0,@(obj,s)(exp(-gamm1(k)*s)),10^(-12));
            SIGMA_Y2_U1(k) = rho_23*GetVolterraCovariance(K2,0,t,t,Ktemp,0,t,t);
        end                
        
        % SIGMA_Y2_U2:
        SIGMA_Y2_U2 = zeros(size(c2,1),1);
        for k=1:size(c2,1)
            Ktemp = KernelFunctionClass(1,0,@(obj,s)(exp(-gamm2(k)*s)),10^(-12));
            SIGMA_Y2_U2(k) = GetVolterraCovariance(K2,0,t,t,Ktemp,0,t,t);
        end               
        
        SIGMA_Y_U = [[SIGMA_Y1_U1;SIGMA_Y1_U2],[SIGMA_Y2_U1;SIGMA_Y2_U2]];
        
        SIGMA = [[SIGMA_Y,SIGMA_Y_U'];...
                 [SIGMA_Y_U,SIGMA_U]];
             
        A = chol_mod(SIGMA,10.^(-(16:-1:6)));
        rvTemp = A*Z(1:m1+m2+2,:);

        % Compute Y's:
        [Y1,Y2] = deal(zeros(nv,N));
        Y1(1,:) = rvTemp(1,:);
        Y2(1,:) = rvTemp(2,:);
        for l=2:size(u,1)
            Y1(l,:) = (c1'.*exp(-gamm1'*(u(l)-t)))*rvTemp(3:2+m1,:);
            Y2(l,:) = (c2'.*exp(-gamm2'*(u(l)-t)))*rvTemp(3+m1:end,:);
        end

        % Compute forward variances:
        drift1 = 0.5*eta^2*(u.^(2*alpha+1) - (u-t).^(2*alpha+1));
        drift2 = 0.5*nu^2*(u.^(2*beta+1) - (u-t).^(2*beta+1));
        xi = xi0*(theta*exp(eta*Y1-drift1) + (1-theta)*exp(nu*Y2-drift2));

        % Run to validate mean:        
        %mu = mean(xi')';
        %se = std(xi')'./sqrt(N);
        %figure;plot([mu-1.96*se,mu,mu+1.96*se]);
        
        % Compute VIX:
        VIX = sqrt((1/DELTA)*(0.5*(xi(1,:) + xi(end,:)) + sum(xi(2:end-1,:),1))*du);
        Vfut = mean(VIX); 
        
        % Strike 1:
        payoffs = max(VIX - strikes(1),0);
        p1_hmf_exact(j,i) = blsimpv(Vfut,strikes(1),0,t,mean(payoffs));
        
        % Strike 2:
        payoffs = max(VIX - strikes(2),0);
        p2_hmf_exact(j,i) = blsimpv(Vfut,strikes(2),0,t,mean(payoffs));
        
        % Strike 3:
        payoffs = max(VIX - strikes(3),0);
        p3_hmf_exact(j,i) = blsimpv(Vfut,strikes(3),0,t,mean(payoffs));

        tt_hmf_exact(j,i) = toc;
        
    end
end


%% Hybrid multifactor approach (scheme):
rng(12345);
n_run = 1;

[p1_hmf,p2_hmf,p3_hmf,tt_hmf] = deal(NaN(n_run,size(nv_test,2)));
for j=1:n_run
    % Sample (V,U1,U2):
    n = 5000;
    paths = model.Simulate(N,n,{'V','U1','U2'},{[t,t+DELTA],[t,t+DELTA],[t,t+DELTA]});    
    
    for i=1:size(nv_test,2)
        nv = nv_test(i);
        disp(['(i,j,nv) = ', '(',num2str(i),',',num2str(j),',',num2str(nv),')']);

        tic;
        
        % VIX integration points:
        du = DELTA/(nv-1);
        u = (t:du:t+DELTA)';        
        
        [Y1,Y2] = deal(zeros(N,nv-1));

        % Compute Y's:
        c1 = paths.U1.c;gamm1 = paths.U1.gamm;
        c2 = paths.U2.c;gamm2 = paths.U2.gamm;

        for l=2:size(u,1)
            Y1(:,l-1) = paths.U1.values(:,:,1)*(c1.*exp(-gamm1*(u(l)-t)));
            Y2(:,l-1) = paths.U2.values(:,:,1)*(c2.*exp(-gamm2*(u(l)-t)));
        end

        % Compute forward variances:
        xi = zeros(N,nv);
        drift1 = 0.5*eta^2*(u(2:end).^(2*alpha+1) - (u(2:end)-t).^(2*alpha+1))';
        drift2 = 0.5*nu^2*(u(2:end).^(2*beta+1) - (u(2:end)-t).^(2*beta+1))';
        xi(:,2:end) = xi0*(theta*exp(Y1-drift1) + (1-theta)*exp(Y2-drift2));
        xi(:,1) = paths.V.values(:,1);

        % Run to validate mean:        
        %mu = mean(xi)';
        %se = std(xi)'./sqrt(N);
        %figure;plot([mu-1.96*se,mu,mu+1.96*se]);
        %drawnow;
        
        % Compute VIX:
        VIX = sqrt((1/DELTA)*(0.5*(xi(:,1) + xi(:,end)) + sum(xi(:,2:end-1),2))*du);
        Vfut = mean(VIX); 
        
        % Strike 1:
        payoffs = max(VIX - strikes(1),0);
        p1_hmf(j,i) = blsimpv(Vfut,strikes(1),0,t,mean(payoffs));
        
        % Strike 2:
        payoffs = max(VIX - strikes(2),0);
        p2_hmf(j,i) = blsimpv(Vfut,strikes(2),0,t,mean(payoffs));
        
        % Strike 3:
        payoffs = max(VIX - strikes(3),0);
        p3_hmf(j,i) = blsimpv(Vfut,strikes(3),0,t,mean(payoffs));

        tt_hmf(j,i) = toc;
    end
end


%% Present results:
[p1_exact;p1_hmf_exact;p1_hmf]
[p2_exact;p2_hmf_exact;p2_hmf]
[p3_exact;p3_hmf_exact;p3_hmf]

round([mean(tt_exact,1);mean(tt_hmf_exact,1);mean(tt_hmf,1)]*1000)


%% Auxiliary function:
function SIGMA = GetCovMatrix(K1,K2,u,t,Npts,isCovM)

        SIGMA = NaN(size(u,1),size(u,1));
        
        if isCovM
            I = zeros(size(u,1)*size(u,1)/2+size(u,1)/2,1);
        else
            I = zeros(size(u,1)*size(u,1),1);
        end
        
        % Compute u-vectors:
        idx1 = repmat((1:size(u,1))',1,size(u,1)); % row index
        idx2 = repmat((1:size(u,1)),size(u,1),1); % column index
        if isCovM
            idxCompute = idx2 <= idx1;
            c1 = u(idx1(idxCompute));
            c2 = u(idx2(idxCompute));            
        else
            c1 = u(idx1(:));
            c2 = u(idx2(:));
        end

        idxSingular = c1 == t | c2 == t;
        
        % Compute integrals (non-singular):
        I(~idxSingular) = SimpsonsIntegration(@(s)(K1.Eval(c1(~idxSingular)-s)...
                                                 .*K2.Eval(c2(~idxSingular)-s)),0,t,Npts);
        
        % Compute integrals (singular):
        
        % We compute int_{t-dt}^t K1(c1-s) K2(c2-s) ds:
        dt = 10^(-4);
        idx1 = c1 == t & c2 == t;
        idx2 = c1 == t & c2 ~= t;
        idx3 = c1 ~= t & c2 == t;
        I(idx1) = (1./(K1.alpha+K2.alpha+1))*dt^(K1.alpha+K2.alpha+1);
        I(idx2) = (1./(K1.alpha+1))*dt^(K1.alpha+1)...
                  *(0.5*(K2.Eval(c2(idx2)-t)+K2.Eval(c2(idx2)-(t-dt))));
        I(idx3) = (1./(K2.alpha+1))*dt^(K2.alpha+1)...
                  *(0.5*(K1.Eval(c1(idx3)-t)+K1.Eval(c1(idx3)-(t-dt))));
        I(idxSingular) = sqrt(2*K1.alpha+1)*sqrt(2*K2.alpha+1)*I(idxSingular);
              
        % We compute int_{0}^{t-dt} K1(c1-s) K2(c2-s) ds for singular integrals:
        I(idxSingular) = I(idxSingular) + SimpsonsIntegration(@(s)(K1.Eval(c1(idxSingular)-s)...
                                                     .*K2.Eval(c2(idxSingular)-s)),0,t-dt,Npts);        
        
        % Reshape and add to covariance matrix:
        if isCovM
            SIGMA(idxCompute) = I;
        else
            SIGMA(:) = I;
        end

        if isCovM
            % Remaining entries follow by symmetry:
            SIGMA(isnan(SIGMA)) = 0;
            SIGMADiagZero =  SIGMA - diag(diag(SIGMA));
            SIGMA = SIGMA + SIGMADiagZero';            
        end
        
end
% 
% %% OLD STUFF:
% %% Load prices:
% load(['C:\Users\sigur\GoogleDrive\university\phd\research\rough_volatility_source_code\',...
%      '05_GitHub\hybrid_multifactor_scheme\project_scripts\for_paper\rbergomi\rHeston_prices.mat']);
% 
% strikes = (80:2.5:110)';
% S_contracts = struct;
% S_contracts.k = log(strikes/model.s0);
% S_contracts.ttm = 0.1*ones(size(strikes));
% S_contracts.mid = iv_rH(:,end);
% S_contracts.bid = S_contracts.mid;
% S_contracts.ask = S_contracts.mid;
% 
% %% Calibrate:
% diffChg = 0.001;
% maxIter = 500;
% options = optimoptions('lsqnonlin',...
%                        'Display','off','DiffMinChange',diffChg,'MaxIter',maxIter,...
%                        'TolX',10^(-10),'TolFun',10^(-10));                 
%                    
% parNames = {'eta_1','eta_2','rho_12','rho_13','zeta'};  
% par0 = {model.eta_1;model.eta_2;model.rho_12;model.rho_13;model.zeta.values}';
% lb = {0;0;-1;-1;0}';
% ub = {30;30;1;1;1}';
% 
% % Calibrate to SPX options:
% model.pricerSettings.n(:) = 5000;
% model.pricerSettings.N = 10000;
% model.pricerSettings.price_estimation.antithetic = false;
% model.pricerSettings.seed = 123;
% result = model.Calibrate(S_contracts,[],parNames,par0,lb,ub,...
%                           'options',options,...
%                           'dispIter',true,...
%                           'error_measure','dist_to_mid',...
%                           'plotFit',true); 
% model.pricerSettings.seed = [];
% 
% % Save model:
% save([fileparts(matlab.desktop.editor.getActiveFilename),'/mixed_rbergomi_calibrated.mat'],'model');
% 
% % Parameters: 
% % eta = 3.85, nu = 3.86, rho_23 = 0.5, theta = 0.70, zeta = 0.15^2.
% model.eta_1 = 3.85;
% model.eta_2 = 3.86;
% model.rho_23 = 0.5;
% model.theta.values = 0.7;
% model.zeta.values = 0.15^2;
