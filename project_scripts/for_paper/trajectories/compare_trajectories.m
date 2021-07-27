%% Description::
% Here we illustrate simulated trajectories using each scheme and that with various settings.

%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
addpath(genpath(project_folder));

%% Main settings
alpha = 0.55;
n = 250;
T = 1;
K = @(t)(t.^(alpha-1));
Kclass = KernelFunctionClass(1,alpha-1,@(obj,t)(1),10^(-12));
kappa_max = floor(n*T);

%% Compute paths:
% Generate underlying numbers:
rng(1234);
W_bold = SampleVolterra(chol_mod(ones(3,3)),{Kclass,Kclass,Kclass},...
                        [0,1,kappa_max],1/n,false,N,n,[]);

%% Run exact simulation:
gamm = 0;
c = 0;
Xpaths_exact = HybridMultifactorScheme(N,n,T,0,0,1,gamm,c,kappa_max,...
                                       'W_bold',W_bold{3},'explicit',false,'K',Kclass);

% Run hybrid-exponential scheme:
[c,gamm] = ApproximateKernel('BM2005','delta',1/n,'T',T,'K',K,'epsilon',0.001,...
                             'error_measure','l2_normalised','n',250);
[Xpaths_exp,Upaths] = HybridMultifactorScheme(N,n,T,0,0,1,gamm,c,kappa_use,...
                                             'W_bold',W_bold{2},'returnU',true,...
                                              'explicit',false,'K',Kclass); 
[c,gamm] = ApproximateKernel('BM2005','delta',1/n,'T',T,'K',K,'epsilon',0.001,...
                             'error_measure','l2_normalised','n',250);                             
[Xpaths_exp_2,Upaths_2] = HybridMultifactorScheme(N,n,T,0,0,1,gamm,c,0,...
                                                 'W_bold',W_bold{1},'returnU',true,...
                                                 'explicit',false,'K',Kclass);

% Plot:
figure;
cm = lines(3);
plot(Xpaths_exact.t,Xpaths_exact.values,'LineWidth',1.15,'color',cm(1,:));hold on;
plot(Xpaths_exp_2.t,Xpaths_exp_2.values,'-','LineWidth',1.15,'color',cm(3,:));hold on;
plot(Xpaths_exp.t,Xpaths_exp.values,'--','LineWidth',1.15,'color',cm(2,:));hold on;
hleg = legend('location','SouthWest');
hleg.Title.String = 'Scheme';
xlabel('$t$','interpreter','latex');
ylabel('$Y(t)$','interpreter','latex');
legend('Exact','Hybrid multifactor ($\kappa = 0$)','Hybrid multifactor ($\kappa = 1$)',...
       'interpreter','latex');
set(get(gca,'YLabel'),'Rotation',0)
ylim([-9,9]);
xlim([-0.02,1.02])


