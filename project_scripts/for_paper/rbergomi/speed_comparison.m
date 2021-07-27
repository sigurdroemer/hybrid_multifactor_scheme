%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));


%% Settings:
H = 0.05;
K = KernelFunctionClass(1,H-0.5,@(obj,t)(1),10^(-12));
T = 1;
n = [16,32,64,128,256,512,1024,2048];
kappa = 1;
N = 10000;

%% Excl. random number generation:
% T, n_total, iter, runtime
nRun = 100;
[res_he,res_tbss,res_tbss_conv2] = deal([]);
temp = NaN(1,4);
m = NaN(size(n));

% Incl. random number generation
for i=1:size(n,2)
    i
    % Random numbers:
    M = sum((0:1:floor(n(i)*T)+1)/n(i) <= T);
    Z = randn((M-1)*N,kappa+1);

    for j=1:nRun
        j

        tic;
        Y = HybridTBSSScheme(N,n(i),T,1,H-0.5,kappa,Z,'fft',[],[],[]);
        tt = toc;
        temp(1) = n(i);
        temp(2) = floor(T*n(i));
        temp(3) = j;
        temp(4) = tt*1000;
        res_tbss = [res_tbss;temp];
        
        tic;
        Y = HybridTBSSScheme(N,n(i),T,1,H-0.5,kappa,Z,'conv_loop',[],[],[]);
        tt = toc;
        temp(1) = n(i);
        temp(2) = floor(T*n(i));
        temp(3) = j;
        temp(4) = tt*1000;
        res_tbss_conv2 = [res_tbss_conv2;temp];   
        
          
        tic;
        % Kernel approximation:
        [c,gamm] = ApproximateKernel('BM2005','H',H,'T',T,'epsilon',0.001,...
                                     'error_measure','l2_normalised',...
                                     'options',[],'delta',1/n(i),'N',n(i));
        c = c.*gamma(H + 0.5);
        m(i) = size(c,1);        
        Xpaths = HybridMultifactorScheme(N,n(i),T,0,0,1,gamm,c,kappa,'K',K,'Z',Z,'explicit',false);
        tt = toc;
        temp(1) = n(i);
        temp(2) = floor(T*n(i));
        temp(3) = j;
        temp(4) = tt*1000;
        res_he = [res_he;temp];        
        
    end
        
end

res_he = array2table(res_he);
res_he.Properties.VariableNames = {'n','n_total','iter','runtime'};

res_tbss = array2table(res_tbss);
res_tbss.Properties.VariableNames = {'n','n_total','iter','runtime'};

res_tbss_conv2 = array2table(res_tbss_conv2);
res_tbss_conv2.Properties.VariableNames = {'n','n_total','iter','runtime'};

stat_he = grpstats(res_he,{'n_total'},{'mean','std'},'DataVars','runtime');
stat_tbss = grpstats(res_tbss,{'n_total'},{'mean','std'},'DataVars','runtime');
stat_tbss_conv2 = grpstats(res_tbss_conv2,{'n_total'},{'mean','std'},'DataVars','runtime');

stat_he.se = stat_he.std_runtime ./ sqrt(nRun);
stat_tbss.se = stat_tbss.std_runtime ./ sqrt(nRun);
stat_tbss_conv2.se = stat_tbss_conv2.std_runtime ./ sqrt(nRun*ones(size(stat_tbss_conv2.std_runtime)));

round([stat_he.mean_runtime,stat_tbss.mean_runtime,stat_tbss_conv2.mean_runtime])
[stat_he.se,stat_tbss.se,stat_tbss_conv2.se]./...
[stat_he.mean_runtime,stat_tbss.mean_runtime,stat_tbss_conv2.mean_runtime]
round([stat_he.se,stat_tbss.se,stat_tbss_conv2.se])
[stat_tbss.mean_runtime./stat_he.mean_runtime,stat_tbss_conv2.mean_runtime./stat_he.mean_runtime]

m

%% Save results:
save([fileparts(matlab.desktop.editor.getActiveFilename),'\comp_times.mat'],...
      'stat_he','stat_tbss','stat_tbss_conv2');

%% Conclusion:
% We conclude that overall the hybrid multifactor scheme is faster than the hybrid TBSS scheme.
% The difference is minor for a moderate number of steps but notably faster at higher number of 
% steps. This is despite both having O(nlog(n)) complexity in the singular case. We explain the
% difference with the fact that practical values of m as demonstrated in section x are 
% comparatively low. Based on figure x we can conclude that the scheme would perform even better
% relative to the hybrid TBSS approach for non-singular kernels (the number of computations for
% the H-TBSS scheme is unchanged for different K's). Finally, as demonstrated by the O(n^2) running
% time of the discrete convolutions we can conclude that the hybrid multifactor scheme is many
% times faster than the H-TBSS method when the diffusion coefficient is state-dependent.
% Overall this favours the hybrid mulitfactor approach.
% 

% TODO: Rerun the last one... (nightly run).
