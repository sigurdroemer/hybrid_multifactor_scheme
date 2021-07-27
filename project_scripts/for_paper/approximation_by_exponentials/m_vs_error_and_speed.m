%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));

%% Experiment:
% Run those where 'm' is an input:
delta = 1/500;
H = 0.10;
n = 500;
T = 1;
m_test = (1:35)';
tt = (delta:delta:T)';
options = optimset('Display','off','MaxFunEval',10^5,'MaxIter',10^5);
n_run = 100;

% Order is: m, runtime, se, error
[res_JaberOptim,res_JaberEuch,res_JaberEuchOptim,res_Optim] = deal(NaN(size(m_test,1),3));
for i=1:size(m_test,1)
    i
    
    % JaberOptim:
    runtime = NaN(n_run,1);
    for j=1:n_run
        tic;
        [c,gamm,K,K_hat,exitflag] = ApproximateKernel('AJ2019optim','H',H,'T',T,'m',m_test(i),...
                                              'options',[],'delta',delta,'options',options); 
        if exitflag < 1
            error('Problem');
        end         
        runtime(j) = toc;
    end
    
    res_JaberOptim(i,1) = size(gamm,1);
    res_JaberOptim(i,2) = mean(runtime);
    res_JaberOptim(i,3) = std(runtime) / n_run;
    res_JaberOptim(i,4) =  sqrt(sum((K(tt)-K_hat(tt)).^2))/sqrt(sum(K(tt).^2));
    
    % JaberEuch
    runtime = NaN(n_run,1);
    for j=1:n_run    
        tic;
        [c,gamm,K,K_hat] = ApproximateKernel('AJEE2019','H',H,'T',T,...
                                             'm',m_test(i),'options',[],'delta',delta);  
        runtime(j) = toc;
    end
    
    res_JaberEuch(i,1) = size(gamm,1);
    res_JaberEuch(i,2) = mean(runtime);
    res_JaberEuch(i,3) = std(runtime)/n_run;
    res_JaberEuch(i,4) =  sqrt(sum((K(tt)-K_hat(tt)).^2))/sqrt(sum(K(tt).^2));
    
    % JaberEuch
    runtime = NaN(n_run,1);
    for j=1:n_run     
        tic;
        [c,gamm,K,K_hat,exitflag] = ApproximateKernel('AJEE2019optim','H',H,'T',T,...
                                         'm',m_test(i),'options',[],'delta',delta,...
                                         'options',options);   
        if exitflag < 1
            error('Problem');
        end                                         
        runtime(j) = toc;
    end
    
    res_JaberEuchOptim(i,1) = size(gamm,1);
    res_JaberEuchOptim(i,2) = mean(runtime);
    res_JaberEuchOptim(i,3) = std(runtime) / n_run;
    res_JaberEuchOptim(i,4) =  sqrt(sum((K(tt)-K_hat(tt)).^2))/sqrt(sum(K(tt).^2));
    
    runtime = NaN(n_run,1);
    for j=1:n_run
        tic;
        % Get initial guess first:
        [c0,gamm0,K,K_hat,exitflag] = ApproximateKernel('AJ2019optim','H',H,'T',T,'m',m_test(i),...
                                               'options',[],'delta',delta,'options',options);
        if exitflag < 1
            error('Problem');
        end
        [c,gamm,K,K_hat,exitflag] = ApproximateKernel('l2optim','H',H,'T',T,'m',m_test(i),...
                                             'options',options,'delta',delta,'n',n,...
                                              'gamm0',gamm0,'options',options);    
                                          
        if exitflag < 1
            error('Problem');
        end
        runtime(j) = toc;
    end
    res_Optim(i,1) = size(gamm,1);
    res_Optim(i,2) = mean(runtime);
    res_Optim(i,3) = std(runtime) / n_run;
    res_Optim(i,4) =  sqrt(sum((K(tt)-K_hat(tt)).^2))/sqrt(sum(K(tt).^2));
    
end


res_beylkin = NaN(size(m_test,1),4);
options = optimset('Display','off','MaxFunEval',10^6,'MaxIter',10^6);
for i=1:size(m_test,1)
    i
    % Beylkin:
    runtime = NaN(n_run,1);
    for j=1:n_run    
    tic;
    [c,gamm,K,K_hat] = ApproximateKernel('BM2005','H',H,'T',T,'m',i,...
                                         'options',options,'delta',delta,'n',n); 
        runtime(j) = toc;
    end
    res_beylkin(i,1) = size(gamm,1);
    res_beylkin(i,2) = mean(runtime);
    res_beylkin(i,3) = std(runtime) / n_run;
    res_beylkin(i,4) =  sqrt(sum((K(tt)-K_hat(tt)).^2))/sqrt(sum(K(tt).^2));
          
end

%% Save results:
save([fileparts(matlab.desktop.editor.getActiveFilename),'\error_runtime_data.mat'],...
    'res_JaberEuch','res_JaberEuchOptim','res_JaberOptim','res_Optim','res_beylkin');

load([fileparts(matlab.desktop.editor.getActiveFilename),'\error_runtime_data.mat']);

%% Check m = 20 case for AJ2019-optim: 
% Debug to get rm value!
[c,gamm,K,K_hat,exitflag] = ApproximateKernel('AJ2019optim','H',H,'T',T,...
                                 'm',20,'options',[],'delta',delta,...
                                 'options',options);

%% Plot (m, error)
figure;
subplot(1,2,1)  
lw = 1;
cm = lines(6);

plot(res_JaberEuch(:,1),res_JaberEuch(:,4),'-','LineWidth',lw,'Color',cm(1,:));hold on;
plot(res_JaberEuchOptim(:,1),res_JaberEuchOptim(:,4),'-','LineWidth',lw,'Color',cm(2,:));hold on;
plot(res_JaberOptim(:,1),res_JaberOptim(:,4),'-','LineWidth',lw,'Color',cm(3,:));hold on;
plot(res_Optim(:,1),res_Optim(:,4),'-x','LineWidth',lw,'Color',cm(4,:));hold on;
plot(res_beylkin(1:15,1),res_beylkin(1:15,4),'-x','LineWidth',lw,'Color',cm(5,:));
hleg = legend('AJEE2018','AJEE2018-optim',...
              'AJ2019-optim',...
              '$l^2$-optim','BM2005','Interpreter','Latex');
title(hleg,'Method');
xlabel('$m$','Interpreter','Latex');%,'FontSize',13);
ylabel('$l^2$-error (normalised)','Interpreter','Latex');%,'FontSize',13);
set(gca, 'YScale', 'log')


%% Plot (m, speed)
subplot(1,2,2);
plot(res_JaberEuch(:,1),res_JaberEuch(:,2),'-x','Color',cm(1,:),'LineWidth',lw);hold on;
plot(res_JaberEuchOptim(:,1),res_JaberEuchOptim(:,2),'-x','Color',cm(2,:),'LineWidth',lw);hold on;
plot(res_JaberOptim(:,1),res_JaberOptim(:,2),'-x','Color',cm(3,:),'LineWidth',lw);hold on;
plot(res_Optim(:,1),res_Optim(:,2),'-x','Color',cm(4,:),'LineWidth',lw);hold on;
plot(res_beylkin(1:15,1),res_beylkin(1:15,2),'-x','Color',cm(5,:),'LineWidth',lw);
hleg = legend('AJEE2019','AJEE2019-optim','AJ2019-optim',...
              '$l^2$-optim','BM2005','Interpreter','Latex');
title(hleg,'Method');%,'FontSize',10);
xlabel('$m$','Interpreter','Latex');%,'FontSize',13);
ylabel('Running time (in seconds)','Interpreter','Latex');%,'FontSize',13,
set(gca, 'YScale', 'log')


%% Standard errors:
([res_JaberEuch(:,3),res_JaberEuchOptim(:,3),res_JaberOptim(:,3),res_Optim(:,3),res_beylkin(:,3)]./...
[res_JaberEuch(:,2),res_JaberEuchOptim(:,2),res_JaberOptim(:,2),res_Optim(:,2),res_beylkin(:,2)])

