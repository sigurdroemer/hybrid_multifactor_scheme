%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));

%% Plot actual fits:
T = 1;
delta = 1/500;
H = 0.1;
t = (delta:delta:T)';
K = @(t)(t.^(H - 0.5) ./ gamma(H + 0.5));

%% (Jaber and Euch, 2019)
m_test = [25;50;100;200];
cm = lines(size(m_test,1)+1);
cm = cm([1,3:end],:);
figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJEE2019','H',H,'T',T,...
                                         'm',m_test(i),'options',[]);
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJEE2019','Interpreter','Latex','FontSize',11);
hold off;
xlim([-0.02,0.25]);
ylim([0.5,9]);
legend boxoff;

figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJEE2019','H',H,'T',T,...
                                         'm',m_test(i),'options',[]);
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJEE2019 (log-log plot)','Interpreter','Latex','FontSize',11);
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ylim([0.5,11]);
legend boxoff;

%% (Jaber and Euch, 2019) Optimize
m_test = [25;50;100;200];
cm = lines(size(m_test,1)+1);
cm = cm([1,3:end],:);
options = optimset('Display','iter','MaxFunEval',10^5,'MaxIter',10^5);
figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJEE2019optim','H',H,'T',T,...
                                         'm',m_test(i),'options',[],...
                                         'options',options);                                        
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJEE2019-optim','Interpreter','Latex','FontSize',11);
hold off;
xlim([-0.02,0.25]);
ylim([0.5,9]);
legend boxoff;

figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJEE2019optim','H',H,'T',T,...
                                         'm',m_test(i),'options',[],...
                                         'options',options);                     
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10,'location','NorthEast');
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJEE2018optim (log-log plot)','Interpreter','Latex','FontSize',11);
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([10^(-3),1])
ylim([0.5,11]);
legend boxoff;

%% (Jaber, 2019) 
m_test = [5;10;15;20];
cm = lines(size(m_test,1)+1);
cm = cm([1,3:end],:);
figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJ2019','H',H,'T',T,'m',m_test(i),...
                                          'options',[],'rm',2.5);                               
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11);
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJ2019','Interpreter','Latex','FontSize',11);
hold off;
xlim([-0.02,0.25]);
ylim([0.5,9]);
legend boxoff;

figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJ2019','H',H,'T',T,'m',m_test(i),...
                                          'options',[],'rm',2.5);                    
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10,'location','SouthWest');
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJ2019 (log-log plot)','Interpreter','Latex','FontSize',11);
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([10^(-3),1])
ylim([0.5,11]);
legend boxoff;

%% (Jaber, 2019) Optimize
m_test = [3;5;10;20];
cm = lines(size(m_test,1)+1);
cm = cm([1,3:end],:);
options = optimset('Display','iter','MaxFunEval',10^5,'MaxIter',10^5);
figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJ2019optim','H',H,'T',T,...
                                         'm',m_test(i),'options',[],'delta',delta,...
                                         'options',options);                                        
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJ2019optim','Interpreter','Latex','FontSize',11);
hold off;
xlim([-0.02,0.25]);
ylim([0.5,9]);
legend boxoff;

figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('AJ2019optim','H',H,'T',T,...
                                         'm',m_test(i),'options',[],'delta',delta,...
                                         'options',options);                     
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10,'location','SouthWest');
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('AJ2019optim (log-log plot)','Interpreter','Latex','FontSize',11);
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([10^(-3),1])
ylim([0.5,11]);
legend boxoff;

%% Beylkin-Monzon 
m_test = (2:5)';
cm = lines(size(m_test,1)+1);
cm = cm([1,3:end],:);
figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('BM2005','H',H,'T',T,'m',m_test(i),...
                                         'delta',delta,'n',500);                                
    m = size(c,1);
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('BM2005','Interpreter','Latex','FontSize',11);
legend boxoff;
hold off;
xlim([-0.02,0.25]);
ylim([0.5,9]);

figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c,gamm,K,K_hat] = ApproximateKernel('BM2005','H',H,'T',T,'m',m_test(i),...
                                         'delta',delta,'n',500);                    
    m = size(c,1);
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10,'location','SouthWest');
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('BM2005 (log-log plot)','Interpreter','Latex','FontSize',11);
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend boxoff;
xlim([10^(-3),1])
ylim([0.5,11]);


%% l2-optimize
c0 = [];
m_test = [2;3;4;5];
cm = lines(size(m_test,1)+1);
cm = cm([1,3:end],:);
options = optimset('Display','off','MaxFunEval',10^5,'MaxIter',10^5);
figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c0,gamm0,K,K_hat,exitflag] = ApproximateKernel('AJ2019optim','H',H,'T',T,'m',m_test(i),...
                                                    'options',[],'delta',delta,'options',options);       
    
    [c,gamm,K,K_hat,exitflag] = ApproximateKernel('l2optim','H',H,'T',T,'m',m_test(i),...
                                         'options',options,'delta',delta,'n',500,...
                                         'gamm0',gamm0,'options',options);                                   
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
legend('Interpreter','Latex','FontSize',10);
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('$l^2$-optim','Interpreter','Latex','FontSize',13);
hold off;
xlim([-0.02,0.25]);
ylim([0.5,9]);
legend boxoff 

figure;
plot(t,K(t),'DisplayName','True kernel','Color','red','LineWidth',1.1);hold on;
for i=1:size(m_test,1)
    i
    [c0,gamm0,K,K_hat,exitflag] = ApproximateKernel('AJ2019optim','H',H,'T',T,'m',m_test(i),...
                                                    'options',[],'delta',delta,'options',options);    
    
    [c,gamm,K,K_hat,exitflag] = ApproximateKernel('l2optim','H',H,'T',T,'m',m_test(i),...
                                         'options',options,'delta',delta,'n',500,...
                                         'gamm0',gamm0,'options',options);              
    plot(t,K_hat(t),'--','DisplayName',['$m=',num2str(m_test(i)),'$'],'Color',cm(i,:),'LineWidth',1.1);hold on;
end
set(gca,'FontSize',11)
hleg = legend('Interpreter','Latex','FontSize',10,'location','SouthWest');
xlabel('$t$','Interpreter','Latex','FontSize',13);
ylabel('Values','Interpreter','Latex','FontSize',13);
title('$l^2$-optim (log-log plot)','Interpreter','Latex','FontSize',13);
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([10^(-3),1])
ylim([0.5,11]);
legend boxoff 



