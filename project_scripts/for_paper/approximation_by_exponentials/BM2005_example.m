%% Initialize Script
clear;
serverRun = false;
project_folder = fileparts(fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename))));
addpath(genpath(project_folder));

%% Set-up:
close all;
epsilon = 10^(-2);
delta = 1/500;
b = 1;
T = b;
N = 250;
H = 0.1;
K = @(t)(t.^(H - 0.5) ./ gamma(H + 0.5));
f = @(x)(K(x.*(b-delta)+delta));
ii = (0:2*N).';
x = ii./(2*N);
h = f(x);
H = hankel(h(1:N+1),h(N+1:end));

% Compute eigenvalues:
L = N+1;
[V, Lambda, flag] = eigs(H,min(size(H,1),L));
[sigmas,ind] = sort(max(diag(Lambda),0),'descend');

figure;
plot((0:19)',sigmas(1:20),'.','MarkerSize',8.5);
set(gca, 'YScale', 'log');
xlabel('Index');
ylabel('Eigenvalue');

% Compute roots and weights:
%m = find(sigmas./sigmas(1) < mu,1);
m_idx = find(sigmas./sqrt(sum(h.^2)) <= epsilon,1);
m = m_idx - 1;
u = flipud(V(:,ind(m_idx)));
gamm_ = roots(u);

%A = (gamm_.').^ii;
%w = pinv(A)*h;

figure;
plot(gamm_,'.');
xlim([-1.15,1.15]);
ylim([-1.15,1.15]);
xlabel('Real part');
ylabel('Imaginary part');

idxKeep = imag(gamm_) == 0 & real(gamm_) >= 0 & real(gamm_) <= 1;
gamm_ = gamm_(idxKeep);
gamm_ = gamm_(1:min(m,size(gamm_,1)));
A = (gamm_.').^ii;
w = pinv(A)*h;

% Transform gamma's back:
t = 2*N*log(gamm_);

% Transform solution to K on [delta,b]:
gamm = -t./(b-delta);
c = w.*exp(gamm.*delta);

imag([c,gamm])
flipud(round(real([gamm,c]),2))

% Compute uniform error:
K_hat = @(t)(sum(exp(-gamm.*t).*c,1));
N_intervals = 2*N;
delta_sample = (T-delta)/N_intervals;
tt = (delta:delta_sample:T)';

[K_hat(tt')',K(tt)]

% Relative error:
sqrt(sum((K_hat(tt')'-K(tt)).^2))./sqrt(sum(K(tt)))
max(abs(K_hat(tt')'-K(tt))./K(tt))

size(c,1)

round(flipud([gamm,c]),2)




