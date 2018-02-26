clear;

N = 10000; % Sample size
G = tf([0 1],[1 -.9],1,'variable','z^-1'); % System
H = tf([1 .7],[1 -.6],1,'variable','z^-1'); % Noise Model

nl = 1; % order of numerator of G, one delay assumed; if not one delay, shift the data
nf = 1; % order of denominator of G
nc = 1; % order of numerator of H
nd = 1; % order of numerator of H

% Signals
e = randn(N,1); % Noise
u = filter(1,[1 -.8],randn(N,1)); % Input

% OE case
y = lsim(G,u)+e; % Output
z = iddata(y,u);
model_wnsf_oe = wnsfOE(z,[nf nl],'UnstructuredOrders',[20 40 60]);

% BJ case
y = lsim(G,u)+lsim(H,e); % Output
z = iddata(y,u);
model_wnsf_bj = wnsfBJ(z,[nf nl nc nd],'UnstructuredOrders',[20 40 60]); % parametric H
model_wnsf_npH = wnsfBJ(z,[nf nl],'UnstructuredOrders',[20 40 60]); % non-parametric H
