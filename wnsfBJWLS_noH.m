function [G,Vmin,it_end,it_best,norm_end] = wnsfBJWLS_noH(eta,R,mf,ml,it,Tol,u,y,H,n_max)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(y);
n = length(eta)/2;
a = eta(1:n);
b = eta(n+1:end);
Qf = toeplitz([0; b(1:end-1)],zeros(1,mf));
Ql = toeplitz([1; a(1:end-1)],[1 zeros(1,ml-1)]);
Q = [-Qf Ql];
theta = Q\b;

V = Inf*ones(it,1);
sys_all = cell(it,1);
theta_prev = theta;
for i=1:it
    Tf = toeplitz([1; theta(1:mf); zeros(n-mf-1,1)],[1 zeros(1,n-1)]);
    Tl = toeplitz([0; theta(mf+1:mf+ml); zeros(n-ml-1,1)],zeros(1,n));
    T = [-Tl Tf];
    [~,R2] = qr(R'\T',0);
    QR = Q'/R2;
    theta = (QR')\((R2'\b));
    F = tf([1 theta(1:mf)'],1,1,'variable','z^-1');
    L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
    sys_all{i,1} = L/F;
    if isstable(L/F)==1
        pe = lsim(1/H,y-lsim(L/F,u));
        V(i) = 1/(N-n_max)*sum(pe(n_max+1:end).^2);
    else
        V(i) = Inf;
    end
    if i==it || norm(theta-theta_prev)/norm(theta)<Tol
        [Vmin,Vmin_i] = min(V);
        G = sys_all{Vmin_i,1};
        it_end = i;
        it_best = Vmin_i;
        norm_end = norm(theta-theta_prev)/norm(theta);
        return
    end
    theta_prev = theta;
end

end

