function [G,H,Vmin,it_end,it_best,norm_end] = wnsfBJWLS(eta,R,Lc,mf,ml,mc,md,it,Tol,u,y,n_max)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(y);
n = length(eta)/2;
a = [1; eta(1:n-1)];
b = [0; eta(n+1:end-1)];
Qf = toeplitz(b,zeros(1,mf));
Ql = toeplitz(a,[1 zeros(1,ml-1)]);
Qc = toeplitz(a,[1 zeros(1,mc-1)]);
Qd = [eye(md); zeros(n-md,md)];
Q = [zeros(n,mf) zeros(n,ml) -Qc         Qd;
     -Qf         Ql          zeros(n,mc) zeros(n,md)];
theta = Q\eta;
F = tf([1 theta(1:mf)'],1,1,'variable','z^-1');
L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
C = tf([1 theta(mf+ml+1:mf+ml+mc)'],1,1,'variable','z^-1');
l_fc = impulse(L/(F*C),n-1);
fi = impulse(1/F,n-1);
ci = impulse(1/C,n-1);

if it==0
    D = tf([1 theta(mf+ml+mc+1:end)'],1,1,'variable','z^-1');
    G = L/F;
    H = C/D;
    it_end = NaN;
    it_best = NaN;
    norm_end = NaN;
    if isstable(D/C) && isstable(L/F) && sum(isnan(theta))==0
        pe = lsim(D/C,y-lsim(L/F,u));
        Vmin = 1/(N-n_max)*sum(pe(n_max+1:end).^2);
    else
        Vmin = Inf;
    end
    return
end

V = Inf*ones(it,1);
sys_all = cell(it,2);
theta_prev = theta;
for i=1:it
    Ti11 = toeplitz(ci,[ci(1) zeros(1,n-1)]);
    Ti22 = toeplitz(fi,[fi(1) zeros(1,n-1)]);
    Ti21 = toeplitz(l_fc,[l_fc(1) zeros(1,n-1)]);
    Ti = [Ti11 zeros(n,n);
          Ti21 Ti22];
%     W = Ti'*(R'*R)*Ti;
%     theta = (Q'*W*Q)\(Q'*W*eta);
    LTiQ = Lc*Ti*Q;
    theta = LTiQ\(Lc*Ti*eta);
    F = tf([1 theta(1:mf)'],1,1,'variable','z^-1');
    L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
    C = tf([1 theta(mf+ml+1:mf+ml+mc)'],1,1,'variable','z^-1');
    D = tf([1 theta(mf+ml+mc+1:end)'],1,1,'variable','z^-1');
    sys_all{i,1} = L/F;
    sys_all{i,2} = C/D;
    if isstable(D/C) && isstable(L/F) && sum(isnan(theta))==0
        pe = lsim(D/C,y-lsim(L/F,u));
        V(i) = 1/(N-n_max)*sum(pe(n_max+1:end).^2);
    else
        V(i) = Inf;
    end
    if i==it || norm(theta-theta_prev)/norm(theta)<Tol || sum(isnan(theta))>0
        [Vmin,Vmin_i] = min(V);
        G = sys_all{Vmin_i,1};
        H = sys_all{Vmin_i,2};
        it_end = i;
        it_best = Vmin_i;
        norm_end = norm(theta-theta_prev)/norm(theta);
        return
    end  
    l_fc = impulse(L/(F*C),n-1);
    fi = impulse(1/F,n-1);
    ci = impulse(1/C,n-1);
    theta_prev = theta;
end

end

