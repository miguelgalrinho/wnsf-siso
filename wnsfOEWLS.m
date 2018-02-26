function [G,Vmin,it_end,it_best,norm_end] = wnsfOEWLS(eta,R,mf,ml,it,Tol,u,y,n_max,theta0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(y);
P = (R'*R);

if n_max>=0
    n = length(eta);
    Q = [-toeplitz([0; eta(1:end-1)],zeros(1,mf)) toeplitz([1; zeros(n-1,1)],[1 zeros(1,ml-1)])];
    if isempty(theta0)==1
        Ti = eye(n);
        ii = 0;
    else
        F = tf([1 theta0(1:mf)'],1,1,'variable','z^-1');
        fi = impulse(1/F,n-1);
        Ti = toeplitz(fi,[1 zeros(1,n-1)]);
        ii = 1;
        theta_prev = theta0;
    end
    V = Inf*ones(it,1);
    sys_all = cell(it,1);
    for i=ii:it
%         if i==0
%             theta = Q\eta;
%         else
        W = Ti'*P*Ti;
        theta = (Q'*W*Q)\(Q'*W*eta);
%             LTiQ = Lc*Ti*Q;
%             theta = LTiQ\(Lc*Ti*eta);
%         end
        F = tf([1 theta(1:mf)'],1,1,'variable','z^-1');
        if i>0
            L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
            sys_all{i,1} = L/F;
            pe = y-lsim(L/F,u);
            V(i) = 1/N*sum(pe.^2);
            if i==it || norm(theta-theta_prev)/norm(theta)<Tol || sum(isnan(theta))>0 || sum(isinf(theta))>0
                [Vmin,Vmin_i] = min(V);
                G = sys_all{Vmin_i,1};
                it_end = i;
                it_best = Vmin_i;
                norm_end = norm(theta-theta_prev)/norm(theta);
                return
            end  
        elseif it==0
            L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
            pe = y-lsim(L/F,u);
            Vmin = 1/N*sum(pe.^2);
            G = L/F;
            it_end = 0;
            it_best = 0;
            norm_end = NaN;
        end
        fi = impulse(1/F,n-1);
        Ti = toeplitz(fi,[fi(1) zeros(1,n-1)]);
        theta_prev = theta;
    end
else
    n = length(eta)/2;
    Q = [-toeplitz([0; eta(1:n-1)],zeros(1,mf)) toeplitz([1; zeros(n-1,1)],[1 zeros(1,ml-1)]) zeros(n,ml);
         -toeplitz([0; eta(n+1:end-1)],zeros(1,mf)) zeros(n,ml) toeplitz([1; zeros(n-1,1)],[1 zeros(1,ml-1)])];
    if isempty(theta0)==1
        Ti = eye(2*n);
        ii = 0;
    else
        F = tf([1 theta0(1:mf)'],1,1,'variable','z^-1');
        fi = impulse(1/F,n-1);
        Ti_aux = toeplitz(fi,[1 zeros(1,n-1)]);
        Ti = blkdiag(Ti_aux,Ti_aux);
        ii = 1;
        theta_prev = theta0;
    end
    V = Inf*ones(it,1);
    sys_all = cell(it,2);
    for i=ii:it
%         if i==0
%             theta = Q\eta;
%         else
            W = Ti'*P*Ti;
            theta = (Q'*W*Q)\(Q'*W*eta);
%         end
        F = tf([1 theta(1:mf)'],1,1,'variable','z^-1');
        if i>0
            L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
            L2 = tf([0 theta(mf+ml+1:end)'],1,1,'variable','z^-1');
            sys_all{i,1} = L/F;
            sys_all{i,2} = L2/F;
            pe = y-lsim(L/F,u)-lsim(L2/F,[1; zeros(N-1,1)]);
            V(i) = 1/N*sum(pe.^2);
            if i==it || norm(theta-theta_prev)/norm(theta)<Tol || sum(isnan(theta))>0 || sum(isinf(theta))>0
                [Vmin,Vmin_i] = min(V);
                G = [sys_all{Vmin_i,1} sys_all{Vmin_i,2}];
                it_end = i;
                it_best = Vmin_i;
                norm_end = norm(theta-theta_prev)/norm(theta);
                return
            end  
        elseif it==0
            L = tf([0 theta(mf+1:mf+ml)'],1,1,'variable','z^-1');
            L2 = tf([0 theta(mf+ml+1:end)'],1,1,'variable','z^-1');
            pe = y-lsim(L/F,u)-lsim(L2/F,[1; zeros(N-1,1)]);
            Vmin = 1/N*sum(pe.^2);
            G = [L/F L2/F];
            it_end = 0;
            it_best = 0;
            norm_end = NaN;
        end
        fi = impulse(1/F,n-1);
        Ti_aux = toeplitz(fi,[fi(1) zeros(1,n-1)]);
        Ti = blkdiag(Ti_aux,Ti_aux);
        theta_prev = theta;
    end
end

end

