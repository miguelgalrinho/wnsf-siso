function [eta,R] = wnsfFIR(u,y,n,InitialConditions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmpi(InitialConditions,'zero')==1
    Phi = toeplitz([0; u(1:end-1)],zeros(1,n));
    eta = Phi\y;
elseif strcmpi(InitialConditions,'truncate')==1
    Phi = toeplitz(u(n:end-1),fliplr(u(1:n)'));
    eta = Phi\y(n+1:end);
else
    N = length(y);
    Phi = [toeplitz([0; u(1:end-1)],zeros(1,n)) toeplitz([1; zeros(N-1,1)],[1 zeros(1,n-1)])];
    eta = Phi\y;
end
[~,R] = qr(Phi,0);

end

