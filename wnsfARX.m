function [eta,R,L] = wnsfARX(u,y,n,InitialConditions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmpi(InitialConditions,'zero')==1
    Phi = [-toeplitz([0; y(1:end-1)],zeros(1,n)) toeplitz([0; u(1:end-1)],zeros(1,n))];
    eta = Phi\y;
else
    Phi = [-toeplitz(y(n:end-1),fliplr(y(1:n)')) toeplitz(u(n:end-1),fliplr(u(1:n)'))];
    eta = Phi\y(n+1:end);
end
[~,R] = qr(Phi,0);
[~,Ri]=qr(fliplr(Phi),0);
L=rot90(Ri,2);


end

