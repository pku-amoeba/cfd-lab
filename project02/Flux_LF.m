function [rhs] = Flux_LF(U,r,gamma)
% Flux_LF calculates the net flux of 1D Euler equations (finite difference based) 
% with the Lax-Friedrichs scheme.
%
% [rhs] = Flux_LF calculates the net flux of 1D Euler equations.
% rhs outputs the net flux around local points, U is convective factor, 
% r is time/space ratio, gamma is gas property constant.
%
% Created    : Apr. 24, 2020; (c) Chen Junlin
% Last update: Apr. 30, 2020
F=zeros(size(U));

P=(gamma-1)*(U(3,:)-0.5*U(2,:).*U(2,:)./U(1,:));

F(1,:)=U(2,:);
F(2,:)=U(2,:).*U(2,:)./U(1,:)+P;
F(3,:)=U(2,:)./U(1,:).*(U(3,:)+P);

rhs=0.5*(circshift(F,[0,-1])-circshift(F,[0,1]))-0.5/r*(circshift(U,[0,-1])-2*U+circshift(U,[0,1]));

end