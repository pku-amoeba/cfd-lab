function [rhs] = Flux_LW(U,r,gamma)
% Flux_LW calculates the net flux of 1D Euler equations (finite difference based) 
% with the two-step Lax-Wendroff scheme , namely, MacCormack's scheme.
%
% [rhs] = Flux_LW calculates the net flux of 1D Euler equations.
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

UM=U-r*(circshift(F,[0,-1])-F);
PM=(gamma-1)*(UM(3,:)-0.5*UM(2,:).*UM(2,:)./UM(1,:));
FM(1,:)=UM(2,:);
FM(2,:)=UM(2,:).*UM(2,:)./UM(1,:)+PM;
FM(3,:)=UM(2,:)./UM(1,:).*(UM(3,:)+PM);

rhs=0.5*(circshift(F,[0,-1])-F+FM-circshift(FM,[0,1]));

end