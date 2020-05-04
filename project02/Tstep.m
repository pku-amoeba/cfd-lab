function [dt]=Tstep(CV,dx,gamma,cfl)
% Tstep calculates the local time-step, which satisfies the CFL condition
%
% Created    : Apr. 30, 2020; (c) Chen Junlin
% Last update: Apr. 30, 2020
RHO=CV(1,:);
U=CV(2,:)./CV(1,:);
P=(gamma-1)*(CV(3,:)-0.5*CV(2,:).*CV(2,:)./CV(1,:));

CS=sqrt(gamma.*P./RHO);
SPRAD=CS+abs(U);

maxvel=max(SPRAD);
dt=cfl*dx/maxvel;

end