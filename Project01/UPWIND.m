function [uo] = UPWIND(uin,f,a,r)
% UPWIND update the 1-D convection equation with the upwind scheme 
% at order one, in conservative form.
%
% [uo] = UPWIND(uin,f,a,r) update the 1-D convection equation with the 
% upwind scheme. uo is the data at time (n+1), uin is the data at time (n), 
% f defines the flux funtion, a is calculated from f by derivation of u, 
% r is the (time step/spatial step) satisfying CFL condition. 
%
% Created    : Apr. 11, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020


upp=circshift(uin,-1);      % u_{j+1} with periodic B.C.
unn=circshift(uin,1);       % u_{j-1} with periodic B.C
up=1/2*(uin+upp);           
un=1/2*(uin+unn);           

% update the data 
uo=uin-r/2.*(f(upp)-f(unn))+1/2.*r.*abs(a(up)).*(upp-uin)-1/2.*r.*abs(a(un)).*(uin-unn);

end