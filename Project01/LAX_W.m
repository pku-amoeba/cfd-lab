function [uo] = LAX_W(uin,f,a,r)
% LAX_W update the 1-D convection equation with the Lax-Wendroff scheme 
% in conservative form.
%
% [uo] = LAX_W(uin,f,a,r) update the 1-D convection equation with the 
% Lax-Wendroff scheme. uo is the data at time (n+1), uin is the data at time (n), 
% f defines the flux funtion, a is calculated from f by derivation of u, 
% r is the (time step/spatial step) satisfying CFL condition. 
%
% Created    : Apr. 11, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020

upp=circshift(uin,-1);      % u_{j+1} with periodic B.C.
unn=circshift(uin,1);       % u_{j-1} with periodic B.C.
up=1/2*(uin+upp);
un=1/2*(uin+unn);

% update the data 
uo=uin-r/2.*(f(upp)-f(unn))+1/2.*r^2.*a(up).^2.*(upp-uin)-1/2.*r^2.*a(un).^2.*(uin-unn);

end