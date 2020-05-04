function [uo] = LAX_F(uin,f,r)
% LAX_F update the 1-D convection equation with the Lax-Friedrichs scheme 
% in conservative form.
%
% [uo] = LAX_F(uin,f,a,r) update the 1-D convection equation with the 
% Lax-Friedrichs scheme. uo is the data at time (n+1), uin is the data at time (n), 
% f defines the flux funtion, r is the (time step/spatial step) 
% satisfying CFL condition. 
%
% Created    : Apr. 11, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020

upp=circshift(uin,-1);      % u_{j+1} with periodic B.C.
unn=circshift(uin,1);       % u_{j-1} with periodic B.C.

% update the data 
uo=1/2*(upp+unn)-r/2*(f(upp)-f(unn));

end