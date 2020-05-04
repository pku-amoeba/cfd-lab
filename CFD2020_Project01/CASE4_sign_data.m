%
% Case#04
% Numerically solve the I.V.P. of a common conservative PDE with sign 
% initial data
%
% Created    : Apr. 11, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020
%
clear 
%% user parameters
% ========== set as appropriate ========== %
Xmin=-2;             % spatial interval (start point)
Xmax=2;          % spatial interval (end point)
Tmax=1.5;           % end of time
N=200;              % number of nodes
CFL=1;              % CFL condition
ss=3;               % numerical scheme (1)CIR (2)Lax-Friedrichs (3)Lax-Wendroff

ff=@(u) 1/4.*(u.^2-1).*(u.^2-4);        % flux function
aa=@(u) 1/2.*u.*(2.*u.^2-5);            % df(u)/du

T=1.2;              % Print data at time step T
% ========== set as appropriate ========== %

%% comput.domain
dx=(Xmax-Xmin)/N; 
x=Xmin:dx:Xmax;
t=0;

%% I.C.
t=0;
uu=-2*sign(x);
u=uu;
%% solver
switch ss
    case 1        
        while t<Tmax
            dt=dx*CFL/max(aa(u));
            t=t+dt;
            u=UPWIND(u,ff,aa,dt/dx);
            if abs(t-T)<dt
                fid=fopen('sign_upwind.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%11.3f %11.4f\r\n',[x;u]);
                fclose(fid);
            end
        end
        
    case 2
        while t<Tmax
            dt=dx*CFL/max(aa(u));            
            t=t+dt;
            u=LAX_F(u,ff,dt/dx);
            if abs(t-T)<dt
                fid=fopen('sign_lf.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%11.3f %11.4f\r\n',[x;u]);
                fclose(fid);
            end
        end
        
    case 3
        while t<Tmax
            dt=dx*CFL/max(aa(u));
            t=t+dt;
            u=LAX_W(u,ff,aa,dt/dx);
             if abs(t-T)<dt
                fid=fopen('sign_lw.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%14.6f %14.6f\r\n',[x;u]);
                fclose(fid);
             end 
        end
end 