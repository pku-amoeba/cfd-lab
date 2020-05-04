%
% Case#01 Propagation of a wave packet
% Calculation of the propagation of a wave packet. The initial data
% consists of a sine wave modulated by a Gaussian in 1-dimension. To obtain
% the pattens of the wave at a sequence of time, the linear advection
% equation with constant coefficient is solved. Numerical schemes including
% 1-D CIR(upwind). Lax-Friedrichs and Lax-Wendroff are applied£¬
%
% Created    : Apr. 10, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020
%
clear
%% user parameters
% ========== set as appropriate ========== %
a=1;            % speed of wave (u_t + a*u_x = 0)
Xmin=0;         % spatial interval (start point)
Xmax=1;         % spatial interval (end point)
Tmax=10;        % end of time
N=200;          % number of nodes
r=0.8;          % r=dt/dx 
ss=1;           % numerical scheme (1)CIR (2)Lax-Friedrichs (3)Lax-Wendroff

ff=@(u) a.*u;   % flux function
aa=@(u) a;      % df(u)/du

T=4;            % Print data at this time step
% ========== set as appropriate ========== %

%% comput.domain
dx=(Xmax-Xmin)/N; 
dt=r*dx;
x=Xmin:dx:Xmax;
t=0;

%% I.C.
u=exp(-100*(x-0.5).^2).*sin(80*x);

%% solver
switch ss
    case 1        
        while t<Tmax
            t=t+dt;
            u=UPWIND(u,ff,aa,r);
            if abs(t-T)<0.00001
                fid=fopen('wave_packet_upwind.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%11.3f %11.4f\r\n',[x;u]);
                fclose(fid);
            end
        end
        
    case 2
        while t<Tmax
            t=t+dt;
            u=LAX_F(u,ff,r);
            if abs(t-T)<0.00001
                fid=fopen('wave_packet_lf.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%11.3f %11.4f\r\n',[x;u]);
                fclose(fid);
            end
        end
        
    case 3
        while t<Tmax
            t=t+dt;
            u=LAX_W(u,ff,aa,r);
             if abs(t-T)<0.00001
                fid=fopen('wave_packet_lw.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%14.6f %14.6f\r\n',[x;u]);
                fclose(fid);
             end 
        end
end

%% Exact solution  
ue=exp(-100.*(x-0.5).^2).*sin(80*x);
fid=fopen('wave_packet_ue.dat','w');
fprintf(fid,'%7s %9s\r\n','x','u');
fprintf(fid,'%14.6f %14.6f\r\n',[x;ue]);
fclose(fid);
