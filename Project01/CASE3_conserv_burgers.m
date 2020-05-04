%
% Case#03 
% Numerically solve the I.V.P. of the Burgers equation in conservative form
% with periodic boundaries 
%
% Created    : Apr. 11, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020
%
clear 
%% user parameters
% ========== set as appropriate ========== %
Xmin=0;             % spatial interval (start point)
Xmax=2*pi;          % spatial interval (end point)
Tmax=2.5;           % end of time
N=500;              % number of nodes
CFL=1;              % CFL condition
ss=3;               % numerical scheme (1)CIR (2)Lax-Friedrichs (3)Lax-Wendroff

ff=@(u) 1/2.*u.^2;  % flux function
aa=@(u) u;          % df(u)/du

T=0.9;              % Print data at time step T
% ========== set as appropriate ========== %

%% comput.domain
dx=(Xmax-Xmin)/N; 
x=Xmin:dx:Xmax;
t=0;

%% I.C.
uu=0.5+sin(x);
u=uu;

%% solver
switch ss
    case 1        
        while t<Tmax
            dt=dx*CFL/max(aa(u));
            t=t+dt;
            u=UPWIND(u,ff,aa,dt/dx);
            if abs(t-T)<dt
                fid=fopen('burgers_upwind.dat','w');
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
                fid=fopen('burgers_lf.dat','w');
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
                fid=fopen('burgers_lw.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%14.6f %14.6f\r\n',[x;u]);
                fclose(fid);
             end 
        end
end 

%% Exact solution 
i=1;
for xx=Xmin:dx:Xmax
        syms u1
        eq=u1-(0.5+sin(xx-u1*T));
        ue(i)=double(solve(eq,u1));
        i=i+1;
end
fid=fopen('burgers_ue.dat','w');
fprintf(fid,'%7s %9s\r\n','x','ue');
fprintf(fid,'%14.6f %14.6f\r\n',[x;ue]);
fclose(fid);