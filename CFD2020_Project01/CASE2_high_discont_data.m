%
% Case#02 Highly discontinuous data
% Numerically solve the I.V.P. of the linear advective equation with
% constant coefficient and highly discontinuous initial data.
%
% Created    : Apr. 10, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020
%
clear 
%% user parameters
% ========== set as appropriate ========== %
a=1;            % speed of wave (u_t + a*u_x = 0)
Xmin=-1;         % spatial interval (start point)
Xmax=1;         % spatial interval (end point)
Tmax=10;        % end of time
N=500;          % number of nodes
r=0.8;          % r=dt/dx 
ss=3;           % numerical scheme (1)CIR (2)Lax-Friedrichs (3)Lax-Wendroff

ff=@(u) a.*u;   % flux function
aa=@(u) a;      % df(u)/du

T=8;            % Print data at time step T
% ========== set as appropriate ========== %

%% comput.domain
dx=(Xmax-Xmin)/N; 
dt=r*dx;
x=Xmin:dx:Xmax;
t=0;

%% I.C.
xi=zeros(size(x));
uu=zeros(size(x));
for i=1:length(x)
    if (x(i)>=-0.7) && (x(i)<=1)
        xi(i)=x(i)-0.3;
    else
        xi(i)=x(i)-0.3+2;
    end
    
    if (xi(i)>=-1) && (xi(i)<-1/3)
        uu(i)=-xi(i)*sin(1.5*pi*xi(i)*xi(i));
    elseif (xi(i)>-1/3) && (xi(i)<1/3)
        uu(i)=abs(sin(2*pi*xi(i)));
    else
        uu(i)=2*xi(i)-1-sin(3*pi*xi(i))/6;
    end
end
u=uu;

%% solver
switch ss
    case 1        
        while t<Tmax
            t=t+dt;
            u=UPWIND(u,ff,aa,r);
            if abs(t-T)<0.00001
                fid=fopen('discont_upwind.dat','w');
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
                fid=fopen('discont_lf.dat','w');
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
                fid=fopen('discont_lw.dat','w');
                fprintf(fid,'%7s %9s\r\n','x','u');
                fprintf(fid,'%14.6f %14.6f\r\n',[x;u]);
                fclose(fid);
             end 
        end
end 

%% Exact solution 
ue=uu;
fid=fopen('discont_ue.dat','w');
fprintf(fid,'%7s %9s\r\n','x','u');
fprintf(fid,'%14.6f %14.6f\r\n',[x;ue]);
fclose(fid);
    

