%
% Case#01 Traditional Riemann Problem for the 1-D Euler equations 
% Initial data, in conservative variables: 
% CVL=(1,0,2.5),      x<0.3; 
% CVR=(0.125,0,0.25), x>=0.3.
% Numerical schemes including Lax-Friedrichs and two-step Lax-Wendroff are applied;
%
% Created    : Apr. 24, 2020; (c) Chen Junlin
% Last update: Apr. 30, 2020
%
clear  
%% user parameters
% ========== set as appropriate ========== %
Xmin=0;         % spatial interval (start point)
Xmax=1;         % spatial interval (end point)
Tmax=0.2;       % end of time
N=100;          % number of nodes
cfl=0.95;        % time step factor
scheme=1;       % numerical scheme (1)Lax-Friedrichs (2)(Two-step)Lax-Wendroff

% ========== set as appropriate ========== %

%% comput.domain
dx=(Xmax-Xmin)/N; 
X=Xmin:dx:Xmax;
t=0;

%% I.C.
% ========== set as specific ========== %
% constant of the gas property
gamma=1.4;      
% conservative variables
CV=(X<0.3).*[1;0;2.5]+(X>=0.3).*[0.125;0;0.25];
% ========== set as specific ========== %

% original variables (prepared for input of the exact solution)
rho1=CV(1,1);
rho4=CV(1,end);
u1=CV(2,1)./CV(1,1);
u4=CV(2,end)./CV(1,end);
p1=(gamma-1)*(CV(3,1)-0.5*CV(2,1).*CV(2,1)./CV(1,1));
p4=(gamma-1)*(CV(3,end)-0.5*CV(2,end).*CV(2,end)./CV(1,end));

%% solver
while t<Tmax
    dt=Tstep(CV,dx,gamma,cfl);
    t=t+dt;
    % Calculate flux
    if scheme==1
            RHS=Flux_LF(CV,dt/dx,gamma);
    elseif scheme==2
            RHS=Flux_LW(CV,dt/dx,gamma);
    end
    % Update conservative variables 
    CV=CV-dt/dx*RHS;
    % Boundary condition
    CV(:,1)=CV(:,2);
    CV(:,end)=CV(:,end-1);
    
    % Original variables, entropy and pressure
    RHO=CV(1,:);
    U=CV(2,:)./CV(1,:);
    P=(gamma-1)*(CV(3,:)-0.5*CV(2,:).*CV(2,:)./CV(1,:));
    E=CV(3,:);
end

%% Exact solution  
tol=1e-6;
% shift the discont' point
xshift=0.3;
[Ue,Pe,Ae]=RiemannExact(p1,rho1,u1,p4,rho4,u4,tol,t,Xmin-xshift,Xmax-xshift,N+1);
RHOe=gamma*Pe./Ae./Ae;
Ee=Pe/(gamma-1)+0.5.*RHOe.*Ue.*Ue;

%% Plot data
if scheme==1
    flag='Lax-Friedrichs';
elseif scheme==2
    flag='MacCormack';
end

figure(1)
% Plot density
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,RHO,'ro', 'MarkerSize',2)
ylabel('Density')
xlabel('x')

hold on
plot(X,RHOe)
legend({flag,'Exact'},'FontSize',5,'Location','NorthEast')
saveas(gcf, ['case1_',flag,'_rho','_CFL',num2str(cfl*100)], 'png')

figure(2)
% Plot velocity
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,U,'ro', 'MarkerSize',2)
ylabel('Velocity')
xlabel('x')

hold on
plot(X,Ue)
saveas(gcf, ['case1_',flag,'_u','_CFL',num2str(cfl*100)], 'png')

figure(3)
% Plot pressure
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,P,'ro', 'MarkerSize',2)
ylabel('Pressure')
xlabel('x')

hold on
plot(X,Pe)
legend({flag,'Exact'},'FontSize',5,'Location','NorthEast')
saveas(gcf, ['case1_',flag,'_p','_CFL',num2str(cfl*100)], 'png')

figure(4)
% Plot entropy
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,E,'ro', 'MarkerSize',2)
ylabel('Entropy')
xlabel('x')

hold on
plot(X,Ee)
saveas(gcf, ['case1_',flag,'_e','_CFL',num2str(cfl*100)], 'png')