%
% Case#02 Riemann problem for 1-D Euler equations with a moving Mach=3 shock
% interactiong with sine waves in density.
%
% Initial data, in original variables: 
% OVL=(3.85714,2.62936,10.33333), x<-4;
% OVR=(1+0.2sin[5x],0,1)         x>=-4.  
%
% Numerical schemes including Lax-Friedrichs and two-step Lax-Wendroff are applied;
%
% Created    : Apr. 30, 2020; (c) Chen Junlin
% Last update: Apr. 30, 2020
%
clear  
%% user parameters
% ========== set as appropriate ========== %
Xmin=-5;        % spatial interval (start point)
Xmax=5;         % spatial interval (end point)
Tmax=1.8;         % end of time
N=400;          % number of nodes
cfl=0.95;       % time step factor
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
% original variables
rho1=3.85714;
RHO4=1+0.2.*sin(5.*X);
u1=2.62936;
u4=0;
p1=10.33333;
p4=1;
% ========== set as specific ========== %
e1=p1./(gamma-1)+0.5.*rho1.*u1.*u1;
E4=p4./(gamma-1)+0.5.*RHO4.*u4.*u4;

CV(1,:)=(X<-4).*rho1+(X>=-4).*RHO4;
CV(2,:)=(X<-4).*rho1.*u1+(X>=-4).*RHO4.*u4;
CV(3,:)=(X<-4).*e1+(X>=-4).*E4;
% conservative variables 

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
    % Update convective variables 
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

% hold on
% plot(X,RHOe)
legend({flag},'FontSize',5,'Location','SouthWest')
saveas(gcf, ['case2_',flag,'_rho'], 'png')

figure(2)
% Plot velocity
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,U,'ro', 'MarkerSize',2)
ylabel('Velocity')
xlabel('x')

% hold on
% plot(X,Ue)
saveas(gcf, ['case2_',flag,'_u'], 'png')

figure(3)
% Plot pressure
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,P,'ro', 'MarkerSize',2)
ylabel('Pressure')
xlabel('x')
% 
% hold on
% plot(X,Pe)
legend({flag},'FontSize',5,'Location','SouthWest')
saveas(gcf, ['case2_',flag,'_p'], 'png')
% 
figure(4)
% Plot entropy
set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(X,E,'ro', 'MarkerSize',2)
ylabel('Entropy')
xlabel('x')
% 
% hold on
% plot(X,Ee)
saveas(gcf, ['case2_',flag,'_e'], 'png')