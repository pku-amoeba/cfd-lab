% Version 2021.01.27: Reduce 6 equations into 5
% Version 2022.02.21: Revise L0(4,1)
% Version 2023.04.18: Use function, auto generating mesh, Blasius profile
%%
clear
ny=251; nyc=ny-1;% N is number of nodes, including upper and lower boundaries
nx=2600;nx1=2500;nz=128;Lx=400;Lx1=350;Lz=35;beta0=5;
Ly=200;
dy1=0.1;aa0=8;
M1=2.8; beta=0; omg=0.41; %beta=-0.3601;
%-----Information of basic flow----------
Ma=M1;
U1=1;U2=0.5;R=(1-U2)/(1+U2);
T2=1;gamma=1.4;Pr=0.72;
%-----Get mesh---------- 
[aa,y]=get_mesh_normal(ny,Ly,dy1,aa0);
yc=get_mesh_center(y,ny);
%-----Get basic flow profile----------
u=get_basic_profile(0.5,yc,200,0.005,1000,1e-5);
T=get_T_profile(u,U2,T2,Pr,M1,gamma);
rho=get_rho_profile(T);
[Du,DT,Drho]=get_derivative_profile(u,T,yc,U2,T2,Pr,Ma,gamma);
%-----Main calculation
% Linear stability matrix
L1=get_L_alpha1(nyc,u,rho,T,M1,gamma);
L0=get_L_alpha0(nyc,yc,omg,beta,rho,T,Du,Drho,DT,M1,gamma);
% Eigenvalue & vector
[eigfun,alph] = eig(L0(1+5:5+5*(nyc-2),1+5:5+5*(nyc-2)),...
                    L1(1+5:5+5*(nyc-2),1+5:5+5*(nyc-2)),'chol','vector');
%%
% Plot discrete spectrum and find the 
flag=0;
for ii=1:length(alph)
    if abs(real(alph(ii))-0.5368)<0.001 && abs(imag(alph(ii))+0.002486)<0.001
        flag=ii;
    end
end
%temp=0+0i;
%for ii=1:length(alph)
%    if real(alph(ii))>0 && imag(alph(ii))<0 && abs(real(alph(ii)))<2 && imag(alph(ii))<imag(temp) 
%        %plot(real(alph(ii)),imag(alph(ii)),'*');hold on
%        temp=alph(ii);
%        flag=ii;
%    end
%end
%alpha_max(jj)=temp;
figure()
plot(real(alph),imag(alph),'*')
hold on
plot(real(alph(flag)),imag(alph(flag)),'o')

%% Plot
figure()
plot(yc(2:end-1),real(eigfun(2:5:end,flag)),'--')
hold on;
plot(yc(2:end-1),imag(eigfun(2:5:end,flag)),'-.')
hold on;
%plot(abs(eigfun(2:5:end,flag)),y(2:end-1))
legend('Real','Imag')
xlabel('$y$','Interpreter','latex')
ylabel('$\hat u$','Interpreter','latex')
% 
ph=rho(2:end-1)'.*eigfun(5:5:end,flag) + eigfun(1:5:end,flag).*T(2:end-1)';
% 
figure()
plot(yc(2:end-1), real(ph),'-')
hold on;
plot(yc(2:end-1), imag(ph),'-.')
hold on;
legend('Real','Imag')
xlabel('$y$','Interpreter','latex')
ylabel('$\hat p$','Interpreter','latex')

%% Convert data to CFD code
output_basic_profile(u,yc)
output_eigenfun(nyc,eigfun,ph,flag)
output_mesh(nx,nx1,ny,nz,Lx,Lx1,Ly,Lz,dy1,aa0,beta0)