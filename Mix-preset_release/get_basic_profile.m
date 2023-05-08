function u=get_basic_profile(ulow,y,eta,dt,K,tol)

% eta = 150;
% dt = 0.005;
N = 2*eta/dt+1;
t = linspace(-eta,eta,N);
% K = 1001;
% tol = 1e-5;

f = zeros(N,K);
g = zeros(N,K);
h = zeros(N,K);

ff = zeros(N,K);
gg = zeros(N,K);
hh = zeros(N,K);

% B.C. - middle
f((N-1)/2+1,:) = 0;    % F
% Initial guess for shooting method
g((N-1)/2+1,:) = ulow;  % F' 
h((N-1)/2+1,1) = 0.0875;    % F''
h((N-1)/2+1,2) = 0.1;

for k=1:K-1
% Shooting to +INF from 0
for i=(N-1)/2+2:N
    f(i,k)=f(i-1,k)+g(i-1,k)*dt;
    g(i,k)=g(i-1,k)+h(i-1,k)*dt;
    h(i,k) = h(i-1,k) - (f(i-1,k)*h(i-1,k)*dt)/2; 
end    
a=sqrt(g(N,k));
% Transform the equation, \yita -> \xi, F -> G
tt=t*a;
dtt=dt*a;
ff(:,k)=f(:,k)/a;
gg(:,k)=g(:,k)/a^2;
hh(:,k)=h(:,k)/a^3;
% Shooting to -INF from 0
for i=(N-1)/2:-1:1
    ff(i,k)=ff(i+1,k)-gg(i+1,k)*dtt;
    gg(i,k)=gg(i+1,k)-hh(i+1,k)*dtt;
    hh(i,k) = hh(i+1,k) + (ff(i+1,k)*hh(i+1,k)*dtt)/2; 
end 


% Compare G'(-INF) with U2
if abs(gg(1,k)-0.5) < tol
    iter = k;
    break;
end
% Choose new F''(0) and iterate
if k == 1 
    continue;
else
    hh((N-1)/2+1,k+1) = hh((N-1)/2+1,k) + (0.5-gg(1,k))*(hh((N-1)/2+1,k)-hh((N-1)/2+1,k-1))/(gg(1,k)-gg(1,k-1)); %>>>>> ??? 
        
    if sign(hh((N-1)/2+1,k+1)) == -1     % ???
        hh((N-1)/2+1,k+1) = -hh((N-1)/2+1,k+1);   % ???
    end                         % ???  >>>>>>
    
end

h((N-1)/2+1,k+1)=a^3*hh((N-1)/2+1,k+1);

end

figure()
plot((gg(:,iter)-ulow)/(1-ulow),tt)
hold on
load zhoudata0
plot(zhoudata0(:,1),zhoudata0(:,2),'o')

% Interpolation
u = interp1(tt,gg(:,iter),y,'spline');
figure()
plot(u)
figure()
plot(u,y,'o')
hold on
plot(gg(:,iter),tt)

end