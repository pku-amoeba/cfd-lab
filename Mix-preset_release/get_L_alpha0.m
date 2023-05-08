function L0=get_L_alpha0(N,y,omg,beta,rho,T,Du,Drho,DT,M1,gamma)
% Matrix for alpha^0
L0=zeros(5*N,5*N);
for k=2:N-1
    % L0(1,1)
    L0(1+5*(k-1),1+5*(k-1))=omg;
    % L0(1,3)
    L0(1+5*(k-1),3+5*(k-1))=1i*Drho(k);
    L0(1+5*(k-1),3+5*(k))=1i*rho(k)/2/(y(k+1)-y(k));
    L0(1+5*(k-1),3+5*(k-2))=-1i*rho(k)/2/(y(k)-y(k-1));
    % L0(1,4)
    L0(1+5*(k-1),4+5*(k-1))=-rho(k)*beta;
    % L0(2,2)
    L0(2+5*(k-1),2+5*(k-1))=omg;
    % L0(2,3)
    L0(2+5*(k -1),3+5*(k-1))=1i*Du(k);
    % L0(3,1)
    L0(3+5*(k-1),1+5*(k-1))=1i*DT(k)/rho(k)/gamma/M1^2;
    L0(3+5*(k-1),1+5*(k))=1i*T(k)/rho(k)/gamma/M1^2/2/(y(k+1)-y(k));
    L0(3+5*(k-1),1+5*(k-2))=-1i*T(k)/rho(k)/gamma/M1^2/2/(y(k)-y(k-1));  
    % L0(3,3)
    L0(3+5*(k-1),3+5*(k-1))=omg;
    % L0(3,5)
    L0(3+5*(k-1),5+5*(k-1))=1i*Drho(k)/rho(k)/gamma/M1^2;
    L0(3+5*(k-1),5+5*(k))=1i/gamma/M1^2/2/(y(k+1)-y(k));
    L0(3+5*(k-1),5+5*(k-2))=-1i/gamma/M1^2/2/(y(k)-y(k-1));
    % L0(4,1)
    %L0(4+5*(k-1),1+5*(k-1))=-1i*beta*T(k)/rho(k)/gamma/M1^2;
    L0(4+5*(k-1),1+5*(k-1))=-beta*T(k)/rho(k)/gamma/M1^2;
    % L0(4,4)
    L0(4+5*(k-1),4+5*(k-1))=omg;
    % L0(4,5)
    L0(4+5*(k-1),5+5*(k-1))=-beta/gamma/M1^2;
    % L0(5,3)
    L0(5+5*(k-1),3+5*(k-1))=1i*DT(k);
    L0(5+5*(k-1),3+5*(k))=1i*(gamma-1)/rho(k)/2/(y(k+1)-y(k));
    L0(5+5*(k-1),3+5*(k-2))=-1i*(gamma-1)/rho(k)/2/(y(k)-y(k-1));   
    % L0(5,4)
    L0(5+5*(k-1),4+5*(k-1))=-(gamma-1)*beta/rho(k);
    % L0(5,5)
    L0(5+5*(k-1),5+5*(k-1))=omg;
end