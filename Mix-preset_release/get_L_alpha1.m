function L1=get_L_alpha1(N,u,rho,T,M1,gamma)

L1=zeros(5*N,5*N);
% Matrix for alpha^1
for k=2:N-1
    L1(1+5*(k-1),1+5*(k-1))=u(k);
    L1(1+5*(k-1),2+5*(k-1))=rho(k);
    L1(2+5*(k-1),1+5*(k-1))=T(k)/rho(k)/gamma/M1^2;
    L1(2+5*(k-1),2+5*(k-1))=u(k);
    L1(2+5*(k-1),5+5*(k-1))=1/gamma/M1^2;
    L1(3+5*(k-1),3+5*(k-1))=u(k);
    L1(4+5*(k-1),4+5*(k-1))=u(k);
    L1(5+5*(k-1),2+5*(k-1))=(gamma-1)/rho(k);
    L1(5+5*(k-1),5+5*(k-1))=u(k);
end

end