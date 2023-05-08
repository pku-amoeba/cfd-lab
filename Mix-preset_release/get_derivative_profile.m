function [Du,DT,Drho]=get_derivative_profile(u,T,y,U2,T2,Pr,Ma,gamma)

Du=zeros(size(u));
for i=2:size(u,2)-1
    Du(i)=(u(i+1)-u(i-1))/(y(i+1)-y(i-1));
end
Du(1)=Du(2);
Du(end)=Du(end-1);

DT=(1-T2).*Du./(1-U2)+1/2*sqrt(Pr)*(gamma-1)*Ma^2.*(-2*u+U2+1).*Du;
Drho=-1./T.^2.*DT;

end