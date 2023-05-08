function [beta,xx]=get_mesh_stream(nx,nx1,Lx,Lx1,beta0)

dx1=Lx1/(nx1-1);
Lx2=Lx-Lx1;
nx2=nx-nx1+1;

fun2 = @(beta) exp(beta/(nx2-1))-1-dx1/Lx2*(exp(beta)-1);
%beta0 = 5;
beta = fsolve(fun2,beta0);

for i=1:nx1
    xx(i)=(i-1)*dx1;
end

for i=1:nx2
    xx(nx1+i-1)=Lx1+Lx2*(exp(beta*(i-1)/(nx2-1))-1)/(exp(beta)-1);
end

plot(xx)

end