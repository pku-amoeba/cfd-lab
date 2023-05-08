function [alph,ii]=get_mesh_normal(ny,Ly,dy1,alph0)

% solve expanding ratio
ny1=(ny-1)/2+1;
fun1 = @(alph) exp(alph/(ny1-1))-1-2.0*dy1/Ly*(exp(alph)-1);
%alph0 = 8;
alph = fsolve(fun1,alph0);

% generate mesh
for j=1:ny1
    gg(j)=0.5*Ly*(exp(alph*(j-1)/(ny1-1))-1)/(exp(alph)-1);
end

for j=1:ny1
    hh(j)=-0.5*Ly*(exp(alph*(j-1)/(ny1-1))-1)/(exp(alph)-1);
end

ii(1:ny1)=hh(end:-1:1);
ii(ny1:ny)=gg;

end