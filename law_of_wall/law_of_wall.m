set(0,'defaultlinelinewidth',3)
set(0,'defaultaxeslinewidth',3);
set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'DefaultLineMarkerSize',12);

yp=0:0.1:1e3;
k=0.41;B=5.2;
lglaw=1/k*log(yp)+B;

lm1=k*yp;
f1=2./(1+(1+4*lm1.^2).^0.5);
up1=zeros(size(yp));
for i=1:length(yp)
   up1(i)=trapz(yp(1:i),f1(1:i));
end

lm2=k*yp.*(1-exp(-yp/26));
f2=2./(1+(1+4*lm2.^2).^0.5);
up2=zeros(size(yp));
for i=1:length(yp)
   up2(i)=trapz(yp(1:i),f2(1:i));
end

figure()
% (1) Prandtl two layer
plot(yp,yp.*(yp<11)+lglaw.*(yp>=11),'-.');hold on;
% (2) Prandtl one layer
plot(yp,up1,'--')
% (3) von Karman three layer
plot(yp,up2,'o','MarkerSize',10,'MarkerIndices',1:10:length(yp))
% (4) van Driest damping function
plot(yp,up2)

set(gca, 'XScale', 'log')
xlim([0.1 1000])
ylim([0 25])
xlabel('$y^+$','Interpreter','latex')
ylabel('$u^+$','Interpreter','latex')
legend('Prandtl two layer','Prandtl one layer','Karman three layer','van Driest')
