%
% PLOT visualize the data
%
% Created    : Apr. 15, 2020; (c) Chen Junlin
% Last update: Apr. 15, 2020
%
clear
% ========== set as specific ========== %
A=importdata('wave_packet_upwind.dat');
x=A.data(:,1);
u=A.data(:,2);

B=importdata('wave_packet_ue.dat');
ue=B.data(:,2);
% ========== set as specific ========== %

set(gcf,'unit','centimeters','position',[10 5 7 5]);
plot(x,u,'ro', 'MarkerSize',2)
ylabel('u')
xlabel('x')

hold on
plot(x,ue)

% ========== set as specific ========== %
%set(gca,'XTick',-2*pi:0.1:1);
%xticklabels({'-1','','','','','-0.5','','','','','0','','','','','0.5','','','','','1'}); 
legend({'UPWIND','Exact'},'FontSize',5,'Location','NorthEast')
saveas(gcf, 'wave_packet_upwind_4s', 'png')
% ========== set as specific ========== %

