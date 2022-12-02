clear all;close all
load example

set(0,'DefaultAxesFontSize',12) 
set(0,'defaultlinelinewidth',2)
set(0,'defaultlinemarkersize',8)
figure;
 ha = tight_subplot(3,1,[.06 .06],[.08 .06],[.07 .07]);
 set(gcf,'Units','normalized');
 set(gcf,'Position',[0 0 0.7 0.8]);
 axes(ha(1))
semilogy(ssh_joint,cdf_joint_2050,ssh_joint,cdf_joint_2100,ssh_joint,cdf_joint_2150);axis tight
title('sea level [mean+extreme]')
legend('2020-2050','2020-2100','2020-2150','location','best')
%xlabel('sea level maximum within planning period [m]')
ylabel('frequency')
 axes(ha(2))
semilogy(ssh_mean,cdf_mean_only_2050,ssh_mean,cdf_mean_only_2100,ssh_mean,cdf_mean_only_2150);axis tight
title('sea level [mean part]')
%legend('2020-2050','2020-2100','2020-2150')
%xlabel('sea level maximum within planning period [m]')
ylabel('frequency')
 axes(ha(3))
 semilogy(ssh_ext,cdf_ext_only_2050,ssh_ext,cdf_ext_only_2100,ssh_ext,cdf_ext_only_2150);axis tight
title('sea level [extreme part]')
%legend('2020-2050','2020-2100','2020-2150')
xlabel('sea level maximum within planning period [m]')
ylabel('frequency')

figure;
 ha = tight_subplot(3,1,[.09 .06],[.08 .06],[.07 .07]);
 set(gcf,'Units','normalized');
 set(gcf,'Position',[0 0 0.7 0.8]);
 axes(ha(1))
plot(ssh_joint,scen_ssh(:,3),'o',ssh_joint,mean_quant_ssh(:,3),'+',ssh_joint,ext_quant_ssh(:,3),'p');axis tight
title('factors controlling sea level maximum [2020-2050]')
legend('mean projection','mean quantile','GEV quantile','location','best')
xlabel('sea level maximum within planning period [m]')
ylabel('average quantile')
axes(ha(2))
plot(ssh_joint,scen_ssh(:,8),'o',ssh_joint,mean_quant_ssh(:,8),'+',ssh_joint,ext_quant_ssh(:,8),'p');axis tight
title('factors controlling sea level maximum [2020-2100]')
%legend('mean projection','mean quantile','GEV quantile')
xlabel('sea level maximum within planning period [m]')
ylabel('average quantile')
axes(ha(3))
plot(ssh_joint,scen_ssh(:,13),'o',ssh_joint,mean_quant_ssh(:,13),'+',ssh_joint,ext_quant_ssh(:,13),'p');axis tight
title('factors controlling sea level maximum [2020-2150]')
%legend('mean projection','mean quantile','GEV quantile')
xlabel('sea level maximum within planning period [m]')
ylabel('average quantile')

figure;
 ha = tight_subplot(2,2,[.09 .08],[.08 .08],[.07 .07]);
 set(gcf,'Units','normalized');
 set(gcf,'Position',[0 0 0.7 0.8]);
 axes(ha(1))
 surf(ssh_joint,ssh_ext,ssh_reached_ext(:,:,13)'./(nr_iter*nr_par));view(2);colorbar;hold on
p1=plot(ssh_joint,ext_of_tot(:,13),'r');
p1.ZData=ones(1,nr_res)*1000;
p1.LineWidth=3;
title('extreme contribution to joint sea level maximum [2020-2150]')
set(ha(1),'Units','normalized')
titleHandle = get( ha(1) ,'Title' );
pos  = get( titleHandle , 'position' );
pos1 = pos + [0 0.2 0]; 
set( titleHandle , 'position' , pos1 );
xlabel('sea level mean+extreme [m]')
ylabel('sea level extreme contribution [m]')
axes(ha(2))
surf(ssh_joint,ssh_mean,ssh_reached_mean(:,:,13)'./(nr_iter*nr_par));view(2);colorbar;hold on
p1=plot(ssh_joint,mean_of_tot(:,13),'r');
p1.ZData=ones(1,nr_res)*1000;
p1.LineWidth=3;
title('    mean contribution to joint sea level maximum [2020-2150]')
set(ha(2),'Units','normalized')
titleHandle = get( ha(2) ,'Title' );
pos  = get( titleHandle , 'position' );
pos1 = pos + [0 .85 0]; 
set( titleHandle , 'position' , pos1 );
xlabel('sea level mean+extreme [m]')
ylabel('sea level mean contribution [m]')
axes(ha(4))
plot(ssh_joint,mean_of_tot(:,3),ssh_joint,mean_of_tot(:,8),ssh_joint,mean_of_tot(:,13));axis tight
title('sea level [mean part]')
legend('2020-2050','2020-2100','2020-2150','location','best')
xlabel('sea level maximum within planning period [m]')
ylabel('average mean contribution [m]')
axes(ha(3))
plot(ssh_joint,ext_of_tot(:,3),ssh_joint,ext_of_tot(:,8),ssh_joint,ext_of_tot(:,13));axis tight
title('sea level [extreme part]')
legend('2020-2050','2020-2100','2020-2150','location','best')
xlabel('sea level maximum within planning period [m]')
ylabel('average extreme contribution [m]')
