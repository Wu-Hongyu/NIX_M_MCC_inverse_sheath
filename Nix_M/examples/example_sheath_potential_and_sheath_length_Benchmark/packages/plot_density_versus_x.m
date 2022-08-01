function [ ] = plot_density_versus_x( rhoe, rhoH, rho, simulation )
% 画图显示 粒子密度-网格节点编号
hold off
plot(-rhoe,'-b')
hold on
plot(rhoH,'--r')
plot(rho,'-.k')
source_region=floor(simulation.source_region/simulation.dx);
line([source_region(1),source_region(1)],[min(rho),max(rhoH)],'linestyle',':','color','k');
line([source_region(2),source_region(2)],[min(rho),max(rhoH)],'linestyle',':','color','k');
%         title('粒子密度分布', 'FontSize', 18);
ylabel('密度 [m^{-3}]');
xlabel('网格节点编号')
axis([0,simulation.num_grid_point,-inf,inf])
h1=plot(NaN,NaN,'-b');
h2=plot(NaN,NaN,'--r');
h3=plot(NaN,NaN,'-.k');
axes2 = axes('position',get(gca,'position'),'visible','off');
L1=legend(axes2,[h1,h2,h3],'e','H','净电荷');
set(L1,'location','northeast');
hold off
end