function [ ] = plot_1u1E_versus_x( u, ti, E, simulation )
% 画图显示 1电势/1电场-网格节点编号
yyaxis left
plot(u,'-b');%某时刻电压值
axis([0,simulation.num_grid_point,-inf,inf])
%         title('电势空间分布', 'FontSize', 18);
xlabel('网格节点编号')
ylabel('\phi [V]');
hold off
yyaxis right
hold off
plot(E,'-.','Color',[0.8500    0.3250    0.0980]);%某时刻电场值
ylabel('E [V/m]');
legend_str1=['当前, 第' num2str(ti) '时步'];
L1=legend(['\phi-' legend_str1],'E-当前');
set(L1,'location','best');
hold off;
end