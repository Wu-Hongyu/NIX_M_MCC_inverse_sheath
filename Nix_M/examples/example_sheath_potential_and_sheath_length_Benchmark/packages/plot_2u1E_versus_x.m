function [ ] = plot_2u1E_versus_x( u, ti, usum, avgSteps,E, simulation )
% 画图显示 2电势/1电场-网格节点编号
yyaxis left
plot(u,'-b');%某时刻电压值
axis([0,simulation.num_grid_point,-inf,inf])
hold on;
plot(usum/(avgSteps+1),'--r')%多个步长电压平均值
%         title('电势空间分布', 'FontSize', 18);
xlabel('网格节点编号')
ylabel('\phi [V]');
hold off
yyaxis right
hold off
plot(E,'-.','Color',[0.8500    0.3250    0.0980]);%某时刻电场值
ylabel('E [V/m]');
legend_str1=['当前, 第' num2str(ti) '时步'];
legend_str2=['\phi-前' num2str(avgSteps+1) '时步平均'];
L1=legend(['\phi-' legend_str1],legend_str2,'E-当前');
% Attention: 2020/11/05 01:23:30 前一版本的图例标错
set(L1,'location','south');
hold off;
end