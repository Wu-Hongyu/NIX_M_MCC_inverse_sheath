function [ ] = plot_no_versus_position( position_e,position_H,simulation )
% 画图显示 粒子编号-粒子位置
scatter(position_H,1:size(position_H),1,'r')%显示离子位置
%         title('粒子空间分布');
hold on
if(~isempty(position_e))
scatter(position_e,1:size(position_e),1,'b')%显示电子位置
end
axis([0 simulation.Lx 0 simulation.num0_macro_e]);
hold on

xlabel('x [m]');
ylabel('粒子编号');
if(~isempty(position_e))
legend('H','e')
else
    legend('H')
end
hold off
end