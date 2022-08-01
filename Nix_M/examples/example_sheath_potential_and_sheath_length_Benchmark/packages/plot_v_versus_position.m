function [ ] = plot_v_versus_position( ve, vH, position_e,position_H,simulation )
% 画图显示 相空间
scatter(position_e,ve(:,1),1,'b')%显示电子位置
axis([0 simulation.Lx -inf inf]);
hold on
scatter(position_H,vH(:,1),1,'r')%显示离子位置
%         title('粒子相空间分布');
xlabel('x [m]');
ylabel('v [m/s]');
legend('e','H')
hold off
end