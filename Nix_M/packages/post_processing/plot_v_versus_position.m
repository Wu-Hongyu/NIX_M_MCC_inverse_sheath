function [ ] = plot_v_versus_position( v, position,simulation ,type,color)
% 画图显示 相空间

scatter(position,v(:,1),1,color)%显示电子位置
axis([0 simulation.Lx -inf inf]);
hold on
%         title('粒子相空间分布');
xlabel('x [m]');
ylabel('v [m/s]');
% legend(type)
hold off
end