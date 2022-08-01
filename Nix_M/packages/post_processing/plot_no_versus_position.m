function [ ] = plot_no_versus_position( position_e,position_H,position_H2p,position_H3p,simulation )
% 画图显示 粒子编号-粒子位置

% 自动改变marker大小
num_macro=max([length(position_e) length(position_H) length(position_H2p)]);
if num_macro>200
    scatter(position_e,1:length(position_e),1,'b')%显示电子位置
else
    size_marker=floor(200/num_macro);
    scatter(position_e,1:length(position_e),size_marker,'b')%显示电子位置
end
axis([0 simulation.Lx 0 num_macro]);
hold on
scatter(position_H,1:length(position_H),1,'r')%显示离子位置
scatter(position_H2p,1:length(position_H2p),1,'g')%显示离子位置
scatter(position_H3p,1:length(position_H3p),1,'m')%显示离子位置
%         title('粒子空间分布');
xlabel('x [m]');
ylabel('粒子编号');
legend({'e','H+','H2+'},'AutoUpdate','off')
hold off
end