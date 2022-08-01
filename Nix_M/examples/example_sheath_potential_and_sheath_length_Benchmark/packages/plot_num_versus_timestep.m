function [ ] = plot_num_versus_timestep( e_num,H_num,H1n_num,dptime )
% 画图显示 粒子数目-时间步
plot(e_num,'-b')
hold on
plot(H_num,'--r')
hold on
plot(H1n_num,'--k')
xlabel(['时间 [' num2str(dptime) '时步]']);
ylabel('粒子数目');
% title('粒子数目时间演化');
legend('e','H','H-')
hold off
end