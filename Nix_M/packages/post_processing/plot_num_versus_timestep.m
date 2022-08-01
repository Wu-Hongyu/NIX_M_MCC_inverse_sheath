function [ ] = plot_num_versus_timestep( e_num,Hp_num,H2p_num,H3p_num,dptime )
% 画图显示 粒子数目-时间步
plot(e_num,'-b')
hold on
plot(Hp_num,'-r')
plot(H2p_num,'-g')
plot(H3p_num,'-m')
xlabel(['时间 [' num2str(dptime) '时步]']);
ylabel('粒子数目');
% title('粒子数目时间演化');
legend('e','Hp','H2p','H3p')
hold off
end