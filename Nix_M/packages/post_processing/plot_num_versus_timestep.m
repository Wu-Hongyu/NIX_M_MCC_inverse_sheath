function [ ] = plot_num_versus_timestep( e_num,Hp_num,H2p_num,H3p_num,dptime )
% ��ͼ��ʾ ������Ŀ-ʱ�䲽
plot(e_num,'-b')
hold on
plot(Hp_num,'-r')
plot(H2p_num,'-g')
plot(H3p_num,'-m')
xlabel(['ʱ�� [' num2str(dptime) 'ʱ��]']);
ylabel('������Ŀ');
% title('������Ŀʱ���ݻ�');
legend('e','Hp','H2p','H3p')
hold off
end