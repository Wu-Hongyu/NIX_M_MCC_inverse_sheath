function [ ] = plot_num_versus_timestep( e_num,H_num,H1n_num,dptime )
% ��ͼ��ʾ ������Ŀ-ʱ�䲽
plot(e_num,'-b')
hold on
plot(H_num,'--r')
hold on
plot(H1n_num,'--k')
xlabel(['ʱ�� [' num2str(dptime) 'ʱ��]']);
ylabel('������Ŀ');
% title('������Ŀʱ���ݻ�');
legend('e','H','H-')
hold off
end