function [ ] = plot_v_versus_position( ve, vH, position_e,position_H,simulation )
% ��ͼ��ʾ ��ռ�
scatter(position_e,ve(:,1),1,'b')%��ʾ����λ��
axis([0 simulation.Lx -inf inf]);
hold on
scatter(position_H,vH(:,1),1,'r')%��ʾ����λ��
%         title('������ռ�ֲ�');
xlabel('x [m]');
ylabel('v [m/s]');
legend('e','H')
hold off
end