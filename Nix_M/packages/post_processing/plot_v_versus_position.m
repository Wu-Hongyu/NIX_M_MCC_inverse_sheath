function [ ] = plot_v_versus_position( v, position,simulation ,type,color)
% ��ͼ��ʾ ��ռ�

scatter(position,v(:,1),1,color)%��ʾ����λ��
axis([0 simulation.Lx -inf inf]);
hold on
%         title('������ռ�ֲ�');
xlabel('x [m]');
ylabel('v [m/s]');
% legend(type)
hold off
end