function [ ] = plot_no_versus_position( position_e,position_H,simulation )
% ��ͼ��ʾ ���ӱ��-����λ��
scatter(position_H,1:size(position_H),1,'r')%��ʾ����λ��
%         title('���ӿռ�ֲ�');
hold on
if(~isempty(position_e))
scatter(position_e,1:size(position_e),1,'b')%��ʾ����λ��
end
axis([0 simulation.Lx 0 simulation.num0_macro_e]);
hold on

xlabel('x [m]');
ylabel('���ӱ��');
if(~isempty(position_e))
legend('H','e')
else
    legend('H')
end
hold off
end