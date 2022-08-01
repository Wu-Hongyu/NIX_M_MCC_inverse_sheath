function [ ] = plot_no_versus_position( position_e,position_H,position_H2p,position_H3p,simulation )
% ��ͼ��ʾ ���ӱ��-����λ��

% �Զ��ı�marker��С
num_macro=max([length(position_e) length(position_H) length(position_H2p)]);
if num_macro>200
    scatter(position_e,1:length(position_e),1,'b')%��ʾ����λ��
else
    size_marker=floor(200/num_macro);
    scatter(position_e,1:length(position_e),size_marker,'b')%��ʾ����λ��
end
axis([0 simulation.Lx 0 num_macro]);
hold on
scatter(position_H,1:length(position_H),1,'r')%��ʾ����λ��
scatter(position_H2p,1:length(position_H2p),1,'g')%��ʾ����λ��
scatter(position_H3p,1:length(position_H3p),1,'m')%��ʾ����λ��
%         title('���ӿռ�ֲ�');
xlabel('x [m]');
ylabel('���ӱ��');
legend({'e','H+','H2+'},'AutoUpdate','off')
hold off
end