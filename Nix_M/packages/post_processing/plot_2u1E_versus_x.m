function [ ] = plot_2u1E_versus_x( u, ti, usum, avgSteps,E, simulation )
% ��ͼ��ʾ 2����/1�糡-����ڵ���
yyaxis left
plot(u,'-b');%ĳʱ�̵�ѹֵ
axis([0,simulation.num_grid_point,-inf,inf])
hold on;
plot(usum/(avgSteps+1),'--r')%���������ѹƽ��ֵ
%         title('���ƿռ�ֲ�', 'FontSize', 18);
xlabel('����ڵ���')
ylabel('\phi [V]');
hold off
yyaxis right
hold off
plot(E,'-.','Color',[0.8500    0.3250    0.0980]);%ĳʱ�̵糡ֵ
ylabel('E [V/m]');
legend_str1=['��ǰ, ��' num2str(ti) 'ʱ��'];
legend_str2=['\phi-ǰ' num2str(avgSteps+1) 'ʱ��ƽ��'];
L1=legend(['\phi-' legend_str1],legend_str2,'E-��ǰ');
% Attention: 2020/11/05 01:23:30 ǰһ�汾��ͼ�����
set(L1,'location','south');
hold off;
end