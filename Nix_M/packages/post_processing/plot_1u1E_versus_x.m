function [ ] = plot_1u1E_versus_x( u, ti, E, simulation )
% ��ͼ��ʾ 1����/1�糡-����ڵ���
yyaxis left
plot(u,'-b');%ĳʱ�̵�ѹֵ
axis([0,simulation.num_grid_point,-inf,inf])
%         title('���ƿռ�ֲ�', 'FontSize', 18);
xlabel('����ڵ���')
ylabel('\phi [V]');
hold off
yyaxis right
hold off
plot(E,'-.','Color',[0.8500    0.3250    0.0980]);%ĳʱ�̵糡ֵ
ylabel('E [V/m]');
legend_str1=['��ǰ, ��' num2str(ti) 'ʱ��'];
L1=legend(['\phi-' legend_str1],'E-��ǰ');
set(L1,'location','best');
hold off;
end