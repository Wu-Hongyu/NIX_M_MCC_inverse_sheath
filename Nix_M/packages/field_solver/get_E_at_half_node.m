function [ E ] = get_E_at_half_node( u, simulation )
%�����ڵ�糡
E=zeros(simulation.num_grid_point+1,1);
E(2:end-1)=diff(-u)/simulation.dx;%����糡ֵ
switch simulation.field_boundaries_type
    case 0 %00����ĩ��һ��߽�
        E(1)=2*E(2)-E(3);%�������
        E(end)=2*E(end-1)-E(end-2);
    case 1 %01���׵�һ�࣬ĩ�ڶ���
        E(1)=2*E(2)-E(3);%�������
        E(end)=simulation.field_boundaries(2);
    case 2 %10���׵ڶ��࣬ĩ��һ��
        E(1)=simulation.field_boundaries(1);
        E(end)=2*E(end-1)-E(end-2);
    case 3 %11�����ڱ߽�����
        E(1)=E(end-1);
        E(end)=E(2);
end
end

