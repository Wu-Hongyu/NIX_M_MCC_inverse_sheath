function [ b ] = get_b( rho,b_extra,simulation )
% ��װ�Ҷ�����
b=rho;
switch simulation.field_boundaries_type
    case 0 %00����ĩ��һ��߽�
        b(1)=b_extra(1);
        b(end)=b_extra(end);
    case 1 %01���׵�һ�࣬ĩ�ڶ���
        b(1)=b_extra(1);
        b(end)=b(end)+b_extra(end);
    case 2 %10���׵ڶ��࣬ĩ��һ��
        b(1)=b(1)+b_extra(1);
        b(end)=b_extra(end);
    case 3 %11�����ڱ߽�
        b(1)=b_extra(1);
        b(end)=b_extra(end);
end
end