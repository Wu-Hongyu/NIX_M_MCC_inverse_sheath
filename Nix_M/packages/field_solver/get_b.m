function [ b ] = get_b( rho,b_extra,simulation )
% 组装右端向量
b=rho;
switch simulation.field_boundaries_type
    case 0 %00，首末第一类边界
        b(1)=b_extra(1);
        b(end)=b_extra(end);
    case 1 %01，首第一类，末第二类
        b(1)=b_extra(1);
        b(end)=b(end)+b_extra(end);
    case 2 %10，首第二类，末第一类
        b(1)=b(1)+b_extra(1);
        b(end)=b_extra(end);
    case 3 %11，周期边界
        b(1)=b_extra(1);
        b(end)=b_extra(end);
end
end