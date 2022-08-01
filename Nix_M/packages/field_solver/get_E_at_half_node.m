function [ E ] = get_E_at_half_node( u, simulation )
%计算半节点电场
E=zeros(simulation.num_grid_point+1,1);
E(2:end-1)=diff(-u)/simulation.dx;%半格点电场值
switch simulation.field_boundaries_type
    case 0 %00，首末第一类边界
        E(1)=2*E(2)-E(3);%线性外插
        E(end)=2*E(end-1)-E(end-2);
    case 1 %01，首第一类，末第二类
        E(1)=2*E(2)-E(3);%线性外插
        E(end)=simulation.field_boundaries(2);
    case 2 %10，首第二类，末第一类
        E(1)=simulation.field_boundaries(1);
        E(end)=2*E(end-1)-E(end-2);
    case 3 %11，周期边界条件
        E(1)=E(end-1);
        E(end)=E(2);
end
end

