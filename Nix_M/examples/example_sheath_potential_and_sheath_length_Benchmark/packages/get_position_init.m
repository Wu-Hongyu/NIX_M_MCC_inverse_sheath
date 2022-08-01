function [ position_init ] = get_position_init( simulation, type_init, dimension_vec )
% 初始化粒子位置
assert(1==dimension_vec(2)) % 目前仅适用于1D
switch type_init
    case 'entire domain uniform+noise'
        %均匀分布在空间Lx内，并给一个随机扰动
        num=dimension_vec(1);
        delta=simulation.Lx/num;
        position_init=(0:delta:simulation.Lx-delta)'+5e-7*rand(num,1);
    case 'entire domain uniform random'
        position_init=simulation.Lx*rand(dimension_vec);
    case 'source region uniform random' 
        source_region_Lx=simulation.source_region(2)-simulation.source_region(1);
        position_init=simulation.source_region(1)+source_region_Lx*rand(dimension_vec);
end
end