function [ simulation ] = get_simulation( type )
% simulation 仿真参数 结构体
switch type
    case 'default' %A1DPIC_Ver1.0.m, 2019Montellano
        simulation.dt=10*10^-12;%时间步长
        simulation.end_time=6E-6;%仿真时间总长
        simulation.all_timesteps=floor(simulation.end_time/simulation.dt);%总循环次数
        simulation.Lx=0.01;%仿真区域长度
        simulation.source_region=[0.1*simulation.Lx,0.4*simulation.Lx];
        simulation.num_grid_point=201;%网格格点数目
        simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
        simulation.n0=1E16;%等离子体密度
        simulation.num0_macro_e=20000;%初始电子数目 Particle e Num
        simulation.num0_macro_H=simulation.num0_macro_e;%初始离子数目Particle H Num
        simulation.weight=simulation.n0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
        simulation.Te=1;%单位eV
        simulation.TH=1;%单位eV
        simulation.field_boundaries_type=0;%电位边界条件类型
        simulation.field_boundaries=[0,0];%电位边界条件值
        simulation.pressure=0.6;%气压 单位Pa
        simulation.Tgas=300;%气温   K
    case 'BACON base'
        error('Not Done');
end

