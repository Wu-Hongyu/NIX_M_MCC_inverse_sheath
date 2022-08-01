function [ simulation ] = get_simulation( type, varargin )
% 生成simulation结构体
% type: 'default', 'auto1'
% varargin: 可变输入参数，通过name-value对输入
% type='auto1': ne0, Te0, end_time, Lx
% simulation: 物理模型参数 结构体
%% 全部给定
switch type
    case 'default' %A1DPIC_Ver1.0.m, 2019Montellano
        %%%%%%% 常用，请勿修改 %%%%%%%%
        % 等离子体参数，以电子为代表
        simulation.ne0=1E16;%等离子体密度
        simulation.Te0=1;%单位eV
        %         TH=1;%单位eV
        % 时空尺度
        simulation.dt=10*10^-12;%时间步长
        simulation.end_time=6E-6;%仿真时间总长
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
        simulation.Lx=0.01;%仿真区域长度
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        simulation.num_grid_point=201;%网格格点数目
        simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
        % 电磁场
        simulation.field_boundaries_type=0;%电位边界条件类型
        simulation.field_boundaries=[0,0];%电位边界条件值
        simulation.B0=0; %代表性磁感应强度 T
        % 粒子
        simulation.particle_boundaries_type=0;%粒子边界条件类型
        %气压
        simulation.pressure=0.6;%气压 单位Pa
        simulation.Tgas=300;%气温   K
        case 'BATMAN_old_MF'
            % TODO：根据扩展腔模拟需要，做修改
        % 等离子体参数，以电子为代表
        simulation.ne0=1E18;%等离子体密度
        simulation.Te0=10;%单位eV
        % 时空尺度
        simulation.dt=10*10^-12;%时间步长
        simulation.end_time=6E-6;%仿真时间总长
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
        simulation.Lx=0.01;%仿真区域长度
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        simulation.num_grid_point=150;%网格格点数目
        simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
        % 电磁场
        simulation.field_boundaries_type=0;%电位边界条件类型
        simulation.field_boundaries=[0,0];%电位边界条件值
        simulation.B0=9e-3; %代表性磁感应强度 T
        % 粒子
        simulation.particle_boundaries_type=0;%粒子边界条件类型
         %气压
        simulation.pressure=0.6;%气压 单位Pa
        simulation.Tgas=293;%气温   K
        %     case 'BACON base'
        
    %% 自动计算
    %TODO：'auto000'这样的type。合并代码
    case 'auto0' %输入ne, Te, Lx, end_time，严格满足离散限制条件,可通过type3的check
        % 解析输入
        p = inputParser;
        validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
        validScalarNotNegNum = @(x) isnumeric(x) && isscalar(x) && ~(x < 0);
        addRequired(p,'type');
        addOptional(p,'ne0',1E16,validScalarPosNum);
        addOptional(p,'Te0',1,validScalarNotNegNum);
        addOptional(p,'end_time',6E-6,validScalarPosNum);
        addOptional(p,'Lx',0.01,validScalarPosNum);
        parse(p,type,varargin{:});
        % p.Results
        % 等离子体参数，以电子为代表
        simulation.ne0=p.Results.ne0;%等离子体密度
        simulation.Te0=p.Results.Te0;%单位eV
        simulation.end_time=p.Results.end_time;%仿真时间总长
        simulation.Lx=p.Results.Lx;%仿真区域长度
        % 时空尺度
        %%%%%%%%%% 时空步长
        % 考虑 蛙跳法求解简谐运动的精度要求 ωp*dt/2<<1
        omega_pe=get_omega_pe(simulation.ne0);
        simulation.dt=(2/omega_pe)*(pi/32); % dt=0.098*(2/ωp)=Tpe/32
        if simulation.Te0>0
            % 考虑数值自加热 dx<λD
            debye_length=get_debye_length(simulation.ne0,simulation.Te0);
            num_grid=ceil(simulation.Lx/(debye_length)); 
            simulation.num_grid_point=num_grid+1;%网格格点数目
            simulation.dx=simulation.Lx/num_grid;
            % 满足自加热限制条件时，CFL条件弱于蛙跳法精度要求，一般均能满足。
        else
            num_grid=20;
            simulation.num_grid_point=num_grid+1;%网格格点数目
            simulation.dx=simulation.Lx/num_grid;   
            % v_th=0，满足CFL条件
        end
%         % 考虑CFL条件 vth_e*Δt<Δx
%         constants=get_constants();
%         vth_e=sqrt(-constants.q_m_ration_e*simulation.Te0);
%         simulation.dt=min([simulation.dt simulation.dx/vth_e]);
        %%%%%%%%%% 其他
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        % 电磁场
        simulation.field_boundaries_type=0;%电位边界条件类型
        simulation.field_boundaries=[0,0];%电位边界条件值
        simulation.B0=0; %代表性磁感应强度 T
        % 粒子
        simulation.particle_boundaries_type=0;%粒子边界条件类型

    case 'auto1' %输入ne, Te, Lx, end_time，大致满足离散限制条件,可通过type2的check
        % 解析输入
        p = inputParser;
        validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
        validScalarNotNegNum = @(x) isnumeric(x) && isscalar(x) && ~(x < 0);
        addRequired(p,'type');
        addOptional(p,'ne0',1E16,validScalarPosNum);
        addOptional(p,'Te0',1,validScalarNotNegNum);
        addOptional(p,'end_time',6E-6,validScalarPosNum);
        addOptional(p,'Lx',0.01,validScalarPosNum);
        parse(p,type,varargin{:});
        % p.Results
        % 等离子体参数，以电子为代表
        simulation.ne0=p.Results.ne0;%等离子体密度
        simulation.Te0=p.Results.Te0;%单位eV
        simulation.end_time=p.Results.end_time;%仿真时间总长
        simulation.Lx=p.Results.Lx;%仿真区域长度
        % 时空尺度
        %%%%%%%%%% 时空步长
        % 考虑 蛙跳法求解简谐运动的精度要求 ωp*dt/2<<1
        omega_pe=get_omega_pe(simulation.ne0);
        simulation.dt=(2/omega_pe)*(pi/20); % dt=0.157*(2/ωp)=Tpe/20
        if simulation.Te0>0
            % 考虑数值自加热 %宽松限制dx<3*λD
            debye_length=get_debye_length(simulation.ne0,simulation.Te0);
            num_grid=ceil(simulation.Lx/(3*debye_length)); 
            simulation.num_grid_point=num_grid+1;%网格格点数目
            simulation.dx=simulation.Lx/num_grid;
            % 满足自加热限制条件时，CFL条件弱于蛙跳法精度要求，一般均能满足。
        else
            num_grid=20;
            simulation.num_grid_point=num_grid+1;%网格格点数目
            simulation.dx=simulation.Lx/num_grid;   
            % v_th=0，满足CFL条件
        end
%         % 考虑CFL条件 vth_e*Δt<Δx
%         constants=get_constants();
%         vth_e=sqrt(-constants.q_m_ration_e*simulation.Te0);
%         simulation.dt=min([simulation.dt simulation.dx/vth_e]);
        %%%%%%%%%% 其他
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        % 电磁场
        simulation.field_boundaries_type=0;%电位边界条件类型
        simulation.field_boundaries=[0,0];%电位边界条件值
        simulation.B0=0; %代表性磁感应强度 T
        % 粒子
        simulation.particle_boundaries_type=0;%粒子边界条件类型
    otherwise
        error('Not Done');
end

