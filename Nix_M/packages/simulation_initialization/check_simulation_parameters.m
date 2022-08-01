function [result] = check_simulation_parameters( simulation, type )
% 检查仿真参数有效性
% 使用: % 运行本函数，根据提示修改simulation参数，再次运行本函数
% type: 检查严格程度。1全部作为警告（不中止程序），3全部报错（中止程序）
% 实现：每种检查抛出指定类型异常，分别try/catch以避免不报错时也只运行一种check
constants=get_constants();
%% 输出results参数到屏幕
result.flag_check=true; %check失败则置零
result.debye_length=get_debye_length(simulation.ne0,simulation.Te0);
result.omega_pe=get_omega_pe(simulation.ne0);
result.T_pe=2*pi/result.omega_pe;
result.N_D=get_N_D(simulation.ne0,simulation.Te0); %德拜球内粒子数目
result.vth_e=sqrt(-constants.q_m_ratio_e*simulation.Te0);
result.flag_magnetized=false;
if simulation.B0>0
    result.omega_ce=-constants.q_m_ratio_e*simulation.B0;
    result.T_ce=2*pi/result.omega_ce;
    result.r_ce=result.vth_e/result.omega_ce;
    if result.r_ce<simulation.Lx
        result.flag_magnetized=true;
    end
end
% result.ppc=

fprintf('[INFO] 时空特征尺度\n')
fprintf('电子等离子体频率 ω_pe = %.2e rad/s\n',result.omega_pe)
fprintf('模拟域长度 Lx = %.2e m;  ',simulation.Lx)
fprintf('德拜长度 λ_De = %.2e m;  ',result.debye_length)
fprintf('空间步长 Δx = %.2e m;  ',simulation.dx)
fprintf('网格节点数 N_node = %d \n',simulation.num_grid_point)
if 3==simulation.field_boundaries_type
    fprintf('无界等离子体：模拟域使用周期边界\n')
else
    fprintf('有界等离子体：模拟域使用非周期边界\n')
end
fprintf('物理时长 t_end = %.2e s;  ',simulation.end_time)
fprintf('电子静电振荡周期 T_pe = %.2e s;  ',result.T_pe)
fprintf('时间步长 Δt = %.2e s;  ',simulation.dt)
fprintf('预计总时步数 N_step = %d \n',simulation.all_timesteps)
fprintf('电子热速度 vth_e = %.2e \n',result.vth_e)
fprintf('德拜球内粒子数目 N_D = %.2e \n',result.N_D)
fprintf('代表性磁感应强度 B0 = %d T\n',simulation.B0)
if simulation.B0>0
    fprintf('拉莫尔频率 ω_ce = %d rad/s;  ',result.omega_ce)
    fprintf('回旋周期 T_ce = %d s;  ',result.T_ce)
    fprintf('拉莫尔半径 r_ce = %d m \n;  ',result.r_ce)
end
% fprintf('宏粒子代表粒子数 weighting = %.2e \n',result.N_D)
% fprintf('初始电子宏粒子数 N_macroe = %d \n',result.N_D)
% fprintf('每单元宏粒子数 PPC = %d \n',result.N_D)

%% 检查物理参数
fprintf('[INFO] 检查物理参数\n')
ratio_Lx_debye_ratio=simulation.Lx/result.debye_length;
fprintf('Lx = %.2e λ_De\n',ratio_Lx_debye_ratio)
try
    if 3==simulation.field_boundaries_type
        if ratio_Lx_debye_ratio<1
            throw(get_exception('warn','Lx < λ_De，无法识别德拜球'))
        end
    else %有界等离子体
        if ratio_Lx_debye_ratio<1
            throw(get_exception('error','Lx < λ_De，等离子体判据要求Lx>>λ_De'))
        elseif ratio_Lx_debye_ratio<10
            %     fprintf('[WARN]: Lx < 10*λ_De\n')
            throw(get_exception('warn','λ_De <= Lx < 10*λ_De，等离子体判据要求Lx>>λ_De'))
        end
    end
catch e
    treate_exception(e,type)
end

try
    if 0==simulation.Te0
        throw(get_exception('warn','Te0=0, 则N_D = 0'))
    elseif result.N_D<1
        throw(get_exception('error','N_D < 1，等离子体判据要求N_D>>>1'))
    elseif result.N_D<100
        throw(get_exception('warn','1 <= N_D < 100，等离子体判据要求N_D>>>1'))
    end
catch e
    treate_exception(e,type)
end
% TODO：计算电子与中性粒子的平均碰撞周期（等离子体判据）。待MCC模块确定后再实现

ratio_end_time_T_pe=simulation.end_time/result.T_pe;
fprintf('t_end = %.2e T_pe\n',ratio_end_time_T_pe)
try
    if ratio_end_time_T_pe<1
        throw(get_exception('warn','t_end < T_pe，无法还原等离子体振荡信息'))
        % 后面时间步长会error，因此此处只warn
    elseif ratio_end_time_T_pe<10
        throw(get_exception('warn','t_end < 10*T_pe，难以还原等离子体振荡信息'))
    end
catch e
    treate_exception(e,type)
end

% 磁化
if result.flag_magnetized
    disp('电子被磁化')
    ratio_end_time_T_ce=simulation.end_time/result.T_ce;
    fprintf('t_end = %.2e T_ce\n',ratio_end_time_T_ce)
    try
        if ratio_end_time_T_ce<1
            throw(get_exception('warn','t_end < T_ce，无法还原拉莫尔运动信息'))
            % 后面时间步长会error，因此此处只warn
        elseif ratio_end_time_T_ce<10
            throw(get_exception('warn','t_end < 10*T_ce，难以还原拉莫尔运动信息'))
        end
    catch e
        treate_exception(e,type)
    end
    
    ratio_Lx_r_ce=simulation.Lx/result.r_ce;
    fprintf('Lx = %.2e r_ce\n',ratio_Lx_r_ce)
    if ratio_Lx_r_ce>10
        fprintf('Lx > 10*r_ce, \n')
    end
end

%% 检查数值参数
fprintf('[INFO] 检查数值参数-离散步长限制\n')
% 考虑数值自加热的空间步长限制条件 dx<≈λD
ratio_dx_debye_length=simulation.dx/result.debye_length;
fprintf('Δx = %.2e λ_De\n',ratio_dx_debye_length)
try
    if 0==simulation.Te0
        throw(get_exception('warn','Te0=0, 则λ_De = 0'))
    elseif ratio_dx_debye_length>3
        throw(get_exception('error','Δx > 3*λ_De，请考虑数值加热'))
    elseif ratio_dx_debye_length>1
        throw(get_exception('warn','Δx > λ_De'))
    end
catch e
    treate_exception(e,type)
end
% 蛙跳法求解电子静电振荡的收敛条件ω_pe*Δt/2<1，精度要求 ωp*dt/2<<1
omega_pe_dt_2=result.omega_pe*simulation.dt/2;
fprintf('ω_pe*Δt/2 = %.2e rad\n',omega_pe_dt_2)
try
    if omega_pe_dt_2>1
        throw(get_exception('error','ω_pe*Δt/2>1 rad，蛙跳法不收敛，存在幅值误差'))
    elseif omega_pe_dt_2>0.1
        throw(get_exception('warn','ω_pe*Δt/2>0.1 rad，请考虑蛙跳法相位误差'))
    end
catch e
    treate_exception(e,type)
end
% 显式时间FDM求解偏微分方程的收敛条件(CFL条件)vth_e*Δt<Δx
ratio_v_dt_dx=result.vth_e*simulation.dt/simulation.dx;
fprintf('vth_e*Δt = %.2e Δx\n',ratio_v_dt_dx)
try
    if ratio_v_dt_dx>1
        throw(get_exception('error','vth_e*Δt>Δx，不满足CFL条件'))
    end
catch e
    treate_exception(e,type)
end
% TODO：result.ppc。待paticle_group确定后再实现
% 密度统计的散粒噪声正比于1/sqrt(ppc)，因此精度要求PPC>10
% TODO：碰撞对时间步长的要求。待MCC模块确定后再实现

% TODO: 断言不应该出现的错误
% assert particle_group[1]参数与simulation中电子参数相等  %也可以放到test_particle_group中

if result.flag_check
    fprintf('OK! \n')
end
end

function treate_exception(exception,type)
switch type
    case 1 %warn与error均warning
        warning(exception.message)
    case 2 % 对于warn包装进warning，对于error则rethrow
        temp_e=get_exception('warn','');
        if exception.identifier==temp_e.identifier
            warning(exception.message)
        else
            result.flag_check=false;
            rethrow(exception)
        end
    case 3 %warn与error都rethrow
        result.flag_check=false;
        rethrow(exception)
end
end