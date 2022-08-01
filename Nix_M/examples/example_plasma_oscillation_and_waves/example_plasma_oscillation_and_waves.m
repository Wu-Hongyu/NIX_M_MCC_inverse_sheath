%1D PIC
% 等离子体振荡与波
% 使用：修改wave_type、Te
% wave_type：1*1 char，其中wave_type(1)
% 0. 并不严格符合等离子体振荡理论模型，但在初始时接近
% 1. 磁场下电子伯恩斯坦波

%% 初始化
clear
close all;

if exist('get_path_init','file')~=2
    addpath('../../packages') % 添加./packages到路径
end
path_list=get_path_init('plasma_waves'); %文件路径
% test_all % 运行全部测试

% 如无说明，则单位为国际单位制
% 全局变量
constants=get_constants();% 全局常数 结构体

%--------仿真参数-----------------------------------------------------------------------------------
wave_type='1'; 
ne0=1E16;
Te0=0; %冷电子
% Te0=1; %热电子
end_time=6e-6;
Lx=0.01;
simulation=get_simulation('auto1',ne0, Te0, end_time, Lx);
% 等离子体波参数
switch str2double(wave_type(1))
    case 0 % 等离子体（静电）振荡
        % Δt<T_pe/10, end_time>10*T_pe
        simulation.B0=0; %无磁场
        result=check_simulation_parameters( simulation, 1 );
        end_time=10*result.T_pe;
        simulation=get_simulation('auto0',ne0, Te0, end_time, Lx);
        simulation.dt=min([simulation.dt,result.T_pe/10]);
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);
    case 1 % 磁场下电子伯恩斯坦波
        % Δt<T_ce/10, end_time>10*T_ce and end_time>T_pe, Lx>2*r_ce
        simulation.B0=1; %均匀恒定磁场
        result=check_simulation_parameters( simulation, 1 );
        end_time=max([12*result.T_ce, result.T_pe]);
        Lx=max([simulation.Lx, 2*result.r_ce]);
        simulation=get_simulation('auto0',ne0, Te0, end_time, Lx);
        simulation.dt=min([simulation.dt,result.T_ce/20]);
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);
        simulation.B0=1;
end
% 无界等离子体
simulation.field_boundaries_type=3;%电位边界条件类型
simulation.particle_boundaries_type=0;
%--------粒子参数-----------------------------------------------------------------------------------
% TODO：部分待加入particle_group
% simulation.num0_macro_e=20000;%初始电子数目 Particle e Num
simulation.num0_macro_e=2000;
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
% 均布静止正离子背景
rho_positive_background=ne0*constants.e;  % 正电荷密度
%--------实时诊断与后处理相关参数---------------------------------------------------------------------------
switch str2double(wave_type(1))
    case 0 % 等离子体（静电）振荡
        % Δt<T_pe/10, end_time>10*T_pe
        display_timesteps=max(floor([simulation.dt,result.T_pe/10]/simulation.dt));
    case 1 % 磁场下电子伯恩斯坦波
        % Δt<T_ce/10, end_time>10*T_ce and end_time>T_pe
%         display_timesteps=max(floor([simulation.dt,result.T_ce/10]/simulation.dt));
        display_timesteps=max(floor([simulation.dt,result.T_ce/20]/simulation.dt));
end
average_timesteps=499;%多少步长内的平均值
node_sample=floor(simulation.num_grid_point/2);

%输出初始日志
log_name = get_log_init( path_list.output, '' );
% diary(log_name) % 重定向方式输出日志
disp(simulation)
result=check_simulation_parameters( simulation, 2 );
disp('')
fprintf('display_timesteps=%d ;每多少步长更新显示 \n',display_timesteps)
fprintf('average_timesteps= %d ;多少步长内的平均值 \n',average_timesteps)
diary off
pause %一般只在设计参数时使用，以查看控制区信息

%--------初始化-------------------------------------------------------------------------------------
% TODO：粒子结构体，与场结构体
%--------粒子数据
%按照Maxwell分布生成宏粒子的速度分布 x y z方向(Particle e velocity)
ve=get_v_init( result.vth_e, 'Maxwellian velocity', [simulation.num0_macro_e,3]);




%%初始速度分布直接规定为-dt/2时刻的值，直接开始蛙跳法
switch str2double(wave_type(1))
    case 0 % 等离子体（静电）振荡 1D1V
        ve(:,2:3)=0;
        Ek=@(v,q_m_ratio) 3*sum(v.*v/(2*abs(q_m_ratio))); %动能，单位eV
    case 1 % 磁场下电子伯恩斯坦波 1D3V
        ve=get_v_init( result.vth_e, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
        Ek=@(v,q_m_ratio) sum(v.*v/(2*abs(q_m_ratio))); %动能，单位eV
end
Eke0=simulation.weight*sum(Ek(ve,constants.q_m_ratio_e)); %初始1D总动能
if 0==Eke0
    Eke0=1;
end

% 均匀随机分布，自带扰动
position_e=get_position_init(simulation, 'entire domain uniform random', [simulation.num0_macro_e,1] );

ae=zeros(simulation.num0_macro_e,3);%加速度初值 电子加速度
%--------场数据
[ u, A, b_extra ] = get_field_init( simulation );%离散泊松方程
rho=zeros(simulation.num_grid_point,1);%净正电荷数密度初值
E=zeros(simulation.num_grid_point+1,1);%半节点电场初值  NG-1+两个墙内的半格点
B=zeros(simulation.num_grid_point+1,1);%半节点磁场初值
%--------诊断数据
j_average_counter=average_timesteps;
usum=0;
record_flag=0; %初始为0
u_record=zeros(simulation.num_grid_point,1);
%--------初始化--------------------------------------------------------------------------------------

h_fig1 = figure('Unit','Normalized','position',...
    [0.02 0.3 0.6 0.6]); %实时诊断使用的大窗口
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------主循环开始-----------------------------------------------------------------------------------

for ti=1:simulation.all_timesteps%主循环
    % for ti=1:30000%主循环
    
    % 主循环中性能考虑
    % 1. 简单代码可以inline。但同样建立函数文件，以便于test；可将test后代码复制修改。
    % 2. 使用专用函数（通过函数名而非传入参数唯一标识），注意传参耗时，不使用全局变量，
    % 尽量避免传值复制，参考 https://www.zhihu.com/question/50408548/answer/120840847
    % 3. 在package的test文件中，测试功能与性能
    
    %------推动粒子-----------------------------------------------------------------------------------------------
    %蛙跳法+磁场矢量分析更新电子位置
    ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ratio_e*simulation.B0*ve(:,3))*simulation.dt;
    ve(:,3)=ve(:,3)+constants.q_m_ratio_e*simulation.B0*simulation.dt*ve(:,1);% By
    position_e=position_e+ve(:,1)*simulation.dt;
    %------推动粒子结束-----------------------------------------------------------------------------------------------
    
    %------粒子越界处理------------------------------------------------------------------------------------------------
    switch simulation.particle_boundaries_type
        case 0 % 周期边界
            % -------周期边界------------------------------------------------------
            tempFlag=(position_e<0);
            position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
            tempFlag=(position_e>simulation.Lx);
            position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
        otherwise
            throw(get_exception('error','no such type of particle boundaries.'))
    end
    %------粒子越界处理结束------------------------------------------------------------------------------------------
    
    %-----------分配电荷----------------------------------------------------------------------------------
    rhoea=zeros(simulation.num_grid_point,1);%网格左右两个点的粒子密度 电子
    rhoeb=zeros(simulation.num_grid_point,1);
    
    enearx=floor(position_e/simulation.dx)+1;%每个电子所在网格的索引 整数格点
    eassigndx=position_e-(enearx-1)*simulation.dx;%每个电子距离左边最近的网格的距离
    Eenearx=floor((position_e+0.5*simulation.dx)/simulation.dx)+1;%每个电子所在半网格的索引 半格点
    Eedx=position_e-(Eenearx*simulation.dx-3/2*simulation.dx);%每个电子距离左边最近的半网格的距离
    
    for j=1:size(position_e)
        rhoea(enearx(j))=rhoea(enearx(j))+(simulation.dx-eassigndx(j));
        rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)+(eassigndx(j));
        %         rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)-e/(dx^2)*(assigndx(j))*weight;
    end
    rhoe=(rhoea+rhoeb)*(-constants.e)/(simulation.dx^2)*simulation.weight;
    rho=rhoe+rho_positive_background;
    %-----------分配电荷结束------------------------------------------------------------
    
    %--------------求解电场---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------求解电场结束---------------------------------------------------------------
    
    % TODO：待获得E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    
    % 电压累加，为  实时诊断中取平均  做准备
    if rem(ti+j_average_counter,display_timesteps)==1
        usum=usum+u;
        j_average_counter=j_average_counter-1;
        if j_average_counter==-1
            j_average_counter=average_timesteps;
        end
    end
    
    %-----------实时诊断----------------------------------------------------------------------------------
    if rem(ti-1,display_timesteps)==0%每隔dptime显示一次
        fprintf('当前时步数：%d/%d \n',ti,simulation.all_timesteps)
        
        display_i=ceil(ti/display_timesteps);
        e_num(display_i)=size(position_e,1);%区域内的电子总量变化
        e_vAve(display_i)=mean(abs(ve(:,1)));%x方向上的平均速率变化
        Eke(display_i)=simulation.weight*sum(Ek(ve,constants.q_m_ratio_e))/Eke0; %归一化1D总动能
        Ep(display_i)=0.5*constants.eps0*simulation.dx*sum(E(2:end-1).*E(2:end-1))/constants.e/Eke0; %归一化1D总静电能，单位eV
        % 粒子电势能
        Phi_e=((simulation.dx-eassigndx).*u(enearx)+eassigndx.*u(enearx+1))/simulation.dx;
        Epe(display_i)=simulation.weight*sum(abs(Phi_e))/Eke0;
        
        % 粒子平均能量
        ave_Eke(display_i)=sum(Ek(ve,constants.q_m_ratio_e))/simulation.num0_macro_e; 
        ave_Epe(display_i)=sum(-Phi_e)/simulation.num0_macro_e;
        
        E_node_sample(display_i)=E(node_sample);
        u_node_sample(display_i)=u(node_sample);
        
        figure(h_fig1)
        subplot(3,2,1);
%         plot_no_ versus_position( position_e,0,simulation );
        hold on
        position_H=(0:simulation.num_grid_point-1)*simulation.dx;
        for i=1:simulation.num_grid_point
             line([position_H(i),position_H(i)],[0,simulation.num0_macro_e],'linestyle','-.','color','r','LineWidth',1);
        end
        hold off
        subplot(3,2,2);
%          plot_num_versus_timestep( e_num,0*e_num,display_timesteps )
        subplot(3,2,3);
        plot_1u1E_versus_x( u, ti, E, simulation )
        %         subplot(3,2,4);
        subplot(3,2,5);%电荷分布
        % 待根据test_energy.m修改
%         % 画图显示 粒子总动能与电场总储能-时间步
%         plot(ave_Eke,'-b')
%         hold on
%         plot(ave_Epe,'--r')
%         plot(ave_Eke+ave_Epe,'-.k')
%         xlabel(['t [' num2str(display_timesteps) '时步]']);
%         ylabel('归一化能量');
%         %         axis([0,inf,0,1.2])
%         L1=legend('Ek','Ep','Ek+Ep');
%         set(L1,'location','east');
%         hold off
        subplot(3,2,6);
        %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%         plot_vineV_versus_position( ve, 0*ve, position_e,0*position_e,constants.q_m_ratio_e, constants.q_m_ratio_e, 'e','', simulation )
%         ylim([-200,200])
        drawnow;
        
        usum=0; %重置
    end
    %-----------实时诊断----------------------------------------------------------------------------------
    % 为主循环终止后的后处理做准备
    if ti==500000 %记录电压的时间步
        record_flag=1;
    end
    if record_flag>0 %从500000一直取到终止，大量时步取平均以降噪
        u_record(:,record_flag)=u;
        record_flag=record_flag+1;
    end
    
    %终止判据：达到总时步数时终止
end
%----主循环结束-----------------------------------------------------------------

now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp([now_str ' 主循环结束'])
diary off
now_str=datestr(now,'HH_MM');

% 终止时实时诊断图像（过程数据的示意）存储
figure(h_fig1)
saveas(gcf,[path_list.output '/终止时实时诊断' now_str '.png'])

%% 频谱分析
% 采样信号信息
dt_sampling=display_timesteps*simulation.dt; %采样时间间隔
f_sampling=1/dt_sampling; %采样频率
L_sampling=length(u_node_sample); %采样点数
% 以最后一个采样物理量实际长度为准，便于中断后做分析
if mod(L_sampling,2) % 需为偶数
    L_sampling=L_sampling-1;
    E_node_sample=E_node_sample(1:L_sampling);
    u_node_sample=u_node_sample(1:L_sampling);
end

display_t=dt_sampling*(1:L_sampling); %采样时间点

figure
subplot(2,1,1) % 采样信号
title('采样信号')
switch str2double(wave_type(1))
    case 0 % 等离子体（静电）振荡
        % Δt<T_pe/10, end_time>10*T_pe
        plot_t=display_t/result.T_pe;
        plot(plot_t,E_node_sample,'-r')
        xlabel('t [2\pi\omega_{pe}^{-1}]');
        ylabel('E_{node sample}');
        for i=1:5
            line([i,i],[min(E_node_sample),max(E_node_sample)],'linestyle',':','color','k');
        end
        hold off
    case 1 % 磁场下电子伯恩斯坦波
        % Δt<T_ce/10, end_time>10*T_ce and end_time>T_pe, Lx>2*r_ce
        plot_t=display_t/result.T_ce;
        plot(plot_t,E_node_sample,'-r')
        xlabel('t [2\pi\omega_{ce}^{-1}]');
        ylabel('E_{node sample}');
        for i=1:5
            line([i,i],[min(E_node_sample),max(E_node_sample)],'linestyle',':','color','k');
        end
        hold off
end

% fft
Y_fft=fft(E_node_sample); %幅度频谱
P2 = abs(Y_fft/L_sampling); %双边幅度频谱
f_plot2=f_sampling*(1:L_sampling)/L_sampling;

P1 = P2(1:L_sampling/2+1);
P1(2:end-1) = 2*P1(2:end-1); %单边幅度频谱
f_plot1 = f_sampling*(0:(L_sampling/2))/L_sampling;

subplot(2,1,2) %幅值频谱
switch str2double(wave_type(1))
    case 0 % 等离子体（静电）振荡
        % Δt<T_pe/10, end_time>10*T_pe
        plot_omega=2*pi*f_plot1/result.omega_pe;
        h_fig4=plot(plot_omega,P1,'b-'); %单边幅度频谱
        xlabel('\omega/\omega_{pe}')
        ylabel('|P1(f)|')
        for i=1:floor(max(plot_omega))
            line([i,i],[0,max(P1)],'linestyle',':','color','k');
        end
        
        
    case 1 % 磁场下电子伯恩斯坦波
        % Δt<T_ce/10, end_time>10*T_ce and end_time>T_pe, Lx>2*r_ce
        yyaxis left;
        % plot(f_plot2,P2) %双边幅度频谱
        % plot(f_plot1,P1) %单边幅度频谱
        % xlabel('f [Hz]')
        % ylabel('|P1(f)|')
        plot_omega=2*pi*f_plot1/result.omega_ce;
        h_fig4=plot(plot_omega,P1,'b-'); %单边幅度频谱
        xlabel('\omega/\omega_{ce}')
        ylabel('|P1(f)|')
        for i=1:floor(max(plot_omega))
            line([i,i],[0,max(P1)],'linestyle',':','color','k');
        end
        yyaxis right;
        h_fig5=plot(plot_omega,log10(P1),'--','Color',[0.8500    0.3250    0.0980]); %单边幅度频谱取对数
        ylabel('lg|P1(f)|')
        legend([h_fig4,h_fig5],'幅值','lg(幅值)')
        % % 人工限制范围，以突出信息
        % yyaxis left;
        % ylim([0,7]) 
end
saveas(gcf,[path_list.output '/信号与频谱' now_str '.png'])

%终止时全部数据存储
save([path_list.output '/final_data' now_str '.mat'])