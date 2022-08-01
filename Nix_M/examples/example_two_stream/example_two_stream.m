%1D PIC
% 双流不稳定性

%% 初始化
clear
close all;

if exist('get_path_init','file')~=2
    addpath('../../packages') % 添加./packages到路径
end
path_list=get_path_init('two_stream'); %文件路径

% 如无说明，则单位为国际单位制
% 全局变量
constants=get_constants();% 全局常数 结构体

%--------仿真参数-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %仿真参数 结构体
% 可通过取消以下代码注释来修改仿真参数。注意一些参数的内在关联。
% simulation.dt=10*10^-12;%时间步长
% simulation.end_time=6E-6;%仿真时间总长
simulation.all_timesteps=1500;%总循环次数
% simulation.Lx=0.01;%仿真区域长度
% simulation.source_region=[0.1*simulation.Lx,0.4*simulation.Lx];
% simulation.num_grid_point=201;%网格格点数目
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
% simulation.ne0=1E16;%等离子体密度
% simulation.num0_macro_e=20000;%初始电子数目 Particle e Num
% simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
% simulation.Te0=1;%e单位eV
simulation.field_boundaries_type=3;%电位边界条件类型
% simulation.field_boundaries=[0,0];%电位边界条件值

simulation.num0_macro_e=20000;%初始电子数目 Particle e Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重

veth=sqrt(constants.e*simulation.Te0/constants.me);%电子温度对应的热速度
% vb=4*veth; %定向速度
vb=0;
num_b=simulation.num0_macro_e/2;
rho_positive_background=simulation.ne0*constants.e; % 等量均布冷离子作为背景
%--------仿真参数-----------------------------------------------------------------------------------

%--------初始化-------------------------------------------------------------------------------------
% TODO：粒子结构体，与场结构体
%--------粒子数据
% 热电子束
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,1]);
ve(1:num_b)=vb+ve(1:num_b);
ve((num_b+1):end)=-vb+ve((num_b+1):end);
% % 冷电子束
% ve=ones(simulation.num0_macro_e,1);
% ve(1:num_b)=vb*ve(1:num_b);
% ve((num_b+1):end)=-vb*ve((num_b+1):end);
%初始速度分布直接规定为-dt/2时刻的值，直接开始蛙跳法
Ek=@(v,q_m_ratio) v.*v*3/(2*abs(q_m_ratio)); %动能，单位eV
Eke0=simulation.weight*sum(Ek(ve,constants.q_m_ration_e)); %初始1D总动能

% 均匀随机分布，自带扰动
position_e=get_position_init(simulation, 'entire domain uniform random', [simulation.num0_macro_e,1] );

ae=zeros(simulation.num0_macro_e,3);%加速度初值 电子加速度
%--------场数据
[ u, A, b_extra ] = get_field_init( simulation );%离散泊松方程
rho=zeros(simulation.num_grid_point,1);%净正电荷数密度初值
E=zeros(simulation.num_grid_point+1,1);%半节点电场初值  NG-1+两个墙内的半格点
B=zeros(simulation.num_grid_point+1,1);%半节点磁场初值
%
%--------初始化--------------------------------------------------------------------------------------

% 后处理相关参数
dptime=10;%每多少步长更新显示
avgSteps=499;%多少步长内的平均值
jj=499;
usum=0;
recordFlag=0;
u_record=zeros(simulation.num_grid_point,1);%暂时不用

log_name = get_log_init( path_list.output,'' );%输出初始日志
h_fig1 = figure('Unit','Normalized','position',...
    [0.02 0.3 0.6 0.6]); %实时诊断使用的大窗口
% 输出动画
% h_fig2=figure;
% video = VideoWriter('two_stream_instability.avi');
% open(video);
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
    % Pev=Pev+ae*dt;  %蛙跳法更新电子位置
    ve(:,1)=ve(:,1)+ae(:,1)*simulation.dt;
    position_e=position_e+ve(:,1)*simulation.dt;
    %------推动粒子结束-----------------------------------------------------------------------------------------------
    
    %------粒子越界处理------------------------------------------------------------------------------------------------
    % -------周期边界------------------------------------------------------
    tempFlag=(position_e<0);
    position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
    tempFlag=(position_e>simulation.Lx);
    position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
    %     参考 谢华生 pic_2strm.m 粒子周期边界条件
    %     xp=xp./L+10.0;
    %     xp=L.*(xp-floor(xp));
    % 回避了条件语句的周期边界，等效于
    %     if xp>0
    %         xp=mod(xp,L) 取余
    %     elseif -10*L<xp<0
    %         xp=L-mod(|xp|,L)
    %     end
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
    %rhoP=rhoP*0;
    rho=rhoe+rho_positive_background;
    %-----------分配电荷结束------------------------------------------------------------
    
    %--------------求解电场---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------求解电场结束---------------------------------------------------------------
    
    % TODO：待获得E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    
    % 电压累加，为  实时后处理中取平均  做准备
    if rem(ti+jj,dptime)==1
        usum=usum+u;
        jj=jj-1;
        if jj==-1
            jj=avgSteps;
        end
    end
    
    
    %-----------实时诊断----------------------------------------------------------------------------------
    if rem(ti-1,dptime)==0%每隔dptime显示一次
        display_i=floor(ti/dptime)+1;
        position_e1=position_e(1:num_b,1);
        position_e2=position_e((num_b+1):end,1);
        ve1=ve(1:num_b,1);
        ve2=ve((num_b+1):end,1);
        
        fprintf('当前时步数：%d \n',ti)
        
        e1_num(display_i)=size(position_e1,1);%区域内的电子总量变化
        e2_num(display_i)=size(position_e2,1);
        e1_vAve(display_i)=mean(abs(ve1(:,1)));%x方向上的平均速率变化
        e2_vAve(display_i)=mean(abs(ve2(:,1)));
        Eke(display_i)=simulation.weight*sum(Ek(ve,constants.q_m_ration_e))/Eke0; %归一化1D总动能
        Ep(display_i)=0.5*constants.eps0*simulation.dx*sum(E(2:end-1).*E(2:end-1))/constants.e/Eke0; %归一化1D总静电能，单位eV
        
        % 粒子电势能
        Phi_e=((simulation.dx-eassigndx).*u(enearx)+eassigndx.*u(enearx+1))/simulation.dx;
        Epe(display_i)=simulation.weight*sum(abs(Phi_e))/Eke0;
        
        figure(h_fig1)
        subplot(3,2,1);
        %         plot_no_versus_position( position_e1,position_e2,simulation );
        % 画图显示 粒子编号-粒子位置
        temp_vec=1:simulation.num0_macro_e/2;
        scatter(position_e1,temp_vec,1,'b')
        axis([0 simulation.Lx 0 simulation.num0_macro_e]);
        hold on
        scatter(position_e1,simulation.num0_macro_e/2+temp_vec,1,'r')
        %         title('粒子空间分布');
        xlabel('x [m]');
        ylabel('粒子编号');
        legend('+x向 e','-x向 e')
        hold off
        subplot(3,2,2);
        plot_num_versus_timestep( e1_num,e2_num,dptime )
        subplot(3,2,3);
        plot_1u1E_versus_x( u, ti, E, simulation )
        subplot(3,2,4);
        plot_Ek_versus_timestep( e1_vAve,e2_vAve, constants.q_m_ration_e, constants.q_m_ration_e, '+x向 e','-x向 e',dptime )
        ylabel('mean Ek [eV]'); %实际上是3倍x向平均动能
        subplot(3,2,5);
        % 能量守恒-谢华生  2020/11/05 17:04:51 失败：为什么不在同一量级？
        % 画图显示 粒子总动能与电场总储能-时间步
        plot(Eke,'-b')
        hold on
        plot(Epe,'--r')
        plot(Eke+Epe,'-.k')
        xlabel(['t [' num2str(dptime) '时步]']);
        ylabel('归一化能量');
%         axis([0,inf,0,3])
        L1=legend('Ek','Ep','Ek+Ep');
        set(L1,'location','east');
        hold off
        subplot(3,2,6);
        plot_vineV_versus_position( ve1, ve2, position_e1,position_e2,constants.q_m_ration_e, constants.q_m_ration_e, '+x向 e','-x向 e', simulation )
        drawnow;
        
%         % 输出动画
%         figure(h_fig2) %使用唯一ID
%         plot_vineV_versus_position( ve1, ve2, position_e1,position_e2,constants.q_m_ration_e, constants.q_m_ration_e, '+x向 e','-x向 e', simulation )
%         M(display_i) = getframe(h_fig2);
%         writeVideo(video,M(display_i) );
%         hold off
        
        usum=0; %重置
    end
    %-----------实时后处理----------------------------------------------------------------------------------
    
    if ti==500000 %稳定判据
        recordFlag=1;
    end
    if recordFlag>0 %稳定后取电压，大量时步取平均以降噪
        u_record(:,recordFlag)=u;
        recordFlag=recordFlag+1;
    end
    
end
%----主循环结束-----------------------------------------------------------------

now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp(now_str)
diary off

figure(h_fig1)
saveas(gcf,[path_list.output '/real_time_display.png'])

% 归一化静电能-时间步
inverse_plasma_frequency=2*pi/sqrt(simulation.ne0*constants.e^2/(constants.eps0*constants.me));
plot_t=(dptime*simulation.dt/inverse_plasma_frequency)*(1:display_i);
figure
Ep2=log(sqrt(Ep));
plot(plot_t,Ep2,'-r')
xlabel('t [2\pi\omega_{pe}^{-1}]');
ylabel('ln(E/Emax)/2');
hold off
saveas(gcf,[path_list.output '/静电能时间演化.png'])

% 输出动画
close(video);
% figure(4)
% movie(M,1)

% save('./final_data.mat')
