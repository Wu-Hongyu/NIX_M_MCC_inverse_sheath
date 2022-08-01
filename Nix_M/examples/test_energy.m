%1D PIC
% 使用：根据物理模型，依次设置各类参数，进行相应初始化，
% 并在主循环中使用相应模块（尤其是pusher和粒子边界条件）

%% 初始化
clear
close all;

if exist('get_path_init','file')~=2
    addpath('./packages') % 添加./packages到路径
end
path_list=get_path_init('now'); %文件路径
% test_all % 运行全部测试

% 如无说明，则单位为国际单位制
% 全局变量
constants=get_constants();% 全局常数 结构体

%--------仿真参数-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %仿真参数 结构体
% 可通过取消以下代码注释来修改仿真参数
% % 等离子体参数，以电子为代表
% simulation.ne0=1E16;%等离子体密度
simulation.Te0=1;%单位eV
% %         TH=1;%单位eV
% % 时空尺度
% simulation.dt=10*10^-12;%时间步长
% simulation.end_time=6E-6;%仿真时间总长
% simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
% simulation.Lx=0.01;%仿真区域长度
% simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
% simulation.num_grid_point=201;%网格格点数目
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
% % 电磁场
simulation.field_boundaries_type=3;%电位边界条件类型
% simulation.field_boundaries=[0,0];%电位边界条件值
% simulation.B0=0; %代表性磁感应强度 T
% 粒子
simulation.particle_boundaries_type=0;%粒子边界条件类型

%--------粒子参数-----------------------------------------------------------------------------------
% TODO：部分待加入particle_group
simulation.TH0=0;%单位eV
simulation.num0_macro_e=20000;%初始电子数目 Particle e Num
simulation.num0_macro_H=simulation.num0_macro_e;%初始离子数目Particle H Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
% TODO: 现在可以通过result=check_simulation_parameters返回值获得result.vth_e
veth=sqrt(-constants.q_m_ration_e*simulation.Te0);%电子温度对应的热速度
vHth=sqrt(constants.q_m_ration_H*simulation.TH0);%离子温度对应的热速度

%--------诊断参数---------------------------------------------------------------------------
display_timesteps=100;%每多少步长更新显示
average_timesteps=499;%多少步长内的平均值

%输出初始日志
log_name = get_log_init( path_list.output, '' );
diary(log_name) % 重定向方式输出日志
disp(simulation)
result_default=check_simulation_parameters( simulation, 1 );
disp('')
fprintf('display_timesteps=%d ;%每多少步长更新显示\r\n',display_timesteps)
fprintf('average_timesteps= %d ;%多少步长内的平均值\r\n',average_timesteps)
diary off

%--------初始化-------------------------------------------------------------------------------------
% TODO：粒子结构体，与场结构体
%--------粒子数据
%按照Maxwell分布生成宏粒子的速度分布 x y z方向(Particle e velocity)
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
%%初始速度分布直接规定为-dt/2时刻的值，直接开始蛙跳法

%均匀分布在空间Lx内，并给一个随机扰动
position_e=get_position_init(simulation, 'entire domain uniform', [simulation.num0_macro_e,1] );
position_H=get_position_init(simulation, 'entire domain uniform', [simulation.num0_macro_H,1] );



ae=zeros(simulation.num0_macro_e,3);%加速度初值 电子加速度
aH=zeros(simulation.num0_macro_H,3);%加速度初值 离子加速度
%--------场数据
[ u, A, b_extra ] = get_field_init( simulation );%离散泊松方程
rho=zeros(simulation.num_grid_point,1);%净正电荷数密度初值
E=zeros(simulation.num_grid_point+1,1);%半节点电场初值  NG-1+两个墙内的半格点
B=zeros(simulation.num_grid_point+1,1);%半节点磁场初值
%--------诊断数据
j_average_counter=average_timesteps;
usum=0;
record_flag=0;
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
    % Pev=Pev+ae*dt;  %蛙跳法更新电子位置
    ve(:,1)=ve(:,1)+ae(:,1)*simulation.dt;
    position_e=position_e+ve(:,1)*simulation.dt;
    
    %     %蛙跳法+磁场矢量分析更新电子位置
    %     ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ration_e*simulation.B0*ve(:,3))*simulation.dt;
    %     ve(:,3)=ve(:,3)+constants.q_m_ration_e*simulation.B0*simulation.dt*ve(:,1);% By
    %     position_e=position_e+ve(:,1)*simulation.dt;
    
    vH(:,1)=vH(:,1)+aH(:,1)*simulation.dt;%蛙跳法更新离子位置
    position_H=position_H+vH(:,1)*simulation.dt;
    
    
    %------推动粒子结束-----------------------------------------------------------------------------------------------
    
    %------粒子越界处理------------------------------------------------------------------------------------------------
    switch simulation.particle_boundaries_type
        case 0 % 周期边界
            % -------周期边界------------------------------------------------------
            tempFlag=(position_e<0);
            position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
            tempFlag=(position_e>simulation.Lx);
            position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
            tempFlag=(position_H<0);
            position_H(tempFlag)=position_H(tempFlag)+simulation.Lx;
            tempFlag=(position_H>simulation.Lx);
            position_H(tempFlag)=position_H(tempFlag)-simulation.Lx;
        case 1 % 成对重注入+热化
            % -------成对重注入源区 pair re-injection into source region--------------
            tempFlag=(position_e<0|position_e>simulation.Lx);
            ve(tempFlag,:)=[];%抹除超出区域的电子
            position_e(tempFlag)=[];
            
            k1=(position_H<0|position_H>simulation.Lx);%取超出范围的离子的编号
            tempsum=sum(k1);
            
            if(tempsum>0)%如果非空
                vH(k1,:)=normrnd(0,vHth,[tempsum,3]);%xyz 重注入离子速度Maxwell分布
                % 重注入粒子位置均匀随机分布在源区内
                position_H(k1)=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
                temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
                position_e=reshape([position_e; temp_matrix],[],1);%合并回原来的数组
                ve=reshape([ve; normrnd(0,veth,[tempsum,3])],[],3);
            end
            
            %--------thermalization 热化  每50步一次
            % 搭配非粒子周期边界使用
            if rem(ti,50)==1
                flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
            end
        otherwise
            throw(get_exception('error','no such type of particle boundaries.'))
    end
    %------粒子越界处理结束------------------------------------------------------------------------------------------
    
    %-----------分配电荷----------------------------------------------------------------------------------
    rhoea=zeros(simulation.num_grid_point,1);%网格左右两个点的粒子密度 电子
    rhoeb=zeros(simulation.num_grid_point,1);
    rhoaH=zeros(simulation.num_grid_point,1);%离子
    rhobH=zeros(simulation.num_grid_point,1);
    
    enearx=floor(position_e/simulation.dx)+1;%每个电子所在网格的索引 整数格点
    eassigndx=position_e-(enearx-1)*simulation.dx;%每个电子距离左边最近的网格的距离
    Eenearx=floor((position_e+0.5*simulation.dx)/simulation.dx)+1;%每个电子所在半网格的索引 半格点
    Eedx=position_e-(Eenearx*simulation.dx-3/2*simulation.dx);%每个电子距离左边最近的半网格的距离
    
    Hnearx=floor(position_H/simulation.dx)+1;
    Hassigndx=position_H-(Hnearx-1)*simulation.dx;
    EHnearx=floor((position_H+0.5*simulation.dx)/simulation.dx)+1;%每个离子所在半网格的索引 半格点
    EHdx=position_H-(EHnearx*simulation.dx-3/2*simulation.dx);
    
    for j=1:size(position_e)
        rhoea(enearx(j))=rhoea(enearx(j))+(simulation.dx-eassigndx(j));
        rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)+(eassigndx(j));
        %         rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)-e/(dx^2)*(assigndx(j))*weight;
    end
    rhoe=(rhoea+rhoeb)*(-constants.e)/(simulation.dx^2)*simulation.weight;
    %rhoP=rhoP*0;
    for j=1:size(position_H)
        rhoaH(Hnearx(j))=rhoaH(Hnearx(j))+(simulation.dx-Hassigndx(j));
        rhobH(Hnearx(j)+1)=rhobH(Hnearx(j)+1)+(Hassigndx(j));
        %         rhobH(Hnearx(j)+1)=rhobH(Hnearx(j)+1)+e/(dx^2)*(Hassigndx(j))*weight;
    end
    rhoH=(rhoaH+rhobH)*constants.e/(simulation.dx^2)*simulation.weight;
    
    rho=rhoe+rhoH;
    %-----------分配电荷结束------------------------------------------------------------
    
    %--------------求解电场---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------求解电场结束---------------------------------------------------------------
    
    % TODO：待获得E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    aH=((simulation.dx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/simulation.dx;
    
    % 电压累加，为  实时诊断中取平均  做准备
    if rem(ti+j_average_counter,display_timesteps)==1 %ti处于display的前average_timesteps范围中
        usum=usum+u;
        j_average_counter=j_average_counter-1;
        if j_average_counter==-1 %ti已超出display的前average_timesteps范围
            j_average_counter=average_timesteps; %重置j_average_counter
        end
    end
    
    %-----------实时诊断----------------------------------------------------------------------------------
    if rem(ti-1,display_timesteps)==0%每隔dptime显示一次
        fprintf('当前时步数：%d/%d \n',ti,simulation.all_timesteps)
        
        display_i=ceil(ti/display_timesteps);
        e_num(display_i)=size(position_e,1);%区域内的电子总量变化
        H_num(display_i)=size(position_H,1);
        e_vAve(display_i)=mean(abs(ve(:,1)));%x方向上的平均速率变化
        H_vAve(display_i)=mean(abs(vH(:,1)));
        
        
        % 粒子平均能量
        Ek=@(v,q_m_ratio) sum(v.*v/(2*abs(q_m_ratio))); %动能，单位eV
        Phi_e=((simulation.dx-eassigndx).*u(enearx)+eassigndx.*u(enearx+1))/simulation.dx;
        Phi_H=((simulation.dx-Hassigndx).*u(Hnearx)+Hassigndx.*u(Hnearx+1))/simulation.dx;
        ave_Eke(display_i)=sum(Ek(ve(:,1),constants.q_m_ration_e))/simulation.num0_macro_e; 
        ave_EkH(display_i)=sum(Ek(vH(:,1),constants.q_m_ration_H))/simulation.num0_macro_e; 
        ave_Ek=ave_EkH+ave_Eke;
        ave_Epe(display_i)=sum(-Phi_e)/simulation.num0_macro_e;
        ave_EpH(display_i)=sum(Phi_H)/simulation.num0_macro_e;
        ave_Ep=ave_Epe+ave_EpH;
        
        figure(h_fig1)
        subplot(3,2,1);
        plot_no_versus_position( position_e,position_H,simulation );
        subplot(3,2,2);
        plot_num_versus_timestep( e_num,H_num,display_timesteps )
        subplot(3,2,3);
        plot_2u1E_versus_x( u, ti, usum, average_timesteps, E, simulation )
        subplot(3,2,4);
%         plot_Ek_versus_timestep( e_vAve,H_vAve, constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', display_timesteps )
%         ylabel('平均Ek [eV]'); %实际上是3倍x向平均动能

% 能量守恒-谢华生  2020/11/05 17:04:51 失败：为什么不在同一量级？
        % 画图显示 粒子总动能与电场总储能-时间步
        plot(ave_Ek-0.5,'b')
        hold on
        plot(ave_Ep,'r')
        plot(ave_Ek+ave_Ep-0.5,'k')
        xlabel(['t [' num2str(display_timesteps) '时步]']);
        ylabel('平均能量');
        %         axis([0,inf,0,1.2])
        L1=legend('Ek','Ep','Ek+Ep');
        set(L1,'location','east');
        hold off



        subplot(3,2,5);%电荷分布
        plot_density_versus_x( rhoe, rhoH, rho, simulation )
        subplot(3,2,6);
        %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
        plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
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

% 终止（人工判断稳定）时电势空间分布
figure
x_final=(0:1:200)*simulation.dx;
u_final=mean(u_record,2);
figure
plot(x_final,u_final,'-b','LineWidth',3)%多个步长电压平均值
axis([0,simulation.Lx,-inf,inf])
%         title('电势空间分布', 'FontSize', 18);
xlabel('x [m]')
ylabel('\phi [V]');
hold on;
line([simulation.Lx-6.66e-4,simulation.Lx],[2.51,2.51],'linestyle','-.','color','r','LineWidth',3);
L1=legend('仿真','解析');
set(L1,'location','south');
set(L1,'AutoUpdate','off');
line([simulation.Lx-6.66e-4,simulation.Lx-6.66e-4],[0,2.51],'linestyle','-.','color','r','LineWidth',3);
saveas(gcf,[path_list.output '/终止时电势分布' now_str '.png'])

%终止时全部数据存储
save([path_list.output '/final_data' now_str '.mat'])

%% 后处理
% 读取存储过程数据