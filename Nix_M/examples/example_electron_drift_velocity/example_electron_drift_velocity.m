%在约化电场下的电子漂移速率测试
%从一端释放，然后经过一段时间后平均所有电子的速度得到漂移速度

%需要修改：
%1.电场用约化电场得到
% 2.气体密度用约化电场得到
% 3.仿真区域应修改较大  文献中250Td的约化电场中，电子的漂移速度为40E4m/s  （发现不用改，实验的漂移长度为10cm）
% 4.步长需要修改  让一个时间步长内的反应概率远小于1

% 先算一个约化电场下的例子


if exist('get_path_init','file')~=2
    addpath('../../packages') % 添加./packages到路径
end
path_list=get_path_init('drift_velocity'); %文件路径

Td_list=[50 75 100 125 150 175 200 250];
% Td_list=[50 75 100 125];
constants=get_constants();% 全局常数 结构体
simulation=get_simulation('default'); %仿真参数 结构体
simulation.num0_macro_e=20000;%初始电子数目 Particle e Num
EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
% Te=[ 10:10:50];
% Te=logspace(log10(0.1),log10(50),50);
% Te=50;

% %         TH=1;%单位eV
% % 时空尺度
simulation.dt=1*10^-12;%时间步长
simulation.end_time=2E-8;%仿真时间总长
simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
% simulation.Lx=0.01;%仿真区域长度
% simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
% simulation.num_grid_point=201;%网格格点数目
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
% % 电磁场
% simulation.field_boundaries_type=0;%电位边界条件类型
% simulation.field_boundaries=[0,0];%电位边界条件值
% simulation.B0=0; %代表性磁感应强度 T
% 粒子
simulation.particle_boundaries_type=0;%粒子边界条件类型
simulation.field_boundaries_type=3;%电位边界条件类型

ddx=simulation.dx;
ddt=simulation.dt;
%--------粒子参数-----------------------------------------------------------------------------------
% TODO：部分待加入particle_group
simulation.TH0=1;%单位eV
simulation.num0_macro_H=simulation.num0_macro_e;%初始离子数目Particle H Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
% TODO: 现在可以通过result=check_simulation_parameters返回值获得result.vth_e
veth=sqrt(-constants.q_m_ratio_e*simulation.Te0);%电子温度对应的热速度
vHth=sqrt(constants.q_m_ratio_Hp*simulation.TH0);%离子温度对应的热速度

for Tdi=1:length(Td_list)
    pressure=400;%Pa
    Td=Td_list(Tdi);
    
    %--------MCC参数-----------------------------------------------------------------------------------
%     [mcc_para,Xsec]=MCC_init(simulation,constants,pressure);
   [mcc_para,Xsec]= MCC_init(simulation,constants,simulation.pressure, '2020LZS');
    E=-Td/(1E21)*mcc_para.n_target;
    
    %--------诊断参数---------------------------------------------------------------------------
    display_timesteps=1000;%每多少步长更新显示
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
    ve=0*get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
    vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
    %%初始速度分布直接规定为-dt/2时刻的值，直接开始蛙跳法
    
    %均匀分布在空间Lx内，并给一个随机扰动
    position_e=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_e,1] );
    position_H=0*get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_H,1] );
    position_e=0.01*position_e;%分布在前1%的位置
    
    ae=zeros(simulation.num0_macro_e,3);%加速度初值 电子加速度
    % aH=zeros(simulation.num0_macro_H,3);%加速度初值 离子加速度
    %--------场数据
    [ u, A, b_extra ] = get_field_init( simulation );%离散泊松方程
    rho=zeros(simulation.num_grid_point,1);%净正电荷数密度初值
    % E=zeros(simulation.num_grid_point+1,1);%半节点电场初值  NG-1+两个墙内的半格点
    B=zeros(simulation.num_grid_point+1,1);%半节点磁场初值
    %--------诊断数据
    j_average_counter=average_timesteps;
    usum=0;
    record_flag=0;
    % u_record=zeros(simulation.num_grid_point,1);
    %--------初始化--------------------------------------------------------------------------------------
    
%     h_fig1 = figure('Unit','Normalized','position',...
%         [0.02 0.3 0.6 0.6]); %实时诊断使用的大窗口
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    collision_energy_count=zeros(1,1);
    
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
        ve(:,1)=ve(:,1)+ae(:,1)*ddt;
        position_e=position_e+ve(:,1)*ddt;
        
        %     %蛙跳法+磁场矢量分析更新电子位置
        %     ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ration_e*simulation.B0*ve(:,3))*ddt;
        %     ve(:,3)=ve(:,3)+constants.q_m_ration_e*simulation.B0*ddt*ve(:,1);% By
        %     position_e=position_e+ve(:,1)*ddt;
        
        %     vH(:,1)=vH(:,1)+aH(:,1)*ddt;%蛙跳法更新离子位置
        %     position_H=position_H+vH(:,1)*ddt;
        
        
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
        
        %------粒子MCC处理------------------------------------------------------------------------------------------
        [collision_index,collision_type]=null_collision( ve,constants,mcc_para,Xsec,'H2' );
        [colli_result]=collision_process(Xsec, mcc_para,ve,constants,collision_index,collision_type ,'H2');
        
        if colli_result.happen==1
            ve(colli_result.ion_index,:)=colli_result.ion_inc_v;
            ve(colli_result.ela_index,:)=colli_result.ela_v;
            ve(colli_result.exc_index,:)=colli_result.exc_v;
            ve(colli_result.att_index,:)=[];
            ve=[ve; colli_result.ion_gen_v ];
            position_e=[position_e;position_e(colli_result.ion_index)];
            %         vH2=[vH2; normrnd(0,vHth,[length(colli_result.ion_gen_v) , 3]);  ];
        end
        
        %------粒子MCC处理结束------------------------------------------------------------------------------------------
        
        %-----------分配电荷----------------------------------------------------------------------------------
        %     rhoea=zeros(simulation.num_grid_point,1);%网格左右两个点的粒子密度 电子
        %     rhoeb=zeros(simulation.num_grid_point,1);
        %     rhoaH=zeros(simulation.num_grid_point,1);%离子
        %     rhobH=zeros(simulation.num_grid_point,1);
        %
        %     enearx=floor(position_e/ddx)+1;%每个电子所在网格的索引 整数格点
        %     eassigndx=position_e-(enearx-1)*ddx;%每个电子距离左边最近的网格的距离
        %     Eenearx=floor((position_e+0.5*ddx)/ddx)+1;%每个电子所在半网格的索引 半格点
        %     Eedx=position_e-(Eenearx*ddx-3/2*ddx);%每个电子距离左边最近的半网格的距离
        %
        %     Hnearx=floor(position_H/ddx)+1;
        %     Hassigndx=position_H-(Hnearx-1)*ddx;
        %     EHnearx=floor((position_H+0.5*ddx)/ddx)+1;%每个离子所在半网格的索引 半格点
        %     EHdx=position_H-(EHnearx*ddx-3/2*ddx);
        %
        %     rhoea=accumarray(enearx,ddx-eassigndx,[simulation.num_grid_point 1]);
        %     rhoeb=accumarray(enearx+1,eassigndx,[simulation.num_grid_point 1]);
        %
        %     rhoe=(rhoea+rhoeb)*(-constants.e)/(ddx^2)*simulation.weight;
        %     %rhoP=rhoP*0;
        %     rhoaH=accumarray(Hnearx,ddx-Hassigndx,[simulation.num_grid_point 1]);
        %     rhobH=accumarray(Hnearx+1,Hassigndx,[simulation.num_grid_point 1]);
        %
        %     rhoH=(rhoaH+rhobH)*constants.e/(ddx^2)*simulation.weight;
        %
        %     rho=rhoe+rhoH;
        %     %-----------分配电荷结束------------------------------------------------------------
        %
        %     %--------------求解电场---------------------------------------------------------------
        %     b=get_b( rho,b_extra,simulation );
        %     u=get_u( u,A,b,'direct inverse');
        %     E=get_E_at_half_node( u, simulation );
        
        %--------------求解电场结束---------------------------------------------------------------
        
        % TODO：待获得E_particle
        ae=-E*constants.e/(constants.me);
        %     aH=((ddx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/ddx;
        
        % 电压累加，为  实时诊断中取平均  做准备
        %     if rem(ti+j_average_counter,display_timesteps)==1 %ti处于display的前average_timesteps范围中
        %         usum=usum+u;
        %         j_average_counter=j_average_counter-1;
        %         if j_average_counter==-1 %ti已超出display的前average_timesteps范围
        %             j_average_counter=average_timesteps; %重置j_average_counter
        %         end
        %     end
        
        %-----------实时诊断----------------------------------------------------------------------------------
        if rem(ti-1,display_timesteps)==0%每隔dptime显示一次
            fprintf('当前总时步数：%d/%d \n',Tdi,length(Td_list))
            fprintf('当前小时步数：%d/%d \n',ti,simulation.all_timesteps)
            fprintf('当前时步数漂移速度：%d m/s \n',mean(ve(:,1)))
            
%             display_i=ceil(ti/display_timesteps);
%             e_num(display_i)=size(position_e,1);%区域内的电子总量变化
%             H_num(display_i)=size(position_H,1);
%             e_vAve(display_i)=mean(abs(ve(:,1)));%x方向上的平均速率变化
%             H_vAve(display_i)=mean(abs(vH(:,1)));
%             
%             figure(h_fig1)
%             subplot(3,2,1);
%             %         plot_no_versus_position( position_e,position_H,simulation );
%             subplot(3,2,2);
%             plot_num_versus_timestep( e_num,H_num,display_timesteps )
%             subplot(3,2,3);
%             %         plot_2u1E_versus_x( u, ti, usum, average_timesteps, E, simulation )
%             subplot(3,2,4);
%             plot_Ek_versus_timestep( e_vAve,H_vAve, constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', display_timesteps )
%             ylabel('平均Ek [eV]'); %实际上是3倍x向平均动能
%             subplot(3,2,5);%电荷分布
%             %         plot_density_versus_x( rhoe, rhoH, rho, simulation )
%             %         subplot(3,2,6);
%             %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%             plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
%             drawnow;
%             
%             usum=0; %重置
        end
        %-----------实时诊断----------------------------------------------------------------------------------
        % 为主循环终止后的后处理做准备
        if ti==15000 %记录电压的时间步
            record_flag=1;
        end
        if record_flag>0 %从15000一直取到终止，大量时步取平均以降噪
            avg_drift_vex(record_flag,1)=mean(ve(:,1));
            record_flag=record_flag+1;
        end
        
        %终止判据：达到总时步数时终止
    end
    drift_velocity(Tdi)=mean(avg_drift_vex);
end

ref_v=[5.861 9.244 13.35 17.27 21.48 25.42 30.45 40.09]*1E4;%m/s
% ref_v=[5.861 9.244 13.35 17.27]*1E4;%m/s
figure
plot(Td_list,drift_velocity,'r')
hold on
scatter(Td_list,ref_v,'k')
xlabel('E/N(Td)') ; ylabel('electron drift velocity(m/s)');
title('electron drift velocity compare to experiment')
legend('MCC drift velocity','Experiment W Roznerski and K Leja 1984')
hold off

save ./examples/example_electron_drift_velocity/MCCdrift_velocity_testResults.mat Td_list drift_velocity ref_v;