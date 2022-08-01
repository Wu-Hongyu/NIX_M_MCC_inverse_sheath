%1D PIC
% 使用：根据物理模型，依次设置各类参数，进行相应初始化，
% 并在主循环中使用相应模块（尤其是pusher和粒子边界条件）

%% 初始化
clear
close all;

Hx=[0 1 10 100 500 1000];%约化磁场

Td=1;%约化电场
pressure=10;%Pa
n0=pressure/(1.381E-23*293);%气体密度 293K温度
E=-Td*1E-21*n0;%V/m



if exist('get_path_init','file')~=2
    addpath('../../packages') % 添加./packages到路径
end
path_list=get_path_init('now'); %文件路径
% test_all % 运行全部测试

for Hi=1:length(Hx)
    B=Hx(Hi)*1E-27*n0;%T
    % 如无说明，则单位为国际单位制
    % 全局变量
    constants=get_constants();% 全局常数 结构体
    
    %--------仿真参数-----------------------------------------------------------------------------------
    simulation=get_simulation('default'); %仿真参数 结构体
    % 可通过取消以下代码注释来修改仿真参数
    % % 等离子体参数，以电子为代表
    % simulation.ne0=1E16;%等离子体密度
    simulation.Te0=0;%单位eV
    % %         TH=1;%单位eV
    % % 时空尺度
    simulation.dt=100*10^-12;%时间步长
    simulation.end_time=2E-6;%仿真时间总长
    simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%总循环次数
    simulation.Lx=0.1;%仿真区域长度
    simulation.source_region=[0.499*simulation.Lx,0.501*simulation.Lx];
    simulation.num_grid_point=101;%网格格点数目
    simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%空间步长
    % % 电磁场
    % simulation.field_boundaries_type=0;%电位边界条件类型
    % simulation.field_boundaries=[0,0];%电位边界条件值
    simulation.B0=B; %代表性磁感应强度 T
    % Bx_posi=0.05;%位置
    % Bdelta=0.01;%形状因子
    % 粒子
    simulation.particle_boundaries_type=0;%粒子边界条件类型 0周期 1成对+热化  2成对恒注入+热化
    simulation.field_boundaries_type=0;%电位边界条件类型 0左右一类边界  3周期边界
    inject_num_pertime=2;%每次注入2个电子和1个H+和1个H2+
    source_density=0;
    source_left_bdidx=ceil((simulation.num_grid_point-1)*simulation.source_region(1)/simulation.Lx);
    source_right_bdidx=floor((simulation.num_grid_point-1)*simulation.source_region(2)/simulation.Lx);
    %气压
    simulation.pressure=pressure;%气压 单位Pa
    
    ddx=simulation.dx;
    ddt=simulation.dt;
    %--------粒子参数-----------------------------------------------------------------------------------
    % TODO：部分待加入particle_group
    simulation.THp=1;%单位eV
    simulation.TH2p=1;%单位eV
    simulation.num0_macro_e=100000;%初始电子数目 Particle e Num
    simulation.num0_macro_Hp=0.5*simulation.num0_macro_e;%初始离子数目Particle H Num
    simulation.num0_macro_H2p=0.5*simulation.num0_macro_e;%初始离子数目Particle H Num
    simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
    simulation.weight_Hp=simulation.weight;
    simulation.weight_H2p=simulation.weight;
    
    % TODO: 现在可以通过result=check_simulation_parameters返回值获得result.vth_e
    veth=sqrt(-constants.q_m_ratio_e*simulation.Te0);%电子温度对应的热速度
    vHth=sqrt(constants.q_m_ratio_Hp*simulation.THp);%离子温度对应的热速度
    
    vH2pth=sqrt(constants.q_m_ratio_H2p*simulation.TH2p);%离子温度对应的热速度
    
    %--------MCC参数-----------------------------------------------------------------------------------
    
    [mcc_para,Xsec]=MCC_init(simulation,constants,simulation.pressure,'1994Ness');
    
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
    if simulation.particle_boundaries_type==2
        simulation.num0_macro_e=0.01* simulation.num0_macro_e;
        simulation.num0_macro_Hp=0.01*simulation.num0_macro_Hp;
        simulation.num0_macro_H2p=0.01*simulation.num0_macro_H2p;
    end
    %按照Maxwell分布生成宏粒子的速度分布 x y z方向(Particle e velocity)
    ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
    vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_Hp,3]);
    vH2p=get_v_init( vH2pth, 'Maxwellian velocity', [simulation.num0_macro_H2p,3]);
    % ve=[1E6*ones(size(ve,1),1) zeros(size(ve,1),1) zeros(size(ve,1),1)];
    
    %%初始速度分布直接规定为-dt/2时刻的值，直接开始蛙跳法
    
    %均匀分布在空间Lx内，并给一个随机扰动
    position_e=get_position_init(simulation, 'source region uniform random', [simulation.num0_macro_e,1] );%source region uniform random
    position_H=get_position_init(simulation, 'source region uniform random', [simulation.num0_macro_Hp,1] );%entire domain uniform+noise
    position_H2p=get_position_init(simulation, 'source region uniform random', [simulation.num0_macro_H2p,1] );
    
    ae=zeros(simulation.num0_macro_e,1);%加速度初值 电子加速度
    aH=zeros(simulation.num0_macro_Hp,1);%加速度初值 离子加速度
    aH2p=zeros(simulation.num0_macro_H2p,1);%加速度初值 离子加速度
    
    if simulation.particle_boundaries_type==2
        simulation.num0_macro_e=100* simulation.num0_macro_e;
        simulation.num0_macro_Hp=100*simulation.num0_macro_Hp;
        simulation.num0_macro_H2p=100*simulation.num0_macro_H2p;
    end
    
    if simulation.particle_boundaries_type==2
        simulation.weight=simulation.ne0*(simulation.source_region(2)-simulation.source_region(1))/simulation.num0_macro_e;%每个宏粒子代表实际粒子的权重
        simulation.weight_Hp=simulation.weight;
        simulation.weight_H2p=simulation.weight;
    end
    
    
    % Bfun=@(position_e) (exp(-(position_e-Bx_posi).^2/(2*Bdelta^2)));
%     Bfun=@(position_e) 1;
    Be=simulation.B0*ones(simulation.num0_macro_e,1);
    
    %--------场数据
    [ u, A, b_extra ] = get_field_init( simulation );%离散泊松方程
    rho=zeros(simulation.num_grid_point,1);%净正电荷数密度初值
    % E=zeros(simulation.num_grid_point+1,1);%半节点电场初值  NG-1+两个墙内的半格点
    % B=zeros(simulation.num_grid_point+1,1);%半节点磁场初值
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
        %     % Pev=Pev+ae*dt;  %蛙跳法更新电子位置
        %     ve(:,1)=ve(:,1)+ae(:,1)*ddt;
        %     position_e=position_e+ve(:,1)*ddt;
        
        %     %蛙跳法+磁场矢量分析更新电子位置
        %     ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ration_e*simulation.B0*ve(:,3))*ddt;
        %     ve(:,3)=ve(:,3)+constants.q_m_ration_e*simulation.B0*ddt*ve(:,1);% By
        
        %%%%%%boris move%%%%%%
        vn(:,1)=ve(:,1)+ae*ddt/2;
        vn(:,3)=ve(:,3);
        t=constants.q_m_ratio_e*Be*ddt/2;
        s=2*t./(1+t.^2);
        
        vp(:,1)=vn(:,1)-(vn(:,3)+vn(:,1).*t).*s;
        vp(:,3)=vn(:,3)+(vn(:,1)-vn(:,3).*t).*s;
        ve(:,1)=vp(:,1)+ae*ddt/2;
        ve(:,3)=vp(:,3);
%         ve=vp+ae*ddt/2;
        clear vp
        %%%%%%%%%%%%%%%%%%
        
        position_e=position_e+ve(:,1)*ddt;
        
        %     vH(:,1)=vH(:,1)+aH(:,1)*ddt;%蛙跳法更新离子位置
        %     position_H=position_H+vH(:,1)*ddt;
        %
        %     vH2p(:,1)=vH2p(:,1)+aH2p(:,1)*ddt;%蛙跳法更新离子位置
        %     position_H2p=position_H2p+vH2p(:,1)*ddt;
        
        %------推动粒子结束-----------------------------------------------------------------------------------------------
        
        %------粒子越界处理------------------------------------------------------------------------------------------------
        switch simulation.particle_boundaries_type
            case 0 % 周期边界
                % -------周期边界------------------------------------------------------
                tempFlag=(position_e<0);
                position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
                tempFlag=(position_e>simulation.Lx);
                position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
                %             tempFlag=(position_H<0);
                %             position_H(tempFlag)=position_H(tempFlag)+simulation.Lx;
                %             tempFlag=(position_H>simulation.Lx);
                %             position_H(tempFlag)=position_H(tempFlag)-simulation.Lx;
                %             tempFlag=(position_H2p<0);
                %             position_H2p(tempFlag)=position_H2p(tempFlag)+simulation.Lx;
                %             tempFlag=(position_H2p>simulation.Lx);
                %             position_H2p(tempFlag)=position_H2p(tempFlag)-simulation.Lx;
            case 1 % 成对重注入+热化
                % -------成对重注入源区 pair re-injection into source region--------------
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%抹除超出区域的电子
                position_e(tempFlag)=[];
                
                k1=(position_H<0|position_H>simulation.Lx);%取超出范围的离子的编号
                tempsum1=sum(k1);
                
                if(tempsum1>0)%如果非空
                    vH(k1,:)=normrnd(0,vHth,[tempsum1,3]);%xyz 重注入离子速度Maxwell分布
                    % 重注入粒子位置均匀随机分布在源区内
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum1,1);
                    position_H(k1)=temp_matrix;
                    position_e=[position_e; temp_matrix];%合并回原来的数组
                    ve=[ve; normrnd(0,veth,[tempsum1,3])];
                end
                
                k2=(position_H2p<0|position_H2p>simulation.Lx);%取超出范围的离子的编号
                tempsum2=sum(k2);
                
                if(tempsum2>0)%如果非空
                    %                 if size(vH2p,1)<=simulation.num0_macro_H2p%%保持H2+的数量不增加
                    vH2p(k2,:)=normrnd(0,vH2pth,[tempsum2,3]);%xyz 重注入离子速度Maxwell分布
                    % 重注入粒子位置均匀随机分布在源区内
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum2,1);
                    position_H2p(k2)=temp_matrix;
                    position_e=[position_e; temp_matrix];%合并回原来的数组
                    ve=[ve; normrnd(0,veth,[tempsum2,3])];
                    %                 else
                    %                     vH2p(k2,:)=[];%抹除超出区域的H2+
                    %                     position_H2p(k2)=[];
                    %                 end
                end
                
                %--------thermalization 热化  每10步一次
                % 搭配非粒子周期边界使用
                if rem(ti,10)==1
                    flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                    ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
                end
                
            case 2 %成对恒注入
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%抹除超出区域的电子
                position_e(tempFlag)=[];
                
                tempFlag2=(position_H<0|position_H>simulation.Lx);
                vH(tempFlag2,:)=[];%抹除超出区域的H+
                position_H(tempFlag2)=[];
                
                tempFlag3=(position_H2p<0|position_H2p>simulation.Lx);
                vH2p(tempFlag3,:)=[];%抹除超出区域的H2+
                position_H2p(tempFlag3)=[];
                
                source_density=0;
                
                if source_density<(simulation.ne0*1.5)
                    if rem(ti,10)==1
                        temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime,1);
                        position_e=[position_e; temp_matrix];%合并回原来的数组
                        ve=[ve; normrnd(0,veth,[inject_num_pertime,3])];
                        
                        %                     temp_matrix1=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime/2,1);
                        position_H=[position_H; temp_matrix(1)];%合并回原来的数组
                        vH=[vH; normrnd(0,vHth,[inject_num_pertime/2,3])];
                        
                        %                     temp_matrix2=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime/2,1);
                        position_H2p=[position_H2p; temp_matrix(2)];%合并回原来的数组
                        vH2p=[vH2p; normrnd(0,vH2pth,[inject_num_pertime/2,3])];
                    end
                end
                
                %--------thermalization 热化  每10步一次
                % 搭配非粒子周期边界使用
                if rem(ti,10)==1
                    flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                    ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
                end
                
            otherwise
                throw(get_exception('error','no such type of particle boundaries.'))
        end
        %------粒子越界处理结束------------------------------------------------------------------------------------------
        
        %------粒子MCC处理------------------------------------------------------------------------------------------
        
        [collision_index,collision_type]=null_collision( ve,constants,mcc_para,Xsec,'H2' );
        [colli_result]=collision_process(Xsec, mcc_para,ve,constants,collision_index,collision_type,'H2' );
        
        if colli_result.happen==1
            ve(colli_result.ion_index,:)=colli_result.ion_inc_v;
            ve(colli_result.ela_index,:)=colli_result.ela_v;
            ve(colli_result.exc_index,:)=colli_result.exc_v;
            %             ve(colli_result.att_index,:)=[];
            %         ve=[ve; colli_result.ion_gen_v ];
            %         position_e=[position_e; position_e(colli_result.ion_index)];
            %
            %         vH2p=[vH2p; normrnd(0,vH2pth,[size(colli_result.ion_gen_v,1) , 3]);  ];
            %         position_H2p=[position_H2p; position_e(colli_result.ion_index)];
        end
        
        
        %------粒子MCC处理结束------------------------------------------------------------------------------------------
        
        %-----------分配电荷----------------------------------------------------------------------------------
        
        %     enearx=floor(position_e/ddx)+1;%每个电子所在网格的索引 整数格点
        %     eassigndx=position_e-(enearx-1)*ddx;%每个电子距离左边最近的网格的距离
        %     Eenearx=floor((position_e+0.5*ddx)/ddx)+1;%每个电子所在半网格的索引 半格点
        %     Eedx=position_e-(Eenearx*ddx-3/2*ddx);%每个电子距离左边最近的半网格的距离
        %
        %     Hnearx=floor(position_H/ddx)+1; %H+
        %     Hassigndx=position_H-(Hnearx-1)*ddx;
        %     EHnearx=floor((position_H+0.5*ddx)/ddx)+1;%每个离子所在半网格的索引 半格点
        %     EHdx=position_H-(EHnearx*ddx-3/2*ddx);
        %
        %     H2nearx=floor(position_H2p/ddx)+1; %H2+
        %     H2assigndx=position_H2p-(H2nearx-1)*ddx;
        %     EH2nearx=floor((position_H2p+0.5*ddx)/ddx)+1;%每个离子所在半网格的索引 半格点
        %     EH2dx=position_H2p-(EH2nearx*ddx-3/2*ddx);
        %
        %     rhoea=accumarray(enearx,ddx-eassigndx,[simulation.num_grid_point 1]);
        %     rhoeb=accumarray(enearx+1,eassigndx,[simulation.num_grid_point 1]);
        %
        %     rhoe=(rhoea+rhoeb)*(-constants.e)/(ddx^2)*simulation.weight;
        %     %rhoP=rhoP*0;
        %     rhoaH=accumarray(Hnearx,ddx-Hassigndx,[simulation.num_grid_point 1]);
        %     rhobH=accumarray(Hnearx+1,Hassigndx,[simulation.num_grid_point 1]);
        %
        %     rhoH=(rhoaH+rhobH)*constants.e/(ddx^2)*simulation.weight_Hp;
        %
        %     rhoaH2p=accumarray(H2nearx,ddx-H2assigndx,[simulation.num_grid_point 1]);
        %     rhobH2p=accumarray(H2nearx+1,H2assigndx,[simulation.num_grid_point 1]);
        %
        %     rhoH2p=(rhoaH2p+rhobH2p)*constants.e/(ddx^2)*simulation.weight_H2p;
        %
        %     rho=rhoe+rhoH+rhoH2p;
        %     %-----------分配电荷结束------------------------------------------------------------
        %
        %     source_density=-mean(rhoe(source_left_bdidx:source_right_bdidx))/constants.e;
        %
        %     %--------------求解电场---------------------------------------------------------------
        %     b=get_b( rho,b_extra,simulation );
        %     u=get_u( u,A,b,'direct inverse');
        %     E=get_E_at_half_node( u, simulation );
        %     E=-1E5*zeros(size(E,1),1);
        %--------------求解电场结束---------------------------------------------------------------
        
        % TODO：待获得E_particle
        %     ae=-((ddx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/ddx;
        %     aH=((ddx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/ddx;
        %     aH2p=((ddx-EH2dx).*E(EH2nearx)+EH2dx.*E(EH2nearx+1))*constants.e/(constants.mH2)/ddx;
        ae=-E*constants.e/(constants.me);
        %     Be=(ddx-Eedx).*B(Eenearx)+Eedx.*B(Eenearx+1);%应该用表达式来写
        Be=simulation.B0*ones(simulation.num0_macro_e,1);
        
        
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
            fprintf('电场方向速度:%d     ExB方向速度:%d \n',mean(ve(:,1)),mean(ve(:,3)))
            
            display_i=ceil(ti/display_timesteps);
            e_num(display_i)=size(position_e,1);%区域内的电子总量变化
            H_num(display_i)=size(position_H,1);
            e_vAve(display_i)=mean(abs(ve(:,1)));%x方向上的平均速率变化
            H_vAve(display_i)=mean(abs(vH(:,1)));
            
            figure(h_fig1)
            subplot(3,2,1);
            hold on
            scatter(ti,mean(ve(:,1)),'b')
            scatter(ti,mean(ve(:,3)),'r')
            xlabel('时间步数');
            ylabel('迁移速度m/s');
            legend('E方向','ExB方向')
            
            %         plot_no_versus_position( position_e(1:floor(size(position_e,1)/1000):end),position_H(1:floor(size(position_H,1)/1000):end),position_H2p(1:floor(size(position_H2p,1)/1000):end),simulation );
            %         plot_no_versus_position( position_e,position_H,position_H2p,simulation );
            subplot(3,2,2);
            plot_num_versus_timestep( e_num,H_num,[],[],display_timesteps )
            subplot(3,2,3);
            plot_2u1E_versus_x( u, ti, usum, average_timesteps, E, simulation )
            subplot(3,2,4);
%             plot_Ek_versus_timestep( e_vAve,H_vAve, constants.q_m_ratio_e, constants.q_m_ratio_Hp, 'e','H', display_timesteps )
%             ylabel('平均Ek [eV]'); %实际上是3倍x向平均动能
%             subplot(3,2,5);%电荷分布
            %         plot_density_versus_x( rhoe, rhoH, rho, simulation )
%             subplot(3,2,6);
            %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%             plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
            drawnow;
            
            usum=0; %重置
        end
        % if rem(ti-1,10)==0
        %     figure(H2)
        %     scatter(position_e(1:floor(size(position_e,1)/20):end),ti*ones(size(position_e(1:floor(size(position_e,1)/20):end),1),1),1,1:size(position_e(1:floor(size(position_e,1)/20):end),1))
        %     hold on
        %     drawnow;
        % end
        
        %-----------实时诊断----------------------------------------------------------------------------------
        % 为主循环终止后的后处理做准备
        if ti==floor(0.95*simulation.all_timesteps) %记录电压的时间步
            record_flag=1;
        end
        if record_flag>0 %从500000一直取到终止，大量时步取平均以降噪
            u_record(:,record_flag)=u;
            avg_ve1(record_flag)=mean(ve(:,1));
            avg_ve3(record_flag)=mean(ve(:,3));
            record_flag=record_flag+1;
        end
        
        %终止判据：达到总时步数时终止
    end
    %----主循环结束-----------------------------------------------------------------
    Wz(Hi)=mean(avg_ve1);
    Wx(Hi)=-mean(avg_ve3);
end
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp([now_str ' 主循环结束'])
diary off
now_str=datestr(now,'HH_MM');

% 终止时实时诊断图像（过程数据的示意）存储
figure(h_fig1)
saveas(gcf,[path_list.output '/终止时实时诊断' now_str '.png'])

Wz_ness=[4.9428 4.9427 4.931 3.972 0.6954 0.19436]*1000;
Wx_ness=[0 0.02443 0.2437 1.963 1.7186 0.96067]*1000;
figure
scatter(1:length(Wz),Wz)
set(gca,'xTickLabel', [0 1 10 100 500 1000])
hold on
scatter(1:length(Wz),Wz_ness)
grid on
xlabel('约化磁场Hx(B/n0)');
ylabel('E方向迁移速度m/s');
legend('Nix-M','Ness1994')
hold off


figure
scatter(1:length(Wx),Wx)
set(gca,'xTickLabel', [0 1 10 100 500 1000])
hold on
scatter(1:length(Wx),Wx_ness)
grid on
xlabel('约化磁场Hx(B/n0)');
ylabel('ExB方向迁移速度m/s');
legend('Nix-M','Ness1994')
hold off

% 终止（人工判断稳定）时电势空间分布
% figure
% x_final=(0:1:(simulation.num_grid_point-1))*ddx;
% u_final=mean(u_record,2);
% figure
% plot(x_final,u_final,'-b','LineWidth',3)%多个步长电压平均值
% axis([0,simulation.Lx,-inf,inf])
% %         title('电势空间分布', 'FontSize', 18);
% xlabel('x [m]')
% ylabel('\phi [V]');
% hold on;
% line([simulation.Lx-6.66e-4,simulation.Lx],[2.51,2.51],'linestyle','-.','color','r','LineWidth',3);
% L1=legend('仿真','解析');
% set(L1,'location','south');
% set(L1,'AutoUpdate','off');
% line([simulation.Lx-6.66e-4,simulation.Lx-6.66e-4],[0,2.51],'linestyle','-.','color','r','LineWidth',3);
% saveas(gcf,[path_list.output '/终止时电势分布' now_str '.png'])

%终止时全部数据存储
% save([path_list.output '/final_data' now_str '.mat'])

%% 后处理
% 读取存储过程数据