%��Լ���糡�µĵ���Ư�����ʲ���
%��һ���ͷţ�Ȼ�󾭹�һ��ʱ���ƽ�����е��ӵ��ٶȵõ�Ư���ٶ�

%��Ҫ�޸ģ�
%1.�糡��Լ���糡�õ�
% 2.�����ܶ���Լ���糡�õ�
% 3.��������Ӧ�޸Ľϴ�  ������250Td��Լ���糡�У����ӵ�Ư���ٶ�Ϊ40E4m/s  �����ֲ��øģ�ʵ���Ư�Ƴ���Ϊ10cm��
% 4.������Ҫ�޸�  ��һ��ʱ�䲽���ڵķ�Ӧ����ԶС��1

% ����һ��Լ���糡�µ�����


if exist('get_path_init','file')~=2
    addpath('../../packages') % ���./packages��·��
end
path_list=get_path_init('drift_velocity'); %�ļ�·��

Td_list=[50 75 100 125 150 175 200 250];
% Td_list=[50 75 100 125];
constants=get_constants();% ȫ�ֳ��� �ṹ��
simulation=get_simulation('default'); %������� �ṹ��
simulation.num0_macro_e=20000;%��ʼ������Ŀ Particle e Num
EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%�ٶȱ�����
% Te=[ 10:10:50];
% Te=logspace(log10(0.1),log10(50),50);
% Te=50;

% %         TH=1;%��λeV
% % ʱ�ճ߶�
simulation.dt=1*10^-12;%ʱ�䲽��
simulation.end_time=2E-8;%����ʱ���ܳ�
simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
% simulation.Lx=0.01;%�������򳤶�
% simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
% simulation.num_grid_point=201;%��������Ŀ
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
% % ��ų�
% simulation.field_boundaries_type=0;%��λ�߽���������
% simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
% simulation.B0=0; %�����ԴŸ�Ӧǿ�� T
% ����
simulation.particle_boundaries_type=0;%���ӱ߽���������
simulation.field_boundaries_type=3;%��λ�߽���������

ddx=simulation.dx;
ddt=simulation.dt;
%--------���Ӳ���-----------------------------------------------------------------------------------
% TODO�����ִ�����particle_group
simulation.TH0=1;%��λeV
simulation.num0_macro_H=simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
% TODO: ���ڿ���ͨ��result=check_simulation_parameters����ֵ���result.vth_e
veth=sqrt(-constants.q_m_ratio_e*simulation.Te0);%�����¶ȶ�Ӧ�����ٶ�
vHth=sqrt(constants.q_m_ratio_Hp*simulation.TH0);%�����¶ȶ�Ӧ�����ٶ�

for Tdi=1:length(Td_list)
    pressure=400;%Pa
    Td=Td_list(Tdi);
    
    %--------MCC����-----------------------------------------------------------------------------------
%     [mcc_para,Xsec]=MCC_init(simulation,constants,pressure);
   [mcc_para,Xsec]= MCC_init(simulation,constants,simulation.pressure, '2020LZS');
    E=-Td/(1E21)*mcc_para.n_target;
    
    %--------��ϲ���---------------------------------------------------------------------------
    display_timesteps=1000;%ÿ���ٲ���������ʾ
    average_timesteps=499;%���ٲ����ڵ�ƽ��ֵ
    
    %�����ʼ��־
    log_name = get_log_init( path_list.output, '' );
    diary(log_name) % �ض���ʽ�����־
    disp(simulation)
    result_default=check_simulation_parameters( simulation, 1 );
    disp('')
    fprintf('display_timesteps=%d ;%ÿ���ٲ���������ʾ\r\n',display_timesteps)
    fprintf('average_timesteps= %d ;%���ٲ����ڵ�ƽ��ֵ\r\n',average_timesteps)
    diary off
    
    %--------��ʼ��-------------------------------------------------------------------------------------
    % TODO�����ӽṹ�壬�볡�ṹ��
    %--------��������
    %����Maxwell�ֲ����ɺ����ӵ��ٶȷֲ� x y z����(Particle e velocity)
    ve=0*get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
    vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
    %%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������
    
    %���ȷֲ��ڿռ�Lx�ڣ�����һ������Ŷ�
    position_e=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_e,1] );
    position_H=0*get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_H,1] );
    position_e=0.01*position_e;%�ֲ���ǰ1%��λ��
    
    ae=zeros(simulation.num0_macro_e,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
    % aH=zeros(simulation.num0_macro_H,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
    %--------������
    [ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
    rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
    % E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
    B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
    %--------�������
    j_average_counter=average_timesteps;
    usum=0;
    record_flag=0;
    % u_record=zeros(simulation.num_grid_point,1);
    %--------��ʼ��--------------------------------------------------------------------------------------
    
%     h_fig1 = figure('Unit','Normalized','position',...
%         [0.02 0.3 0.6 0.6]); %ʵʱ���ʹ�õĴ󴰿�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    collision_energy_count=zeros(1,1);
    
    %--------��ѭ����ʼ-----------------------------------------------------------------------------------
    
    for ti=1:simulation.all_timesteps%��ѭ��
        % for ti=1:30000%��ѭ��
        
        % ��ѭ�������ܿ���
        % 1. �򵥴������inline����ͬ�����������ļ����Ա���test���ɽ�test����븴���޸ġ�
        % 2. ʹ��ר�ú�����ͨ�����������Ǵ������Ψһ��ʶ����ע�⴫�κ�ʱ����ʹ��ȫ�ֱ�����
        % �������⴫ֵ���ƣ��ο� https://www.zhihu.com/question/50408548/answer/120840847
        % 3. ��package��test�ļ��У����Թ���������
        
        %------�ƶ�����-----------------------------------------------------------------------------------------------
        % Pev=Pev+ae*dt;  %���������µ���λ��
        ve(:,1)=ve(:,1)+ae(:,1)*ddt;
        position_e=position_e+ve(:,1)*ddt;
        
        %     %������+�ų�ʸ���������µ���λ��
        %     ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ration_e*simulation.B0*ve(:,3))*ddt;
        %     ve(:,3)=ve(:,3)+constants.q_m_ration_e*simulation.B0*ddt*ve(:,1);% By
        %     position_e=position_e+ve(:,1)*ddt;
        
        %     vH(:,1)=vH(:,1)+aH(:,1)*ddt;%��������������λ��
        %     position_H=position_H+vH(:,1)*ddt;
        
        
        %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
        
        %------����Խ�紦��------------------------------------------------------------------------------------------------
        switch simulation.particle_boundaries_type
            case 0 % ���ڱ߽�
                % -------���ڱ߽�------------------------------------------------------
                tempFlag=(position_e<0);
                position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
                tempFlag=(position_e>simulation.Lx);
                position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
                tempFlag=(position_H<0);
                position_H(tempFlag)=position_H(tempFlag)+simulation.Lx;
                tempFlag=(position_H>simulation.Lx);
                position_H(tempFlag)=position_H(tempFlag)-simulation.Lx;
            case 1 % �ɶ���ע��+�Ȼ�
                % -------�ɶ���ע��Դ�� pair re-injection into source region--------------
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%Ĩ����������ĵ���
                position_e(tempFlag)=[];
                
                k1=(position_H<0|position_H>simulation.Lx);%ȡ������Χ�����ӵı��
                tempsum=sum(k1);
                
                if(tempsum>0)%����ǿ�
                    vH(k1,:)=normrnd(0,vHth,[tempsum,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
                    % ��ע������λ�þ�������ֲ���Դ����
                    position_H(k1)=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
                    position_e=reshape([position_e; temp_matrix],[],1);%�ϲ���ԭ��������
                    ve=reshape([ve; normrnd(0,veth,[tempsum,3])],[],3);
                end
                
                %--------thermalization �Ȼ�  ÿ50��һ��
                % ������������ڱ߽�ʹ��
                if rem(ti,50)==1
                    flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                    ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
                end
            otherwise
                throw(get_exception('error','no such type of particle boundaries.'))
        end
        %------����Խ�紦�����------------------------------------------------------------------------------------------
        
        %------����MCC����------------------------------------------------------------------------------------------
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
        
        %------����MCC�������------------------------------------------------------------------------------------------
        
        %-----------������----------------------------------------------------------------------------------
        %     rhoea=zeros(simulation.num_grid_point,1);%��������������������ܶ� ����
        %     rhoeb=zeros(simulation.num_grid_point,1);
        %     rhoaH=zeros(simulation.num_grid_point,1);%����
        %     rhobH=zeros(simulation.num_grid_point,1);
        %
        %     enearx=floor(position_e/ddx)+1;%ÿ������������������� �������
        %     eassigndx=position_e-(enearx-1)*ddx;%ÿ�����Ӿ���������������ľ���
        %     Eenearx=floor((position_e+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
        %     Eedx=position_e-(Eenearx*ddx-3/2*ddx);%ÿ�����Ӿ����������İ�����ľ���
        %
        %     Hnearx=floor(position_H/ddx)+1;
        %     Hassigndx=position_H-(Hnearx-1)*ddx;
        %     EHnearx=floor((position_H+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
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
        %     %-----------�����ɽ���------------------------------------------------------------
        %
        %     %--------------���糡---------------------------------------------------------------
        %     b=get_b( rho,b_extra,simulation );
        %     u=get_u( u,A,b,'direct inverse');
        %     E=get_E_at_half_node( u, simulation );
        
        %--------------���糡����---------------------------------------------------------------
        
        % TODO�������E_particle
        ae=-E*constants.e/(constants.me);
        %     aH=((ddx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/ddx;
        
        % ��ѹ�ۼӣ�Ϊ  ʵʱ�����ȡƽ��  ��׼��
        %     if rem(ti+j_average_counter,display_timesteps)==1 %ti����display��ǰaverage_timesteps��Χ��
        %         usum=usum+u;
        %         j_average_counter=j_average_counter-1;
        %         if j_average_counter==-1 %ti�ѳ���display��ǰaverage_timesteps��Χ
        %             j_average_counter=average_timesteps; %����j_average_counter
        %         end
        %     end
        
        %-----------ʵʱ���----------------------------------------------------------------------------------
        if rem(ti-1,display_timesteps)==0%ÿ��dptime��ʾһ��
            fprintf('��ǰ��ʱ������%d/%d \n',Tdi,length(Td_list))
            fprintf('��ǰСʱ������%d/%d \n',ti,simulation.all_timesteps)
            fprintf('��ǰʱ����Ư���ٶȣ�%d m/s \n',mean(ve(:,1)))
            
%             display_i=ceil(ti/display_timesteps);
%             e_num(display_i)=size(position_e,1);%�����ڵĵ��������仯
%             H_num(display_i)=size(position_H,1);
%             e_vAve(display_i)=mean(abs(ve(:,1)));%x�����ϵ�ƽ�����ʱ仯
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
%             ylabel('ƽ��Ek [eV]'); %ʵ������3��x��ƽ������
%             subplot(3,2,5);%��ɷֲ�
%             %         plot_density_versus_x( rhoe, rhoH, rho, simulation )
%             %         subplot(3,2,6);
%             %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%             plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
%             drawnow;
%             
%             usum=0; %����
        end
        %-----------ʵʱ���----------------------------------------------------------------------------------
        % Ϊ��ѭ����ֹ��ĺ�����׼��
        if ti==15000 %��¼��ѹ��ʱ�䲽
            record_flag=1;
        end
        if record_flag>0 %��15000һֱȡ����ֹ������ʱ��ȡƽ���Խ���
            avg_drift_vex(record_flag,1)=mean(ve(:,1));
            record_flag=record_flag+1;
        end
        
        %��ֹ�оݣ��ﵽ��ʱ����ʱ��ֹ
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