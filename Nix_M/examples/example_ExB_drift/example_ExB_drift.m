%1D PIC
% ʹ�ã���������ģ�ͣ��������ø��������������Ӧ��ʼ����
% ������ѭ����ʹ����Ӧģ�飨������pusher�����ӱ߽�������

%% ��ʼ��
clear
close all;

Hx=[0 1 10 100 500 1000];%Լ���ų�

Td=1;%Լ���糡
pressure=10;%Pa
n0=pressure/(1.381E-23*293);%�����ܶ� 293K�¶�
E=-Td*1E-21*n0;%V/m



if exist('get_path_init','file')~=2
    addpath('../../packages') % ���./packages��·��
end
path_list=get_path_init('now'); %�ļ�·��
% test_all % ����ȫ������

for Hi=1:length(Hx)
    B=Hx(Hi)*1E-27*n0;%T
    % ����˵������λΪ���ʵ�λ��
    % ȫ�ֱ���
    constants=get_constants();% ȫ�ֳ��� �ṹ��
    
    %--------�������-----------------------------------------------------------------------------------
    simulation=get_simulation('default'); %������� �ṹ��
    % ��ͨ��ȡ�����´���ע�����޸ķ������
    % % ��������������Ե���Ϊ����
    % simulation.ne0=1E16;%���������ܶ�
    simulation.Te0=0;%��λeV
    % %         TH=1;%��λeV
    % % ʱ�ճ߶�
    simulation.dt=100*10^-12;%ʱ�䲽��
    simulation.end_time=2E-6;%����ʱ���ܳ�
    simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
    simulation.Lx=0.1;%�������򳤶�
    simulation.source_region=[0.499*simulation.Lx,0.501*simulation.Lx];
    simulation.num_grid_point=101;%��������Ŀ
    simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
    % % ��ų�
    % simulation.field_boundaries_type=0;%��λ�߽���������
    % simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
    simulation.B0=B; %�����ԴŸ�Ӧǿ�� T
    % Bx_posi=0.05;%λ��
    % Bdelta=0.01;%��״����
    % ����
    simulation.particle_boundaries_type=0;%���ӱ߽��������� 0���� 1�ɶ�+�Ȼ�  2�ɶԺ�ע��+�Ȼ�
    simulation.field_boundaries_type=0;%��λ�߽��������� 0����һ��߽�  3���ڱ߽�
    inject_num_pertime=2;%ÿ��ע��2�����Ӻ�1��H+��1��H2+
    source_density=0;
    source_left_bdidx=ceil((simulation.num_grid_point-1)*simulation.source_region(1)/simulation.Lx);
    source_right_bdidx=floor((simulation.num_grid_point-1)*simulation.source_region(2)/simulation.Lx);
    %��ѹ
    simulation.pressure=pressure;%��ѹ ��λPa
    
    ddx=simulation.dx;
    ddt=simulation.dt;
    %--------���Ӳ���-----------------------------------------------------------------------------------
    % TODO�����ִ�����particle_group
    simulation.THp=1;%��λeV
    simulation.TH2p=1;%��λeV
    simulation.num0_macro_e=100000;%��ʼ������Ŀ Particle e Num
    simulation.num0_macro_Hp=0.5*simulation.num0_macro_e;%��ʼ������ĿParticle H Num
    simulation.num0_macro_H2p=0.5*simulation.num0_macro_e;%��ʼ������ĿParticle H Num
    simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
    simulation.weight_Hp=simulation.weight;
    simulation.weight_H2p=simulation.weight;
    
    % TODO: ���ڿ���ͨ��result=check_simulation_parameters����ֵ���result.vth_e
    veth=sqrt(-constants.q_m_ratio_e*simulation.Te0);%�����¶ȶ�Ӧ�����ٶ�
    vHth=sqrt(constants.q_m_ratio_Hp*simulation.THp);%�����¶ȶ�Ӧ�����ٶ�
    
    vH2pth=sqrt(constants.q_m_ratio_H2p*simulation.TH2p);%�����¶ȶ�Ӧ�����ٶ�
    
    %--------MCC����-----------------------------------------------------------------------------------
    
    [mcc_para,Xsec]=MCC_init(simulation,constants,simulation.pressure,'1994Ness');
    
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
    if simulation.particle_boundaries_type==2
        simulation.num0_macro_e=0.01* simulation.num0_macro_e;
        simulation.num0_macro_Hp=0.01*simulation.num0_macro_Hp;
        simulation.num0_macro_H2p=0.01*simulation.num0_macro_H2p;
    end
    %����Maxwell�ֲ����ɺ����ӵ��ٶȷֲ� x y z����(Particle e velocity)
    ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
    vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_Hp,3]);
    vH2p=get_v_init( vH2pth, 'Maxwellian velocity', [simulation.num0_macro_H2p,3]);
    % ve=[1E6*ones(size(ve,1),1) zeros(size(ve,1),1) zeros(size(ve,1),1)];
    
    %%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������
    
    %���ȷֲ��ڿռ�Lx�ڣ�����һ������Ŷ�
    position_e=get_position_init(simulation, 'source region uniform random', [simulation.num0_macro_e,1] );%source region uniform random
    position_H=get_position_init(simulation, 'source region uniform random', [simulation.num0_macro_Hp,1] );%entire domain uniform+noise
    position_H2p=get_position_init(simulation, 'source region uniform random', [simulation.num0_macro_H2p,1] );
    
    ae=zeros(simulation.num0_macro_e,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
    aH=zeros(simulation.num0_macro_Hp,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
    aH2p=zeros(simulation.num0_macro_H2p,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
    
    if simulation.particle_boundaries_type==2
        simulation.num0_macro_e=100* simulation.num0_macro_e;
        simulation.num0_macro_Hp=100*simulation.num0_macro_Hp;
        simulation.num0_macro_H2p=100*simulation.num0_macro_H2p;
    end
    
    if simulation.particle_boundaries_type==2
        simulation.weight=simulation.ne0*(simulation.source_region(2)-simulation.source_region(1))/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
        simulation.weight_Hp=simulation.weight;
        simulation.weight_H2p=simulation.weight;
    end
    
    
    % Bfun=@(position_e) (exp(-(position_e-Bx_posi).^2/(2*Bdelta^2)));
%     Bfun=@(position_e) 1;
    Be=simulation.B0*ones(simulation.num0_macro_e,1);
    
    %--------������
    [ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
    rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
    % E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
    % B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
    %--------�������
    j_average_counter=average_timesteps;
    usum=0;
    record_flag=0;
    u_record=zeros(simulation.num_grid_point,1);
    %--------��ʼ��--------------------------------------------------------------------------------------
    
    h_fig1 = figure('Unit','Normalized','position',...
        [0.02 0.3 0.6 0.6]); %ʵʱ���ʹ�õĴ󴰿�
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %--------��ѭ����ʼ-----------------------------------------------------------------------------------
    
    for ti=1:simulation.all_timesteps%��ѭ��
        % for ti=1:30000%��ѭ��
        
        % ��ѭ�������ܿ���
        % 1. �򵥴������inline����ͬ�����������ļ����Ա���test���ɽ�test����븴���޸ġ�
        % 2. ʹ��ר�ú�����ͨ�����������Ǵ������Ψһ��ʶ����ע�⴫�κ�ʱ����ʹ��ȫ�ֱ�����
        % �������⴫ֵ���ƣ��ο� https://www.zhihu.com/question/50408548/answer/120840847
        % 3. ��package��test�ļ��У����Թ���������
        
        %------�ƶ�����-----------------------------------------------------------------------------------------------
        %     % Pev=Pev+ae*dt;  %���������µ���λ��
        %     ve(:,1)=ve(:,1)+ae(:,1)*ddt;
        %     position_e=position_e+ve(:,1)*ddt;
        
        %     %������+�ų�ʸ���������µ���λ��
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
        
        %     vH(:,1)=vH(:,1)+aH(:,1)*ddt;%��������������λ��
        %     position_H=position_H+vH(:,1)*ddt;
        %
        %     vH2p(:,1)=vH2p(:,1)+aH2p(:,1)*ddt;%��������������λ��
        %     position_H2p=position_H2p+vH2p(:,1)*ddt;
        
        %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
        
        %------����Խ�紦��------------------------------------------------------------------------------------------------
        switch simulation.particle_boundaries_type
            case 0 % ���ڱ߽�
                % -------���ڱ߽�------------------------------------------------------
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
            case 1 % �ɶ���ע��+�Ȼ�
                % -------�ɶ���ע��Դ�� pair re-injection into source region--------------
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%Ĩ����������ĵ���
                position_e(tempFlag)=[];
                
                k1=(position_H<0|position_H>simulation.Lx);%ȡ������Χ�����ӵı��
                tempsum1=sum(k1);
                
                if(tempsum1>0)%����ǿ�
                    vH(k1,:)=normrnd(0,vHth,[tempsum1,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
                    % ��ע������λ�þ�������ֲ���Դ����
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum1,1);
                    position_H(k1)=temp_matrix;
                    position_e=[position_e; temp_matrix];%�ϲ���ԭ��������
                    ve=[ve; normrnd(0,veth,[tempsum1,3])];
                end
                
                k2=(position_H2p<0|position_H2p>simulation.Lx);%ȡ������Χ�����ӵı��
                tempsum2=sum(k2);
                
                if(tempsum2>0)%����ǿ�
                    %                 if size(vH2p,1)<=simulation.num0_macro_H2p%%����H2+������������
                    vH2p(k2,:)=normrnd(0,vH2pth,[tempsum2,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
                    % ��ע������λ�þ�������ֲ���Դ����
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum2,1);
                    position_H2p(k2)=temp_matrix;
                    position_e=[position_e; temp_matrix];%�ϲ���ԭ��������
                    ve=[ve; normrnd(0,veth,[tempsum2,3])];
                    %                 else
                    %                     vH2p(k2,:)=[];%Ĩ�����������H2+
                    %                     position_H2p(k2)=[];
                    %                 end
                end
                
                %--------thermalization �Ȼ�  ÿ10��һ��
                % ������������ڱ߽�ʹ��
                if rem(ti,10)==1
                    flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                    ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
                end
                
            case 2 %�ɶԺ�ע��
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%Ĩ����������ĵ���
                position_e(tempFlag)=[];
                
                tempFlag2=(position_H<0|position_H>simulation.Lx);
                vH(tempFlag2,:)=[];%Ĩ�����������H+
                position_H(tempFlag2)=[];
                
                tempFlag3=(position_H2p<0|position_H2p>simulation.Lx);
                vH2p(tempFlag3,:)=[];%Ĩ�����������H2+
                position_H2p(tempFlag3)=[];
                
                source_density=0;
                
                if source_density<(simulation.ne0*1.5)
                    if rem(ti,10)==1
                        temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime,1);
                        position_e=[position_e; temp_matrix];%�ϲ���ԭ��������
                        ve=[ve; normrnd(0,veth,[inject_num_pertime,3])];
                        
                        %                     temp_matrix1=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime/2,1);
                        position_H=[position_H; temp_matrix(1)];%�ϲ���ԭ��������
                        vH=[vH; normrnd(0,vHth,[inject_num_pertime/2,3])];
                        
                        %                     temp_matrix2=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime/2,1);
                        position_H2p=[position_H2p; temp_matrix(2)];%�ϲ���ԭ��������
                        vH2p=[vH2p; normrnd(0,vH2pth,[inject_num_pertime/2,3])];
                    end
                end
                
                %--------thermalization �Ȼ�  ÿ10��һ��
                % ������������ڱ߽�ʹ��
                if rem(ti,10)==1
                    flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                    ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
                end
                
            otherwise
                throw(get_exception('error','no such type of particle boundaries.'))
        end
        %------����Խ�紦�����------------------------------------------------------------------------------------------
        
        %------����MCC����------------------------------------------------------------------------------------------
        
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
        
        
        %------����MCC�������------------------------------------------------------------------------------------------
        
        %-----------������----------------------------------------------------------------------------------
        
        %     enearx=floor(position_e/ddx)+1;%ÿ������������������� �������
        %     eassigndx=position_e-(enearx-1)*ddx;%ÿ�����Ӿ���������������ľ���
        %     Eenearx=floor((position_e+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
        %     Eedx=position_e-(Eenearx*ddx-3/2*ddx);%ÿ�����Ӿ����������İ�����ľ���
        %
        %     Hnearx=floor(position_H/ddx)+1; %H+
        %     Hassigndx=position_H-(Hnearx-1)*ddx;
        %     EHnearx=floor((position_H+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
        %     EHdx=position_H-(EHnearx*ddx-3/2*ddx);
        %
        %     H2nearx=floor(position_H2p/ddx)+1; %H2+
        %     H2assigndx=position_H2p-(H2nearx-1)*ddx;
        %     EH2nearx=floor((position_H2p+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
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
        %     %-----------�����ɽ���------------------------------------------------------------
        %
        %     source_density=-mean(rhoe(source_left_bdidx:source_right_bdidx))/constants.e;
        %
        %     %--------------���糡---------------------------------------------------------------
        %     b=get_b( rho,b_extra,simulation );
        %     u=get_u( u,A,b,'direct inverse');
        %     E=get_E_at_half_node( u, simulation );
        %     E=-1E5*zeros(size(E,1),1);
        %--------------���糡����---------------------------------------------------------------
        
        % TODO�������E_particle
        %     ae=-((ddx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/ddx;
        %     aH=((ddx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/ddx;
        %     aH2p=((ddx-EH2dx).*E(EH2nearx)+EH2dx.*E(EH2nearx+1))*constants.e/(constants.mH2)/ddx;
        ae=-E*constants.e/(constants.me);
        %     Be=(ddx-Eedx).*B(Eenearx)+Eedx.*B(Eenearx+1);%Ӧ���ñ��ʽ��д
        Be=simulation.B0*ones(simulation.num0_macro_e,1);
        
        
        % ��ѹ�ۼӣ�Ϊ  ʵʱ�����ȡƽ��  ��׼��
        if rem(ti+j_average_counter,display_timesteps)==1 %ti����display��ǰaverage_timesteps��Χ��
            usum=usum+u;
            j_average_counter=j_average_counter-1;
            if j_average_counter==-1 %ti�ѳ���display��ǰaverage_timesteps��Χ
                j_average_counter=average_timesteps; %����j_average_counter
            end
        end
        
        %-----------ʵʱ���----------------------------------------------------------------------------------
        if rem(ti-1,display_timesteps)==0%ÿ��dptime��ʾһ��
            fprintf('��ǰʱ������%d/%d \n',ti,simulation.all_timesteps)
            fprintf('�糡�����ٶ�:%d     ExB�����ٶ�:%d \n',mean(ve(:,1)),mean(ve(:,3)))
            
            display_i=ceil(ti/display_timesteps);
            e_num(display_i)=size(position_e,1);%�����ڵĵ��������仯
            H_num(display_i)=size(position_H,1);
            e_vAve(display_i)=mean(abs(ve(:,1)));%x�����ϵ�ƽ�����ʱ仯
            H_vAve(display_i)=mean(abs(vH(:,1)));
            
            figure(h_fig1)
            subplot(3,2,1);
            hold on
            scatter(ti,mean(ve(:,1)),'b')
            scatter(ti,mean(ve(:,3)),'r')
            xlabel('ʱ�䲽��');
            ylabel('Ǩ���ٶ�m/s');
            legend('E����','ExB����')
            
            %         plot_no_versus_position( position_e(1:floor(size(position_e,1)/1000):end),position_H(1:floor(size(position_H,1)/1000):end),position_H2p(1:floor(size(position_H2p,1)/1000):end),simulation );
            %         plot_no_versus_position( position_e,position_H,position_H2p,simulation );
            subplot(3,2,2);
            plot_num_versus_timestep( e_num,H_num,[],[],display_timesteps )
            subplot(3,2,3);
            plot_2u1E_versus_x( u, ti, usum, average_timesteps, E, simulation )
            subplot(3,2,4);
%             plot_Ek_versus_timestep( e_vAve,H_vAve, constants.q_m_ratio_e, constants.q_m_ratio_Hp, 'e','H', display_timesteps )
%             ylabel('ƽ��Ek [eV]'); %ʵ������3��x��ƽ������
%             subplot(3,2,5);%��ɷֲ�
            %         plot_density_versus_x( rhoe, rhoH, rho, simulation )
%             subplot(3,2,6);
            %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%             plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
            drawnow;
            
            usum=0; %����
        end
        % if rem(ti-1,10)==0
        %     figure(H2)
        %     scatter(position_e(1:floor(size(position_e,1)/20):end),ti*ones(size(position_e(1:floor(size(position_e,1)/20):end),1),1),1,1:size(position_e(1:floor(size(position_e,1)/20):end),1))
        %     hold on
        %     drawnow;
        % end
        
        %-----------ʵʱ���----------------------------------------------------------------------------------
        % Ϊ��ѭ����ֹ��ĺ�����׼��
        if ti==floor(0.95*simulation.all_timesteps) %��¼��ѹ��ʱ�䲽
            record_flag=1;
        end
        if record_flag>0 %��500000һֱȡ����ֹ������ʱ��ȡƽ���Խ���
            u_record(:,record_flag)=u;
            avg_ve1(record_flag)=mean(ve(:,1));
            avg_ve3(record_flag)=mean(ve(:,3));
            record_flag=record_flag+1;
        end
        
        %��ֹ�оݣ��ﵽ��ʱ����ʱ��ֹ
    end
    %----��ѭ������-----------------------------------------------------------------
    Wz(Hi)=mean(avg_ve1);
    Wx(Hi)=-mean(avg_ve3);
end
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp([now_str ' ��ѭ������'])
diary off
now_str=datestr(now,'HH_MM');

% ��ֹʱʵʱ���ͼ�񣨹������ݵ�ʾ�⣩�洢
figure(h_fig1)
saveas(gcf,[path_list.output '/��ֹʱʵʱ���' now_str '.png'])

Wz_ness=[4.9428 4.9427 4.931 3.972 0.6954 0.19436]*1000;
Wx_ness=[0 0.02443 0.2437 1.963 1.7186 0.96067]*1000;
figure
scatter(1:length(Wz),Wz)
set(gca,'xTickLabel', [0 1 10 100 500 1000])
hold on
scatter(1:length(Wz),Wz_ness)
grid on
xlabel('Լ���ų�Hx(B/n0)');
ylabel('E����Ǩ���ٶ�m/s');
legend('Nix-M','Ness1994')
hold off


figure
scatter(1:length(Wx),Wx)
set(gca,'xTickLabel', [0 1 10 100 500 1000])
hold on
scatter(1:length(Wx),Wx_ness)
grid on
xlabel('Լ���ų�Hx(B/n0)');
ylabel('ExB����Ǩ���ٶ�m/s');
legend('Nix-M','Ness1994')
hold off

% ��ֹ���˹��ж��ȶ���ʱ���ƿռ�ֲ�
% figure
% x_final=(0:1:(simulation.num_grid_point-1))*ddx;
% u_final=mean(u_record,2);
% figure
% plot(x_final,u_final,'-b','LineWidth',3)%���������ѹƽ��ֵ
% axis([0,simulation.Lx,-inf,inf])
% %         title('���ƿռ�ֲ�', 'FontSize', 18);
% xlabel('x [m]')
% ylabel('\phi [V]');
% hold on;
% line([simulation.Lx-6.66e-4,simulation.Lx],[2.51,2.51],'linestyle','-.','color','r','LineWidth',3);
% L1=legend('����','����');
% set(L1,'location','south');
% set(L1,'AutoUpdate','off');
% line([simulation.Lx-6.66e-4,simulation.Lx-6.66e-4],[0,2.51],'linestyle','-.','color','r','LineWidth',3);
% saveas(gcf,[path_list.output '/��ֹʱ���Ʒֲ�' now_str '.png'])

%��ֹʱȫ�����ݴ洢
% save([path_list.output '/final_data' now_str '.mat'])

%% ����
% ��ȡ�洢��������