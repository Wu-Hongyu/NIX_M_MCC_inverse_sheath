%1D PIC
% ʹ�ã���������ģ�ͣ��������ø��������������Ӧ��ʼ����
% ������ѭ����ʹ����Ӧģ�飨������pusher�����ӱ߽�������

%% ��ʼ��
clear
close all;

if exist('get_path_init','file')~=2
    addpath('./packages') % ���./packages��·��
end
path_list=get_path_init('now'); %�ļ�·��
% test_all % ����ȫ������

% ����˵������λΪ���ʵ�λ��
% ȫ�ֱ���
constants=get_constants();% ȫ�ֳ��� �ṹ��

%--------�������-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %������� �ṹ��
% ��ͨ��ȡ�����´���ע�����޸ķ������
% % ��������������Ե���Ϊ����
% simulation.ne0=1E16;%���������ܶ�
simulation.Te0=1;%��λeV
% %         TH=1;%��λeV
% % ʱ�ճ߶�
% simulation.dt=10*10^-12;%ʱ�䲽��
% simulation.end_time=6E-6;%����ʱ���ܳ�
% simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
% simulation.Lx=0.01;%�������򳤶�
% simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
% simulation.num_grid_point=201;%��������Ŀ
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
% % ��ų�
simulation.field_boundaries_type=3;%��λ�߽���������
% simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
% simulation.B0=0; %�����ԴŸ�Ӧǿ�� T
% ����
simulation.particle_boundaries_type=0;%���ӱ߽���������

%--------���Ӳ���-----------------------------------------------------------------------------------
% TODO�����ִ�����particle_group
simulation.TH0=0;%��λeV
simulation.num0_macro_e=20000;%��ʼ������Ŀ Particle e Num
simulation.num0_macro_H=simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
% TODO: ���ڿ���ͨ��result=check_simulation_parameters����ֵ���result.vth_e
veth=sqrt(-constants.q_m_ration_e*simulation.Te0);%�����¶ȶ�Ӧ�����ٶ�
vHth=sqrt(constants.q_m_ration_H*simulation.TH0);%�����¶ȶ�Ӧ�����ٶ�

%--------��ϲ���---------------------------------------------------------------------------
display_timesteps=100;%ÿ���ٲ���������ʾ
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
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
%%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������

%���ȷֲ��ڿռ�Lx�ڣ�����һ������Ŷ�
position_e=get_position_init(simulation, 'entire domain uniform', [simulation.num0_macro_e,1] );
position_H=get_position_init(simulation, 'entire domain uniform', [simulation.num0_macro_H,1] );



ae=zeros(simulation.num0_macro_e,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
aH=zeros(simulation.num0_macro_H,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
%--------������
[ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
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
    % Pev=Pev+ae*dt;  %���������µ���λ��
    ve(:,1)=ve(:,1)+ae(:,1)*simulation.dt;
    position_e=position_e+ve(:,1)*simulation.dt;
    
    %     %������+�ų�ʸ���������µ���λ��
    %     ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ration_e*simulation.B0*ve(:,3))*simulation.dt;
    %     ve(:,3)=ve(:,3)+constants.q_m_ration_e*simulation.B0*simulation.dt*ve(:,1);% By
    %     position_e=position_e+ve(:,1)*simulation.dt;
    
    vH(:,1)=vH(:,1)+aH(:,1)*simulation.dt;%��������������λ��
    position_H=position_H+vH(:,1)*simulation.dt;
    
    
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
    
    %-----------������----------------------------------------------------------------------------------
    rhoea=zeros(simulation.num_grid_point,1);%��������������������ܶ� ����
    rhoeb=zeros(simulation.num_grid_point,1);
    rhoaH=zeros(simulation.num_grid_point,1);%����
    rhobH=zeros(simulation.num_grid_point,1);
    
    enearx=floor(position_e/simulation.dx)+1;%ÿ������������������� �������
    eassigndx=position_e-(enearx-1)*simulation.dx;%ÿ�����Ӿ���������������ľ���
    Eenearx=floor((position_e+0.5*simulation.dx)/simulation.dx)+1;%ÿ���������ڰ���������� ����
    Eedx=position_e-(Eenearx*simulation.dx-3/2*simulation.dx);%ÿ�����Ӿ����������İ�����ľ���
    
    Hnearx=floor(position_H/simulation.dx)+1;
    Hassigndx=position_H-(Hnearx-1)*simulation.dx;
    EHnearx=floor((position_H+0.5*simulation.dx)/simulation.dx)+1;%ÿ���������ڰ���������� ����
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
    %-----------�����ɽ���------------------------------------------------------------
    
    %--------------���糡---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------���糡����---------------------------------------------------------------
    
    % TODO�������E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    aH=((simulation.dx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/simulation.dx;
    
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
        
        display_i=ceil(ti/display_timesteps);
        e_num(display_i)=size(position_e,1);%�����ڵĵ��������仯
        H_num(display_i)=size(position_H,1);
        e_vAve(display_i)=mean(abs(ve(:,1)));%x�����ϵ�ƽ�����ʱ仯
        H_vAve(display_i)=mean(abs(vH(:,1)));
        
        
        % ����ƽ������
        Ek=@(v,q_m_ratio) sum(v.*v/(2*abs(q_m_ratio))); %���ܣ���λeV
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
%         ylabel('ƽ��Ek [eV]'); %ʵ������3��x��ƽ������

% �����غ�-л����  2020/11/05 17:04:51 ʧ�ܣ�Ϊʲô����ͬһ������
        % ��ͼ��ʾ �����ܶ�����糡�ܴ���-ʱ�䲽
        plot(ave_Ek-0.5,'b')
        hold on
        plot(ave_Ep,'r')
        plot(ave_Ek+ave_Ep-0.5,'k')
        xlabel(['t [' num2str(display_timesteps) 'ʱ��]']);
        ylabel('ƽ������');
        %         axis([0,inf,0,1.2])
        L1=legend('Ek','Ep','Ek+Ep');
        set(L1,'location','east');
        hold off



        subplot(3,2,5);%��ɷֲ�
        plot_density_versus_x( rhoe, rhoH, rho, simulation )
        subplot(3,2,6);
        %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
        plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
        drawnow;
        
        usum=0; %����
    end
    %-----------ʵʱ���----------------------------------------------------------------------------------
    % Ϊ��ѭ����ֹ��ĺ�����׼��
    if ti==500000 %��¼��ѹ��ʱ�䲽
        record_flag=1;
    end
    if record_flag>0 %��500000һֱȡ����ֹ������ʱ��ȡƽ���Խ���
        u_record(:,record_flag)=u;
        record_flag=record_flag+1;
    end
    
    %��ֹ�оݣ��ﵽ��ʱ����ʱ��ֹ
end
%----��ѭ������-----------------------------------------------------------------

now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp([now_str ' ��ѭ������'])
diary off
now_str=datestr(now,'HH_MM');

% ��ֹʱʵʱ���ͼ�񣨹������ݵ�ʾ�⣩�洢
figure(h_fig1)
saveas(gcf,[path_list.output '/��ֹʱʵʱ���' now_str '.png'])

% ��ֹ���˹��ж��ȶ���ʱ���ƿռ�ֲ�
figure
x_final=(0:1:200)*simulation.dx;
u_final=mean(u_record,2);
figure
plot(x_final,u_final,'-b','LineWidth',3)%���������ѹƽ��ֵ
axis([0,simulation.Lx,-inf,inf])
%         title('���ƿռ�ֲ�', 'FontSize', 18);
xlabel('x [m]')
ylabel('\phi [V]');
hold on;
line([simulation.Lx-6.66e-4,simulation.Lx],[2.51,2.51],'linestyle','-.','color','r','LineWidth',3);
L1=legend('����','����');
set(L1,'location','south');
set(L1,'AutoUpdate','off');
line([simulation.Lx-6.66e-4,simulation.Lx-6.66e-4],[0,2.51],'linestyle','-.','color','r','LineWidth',3);
saveas(gcf,[path_list.output '/��ֹʱ���Ʒֲ�' now_str '.png'])

%��ֹʱȫ�����ݴ洢
save([path_list.output '/final_data' now_str '.mat'])

%% ����
% ��ȡ�洢��������