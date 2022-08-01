%1D PIC
% ˫�����ȶ���

%% ��ʼ��
clear
close all;

if exist('get_path_init','file')~=2
    addpath('../../packages') % ���./packages��·��
end
path_list=get_path_init('two_stream'); %�ļ�·��

% ����˵������λΪ���ʵ�λ��
% ȫ�ֱ���
constants=get_constants();% ȫ�ֳ��� �ṹ��

%--------�������-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %������� �ṹ��
% ��ͨ��ȡ�����´���ע�����޸ķ��������ע��һЩ���������ڹ�����
% simulation.dt=10*10^-12;%ʱ�䲽��
% simulation.end_time=6E-6;%����ʱ���ܳ�
simulation.all_timesteps=1500;%��ѭ������
% simulation.Lx=0.01;%�������򳤶�
% simulation.source_region=[0.1*simulation.Lx,0.4*simulation.Lx];
% simulation.num_grid_point=201;%��������Ŀ
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
% simulation.ne0=1E16;%���������ܶ�
% simulation.num0_macro_e=20000;%��ʼ������Ŀ Particle e Num
% simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
% simulation.Te0=1;%e��λeV
simulation.field_boundaries_type=3;%��λ�߽���������
% simulation.field_boundaries=[0,0];%��λ�߽�����ֵ

simulation.num0_macro_e=20000;%��ʼ������Ŀ Particle e Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��

veth=sqrt(constants.e*simulation.Te0/constants.me);%�����¶ȶ�Ӧ�����ٶ�
% vb=4*veth; %�����ٶ�
vb=0;
num_b=simulation.num0_macro_e/2;
rho_positive_background=simulation.ne0*constants.e; % ����������������Ϊ����
%--------�������-----------------------------------------------------------------------------------

%--------��ʼ��-------------------------------------------------------------------------------------
% TODO�����ӽṹ�壬�볡�ṹ��
%--------��������
% �ȵ�����
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,1]);
ve(1:num_b)=vb+ve(1:num_b);
ve((num_b+1):end)=-vb+ve((num_b+1):end);
% % �������
% ve=ones(simulation.num0_macro_e,1);
% ve(1:num_b)=vb*ve(1:num_b);
% ve((num_b+1):end)=-vb*ve((num_b+1):end);
%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������
Ek=@(v,q_m_ratio) v.*v*3/(2*abs(q_m_ratio)); %���ܣ���λeV
Eke0=simulation.weight*sum(Ek(ve,constants.q_m_ration_e)); %��ʼ1D�ܶ���

% ��������ֲ����Դ��Ŷ�
position_e=get_position_init(simulation, 'entire domain uniform random', [simulation.num0_macro_e,1] );

ae=zeros(simulation.num0_macro_e,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
%--------������
[ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
%
%--------��ʼ��--------------------------------------------------------------------------------------

% ������ز���
dptime=10;%ÿ���ٲ���������ʾ
avgSteps=499;%���ٲ����ڵ�ƽ��ֵ
jj=499;
usum=0;
recordFlag=0;
u_record=zeros(simulation.num_grid_point,1);%��ʱ����

log_name = get_log_init( path_list.output,'' );%�����ʼ��־
h_fig1 = figure('Unit','Normalized','position',...
    [0.02 0.3 0.6 0.6]); %ʵʱ���ʹ�õĴ󴰿�
% �������
% h_fig2=figure;
% video = VideoWriter('two_stream_instability.avi');
% open(video);
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
    %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
    
    %------����Խ�紦��------------------------------------------------------------------------------------------------
    % -------���ڱ߽�------------------------------------------------------
    tempFlag=(position_e<0);
    position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
    tempFlag=(position_e>simulation.Lx);
    position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
    %     �ο� л���� pic_2strm.m �������ڱ߽�����
    %     xp=xp./L+10.0;
    %     xp=L.*(xp-floor(xp));
    % �ر��������������ڱ߽磬��Ч��
    %     if xp>0
    %         xp=mod(xp,L) ȡ��
    %     elseif -10*L<xp<0
    %         xp=L-mod(|xp|,L)
    %     end
    %------����Խ�紦�����------------------------------------------------------------------------------------------
    
    
    %-----------������----------------------------------------------------------------------------------
    rhoea=zeros(simulation.num_grid_point,1);%��������������������ܶ� ����
    rhoeb=zeros(simulation.num_grid_point,1);
    
    enearx=floor(position_e/simulation.dx)+1;%ÿ������������������� �������
    eassigndx=position_e-(enearx-1)*simulation.dx;%ÿ�����Ӿ���������������ľ���
    Eenearx=floor((position_e+0.5*simulation.dx)/simulation.dx)+1;%ÿ���������ڰ���������� ����
    Eedx=position_e-(Eenearx*simulation.dx-3/2*simulation.dx);%ÿ�����Ӿ����������İ�����ľ���
    
    for j=1:size(position_e)
        rhoea(enearx(j))=rhoea(enearx(j))+(simulation.dx-eassigndx(j));
        rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)+(eassigndx(j));
        %         rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)-e/(dx^2)*(assigndx(j))*weight;
    end
    rhoe=(rhoea+rhoeb)*(-constants.e)/(simulation.dx^2)*simulation.weight;
    %rhoP=rhoP*0;
    rho=rhoe+rho_positive_background;
    %-----------�����ɽ���------------------------------------------------------------
    
    %--------------���糡---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------���糡����---------------------------------------------------------------
    
    % TODO�������E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    
    % ��ѹ�ۼӣ�Ϊ  ʵʱ������ȡƽ��  ��׼��
    if rem(ti+jj,dptime)==1
        usum=usum+u;
        jj=jj-1;
        if jj==-1
            jj=avgSteps;
        end
    end
    
    
    %-----------ʵʱ���----------------------------------------------------------------------------------
    if rem(ti-1,dptime)==0%ÿ��dptime��ʾһ��
        display_i=floor(ti/dptime)+1;
        position_e1=position_e(1:num_b,1);
        position_e2=position_e((num_b+1):end,1);
        ve1=ve(1:num_b,1);
        ve2=ve((num_b+1):end,1);
        
        fprintf('��ǰʱ������%d \n',ti)
        
        e1_num(display_i)=size(position_e1,1);%�����ڵĵ��������仯
        e2_num(display_i)=size(position_e2,1);
        e1_vAve(display_i)=mean(abs(ve1(:,1)));%x�����ϵ�ƽ�����ʱ仯
        e2_vAve(display_i)=mean(abs(ve2(:,1)));
        Eke(display_i)=simulation.weight*sum(Ek(ve,constants.q_m_ration_e))/Eke0; %��һ��1D�ܶ���
        Ep(display_i)=0.5*constants.eps0*simulation.dx*sum(E(2:end-1).*E(2:end-1))/constants.e/Eke0; %��һ��1D�ܾ����ܣ���λeV
        
        % ���ӵ�����
        Phi_e=((simulation.dx-eassigndx).*u(enearx)+eassigndx.*u(enearx+1))/simulation.dx;
        Epe(display_i)=simulation.weight*sum(abs(Phi_e))/Eke0;
        
        figure(h_fig1)
        subplot(3,2,1);
        %         plot_no_versus_position( position_e1,position_e2,simulation );
        % ��ͼ��ʾ ���ӱ��-����λ��
        temp_vec=1:simulation.num0_macro_e/2;
        scatter(position_e1,temp_vec,1,'b')
        axis([0 simulation.Lx 0 simulation.num0_macro_e]);
        hold on
        scatter(position_e1,simulation.num0_macro_e/2+temp_vec,1,'r')
        %         title('���ӿռ�ֲ�');
        xlabel('x [m]');
        ylabel('���ӱ��');
        legend('+x�� e','-x�� e')
        hold off
        subplot(3,2,2);
        plot_num_versus_timestep( e1_num,e2_num,dptime )
        subplot(3,2,3);
        plot_1u1E_versus_x( u, ti, E, simulation )
        subplot(3,2,4);
        plot_Ek_versus_timestep( e1_vAve,e2_vAve, constants.q_m_ration_e, constants.q_m_ration_e, '+x�� e','-x�� e',dptime )
        ylabel('mean Ek [eV]'); %ʵ������3��x��ƽ������
        subplot(3,2,5);
        % �����غ�-л����  2020/11/05 17:04:51 ʧ�ܣ�Ϊʲô����ͬһ������
        % ��ͼ��ʾ �����ܶ�����糡�ܴ���-ʱ�䲽
        plot(Eke,'-b')
        hold on
        plot(Epe,'--r')
        plot(Eke+Epe,'-.k')
        xlabel(['t [' num2str(dptime) 'ʱ��]']);
        ylabel('��һ������');
%         axis([0,inf,0,3])
        L1=legend('Ek','Ep','Ek+Ep');
        set(L1,'location','east');
        hold off
        subplot(3,2,6);
        plot_vineV_versus_position( ve1, ve2, position_e1,position_e2,constants.q_m_ration_e, constants.q_m_ration_e, '+x�� e','-x�� e', simulation )
        drawnow;
        
%         % �������
%         figure(h_fig2) %ʹ��ΨһID
%         plot_vineV_versus_position( ve1, ve2, position_e1,position_e2,constants.q_m_ration_e, constants.q_m_ration_e, '+x�� e','-x�� e', simulation )
%         M(display_i) = getframe(h_fig2);
%         writeVideo(video,M(display_i) );
%         hold off
        
        usum=0; %����
    end
    %-----------ʵʱ����----------------------------------------------------------------------------------
    
    if ti==500000 %�ȶ��о�
        recordFlag=1;
    end
    if recordFlag>0 %�ȶ���ȡ��ѹ������ʱ��ȡƽ���Խ���
        u_record(:,recordFlag)=u;
        recordFlag=recordFlag+1;
    end
    
end
%----��ѭ������-----------------------------------------------------------------

now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp(now_str)
diary off

figure(h_fig1)
saveas(gcf,[path_list.output '/real_time_display.png'])

% ��һ��������-ʱ�䲽
inverse_plasma_frequency=2*pi/sqrt(simulation.ne0*constants.e^2/(constants.eps0*constants.me));
plot_t=(dptime*simulation.dt/inverse_plasma_frequency)*(1:display_i);
figure
Ep2=log(sqrt(Ep));
plot(plot_t,Ep2,'-r')
xlabel('t [2\pi\omega_{pe}^{-1}]');
ylabel('ln(E/Emax)/2');
hold off
saveas(gcf,[path_list.output '/������ʱ���ݻ�.png'])

% �������
close(video);
% figure(4)
% movie(M,1)

% save('./final_data.mat')
