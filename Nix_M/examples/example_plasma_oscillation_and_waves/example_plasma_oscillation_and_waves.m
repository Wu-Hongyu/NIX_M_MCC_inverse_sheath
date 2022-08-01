%1D PIC
% �����������벨
% ʹ�ã��޸�wave_type��Te
% wave_type��1*1 char������wave_type(1)
% 0. �����ϸ���ϵ�������������ģ�ͣ����ڳ�ʼʱ�ӽ�
% 1. �ų��µ��Ӳ���˹̹��

%% ��ʼ��
clear
close all;

if exist('get_path_init','file')~=2
    addpath('../../packages') % ���./packages��·��
end
path_list=get_path_init('plasma_waves'); %�ļ�·��
% test_all % ����ȫ������

% ����˵������λΪ���ʵ�λ��
% ȫ�ֱ���
constants=get_constants();% ȫ�ֳ��� �ṹ��

%--------�������-----------------------------------------------------------------------------------
wave_type='1'; 
ne0=1E16;
Te0=0; %�����
% Te0=1; %�ȵ���
end_time=6e-6;
Lx=0.01;
simulation=get_simulation('auto1',ne0, Te0, end_time, Lx);
% �������岨����
switch str2double(wave_type(1))
    case 0 % �������壨���磩��
        % ��t<T_pe/10, end_time>10*T_pe
        simulation.B0=0; %�޴ų�
        result=check_simulation_parameters( simulation, 1 );
        end_time=10*result.T_pe;
        simulation=get_simulation('auto0',ne0, Te0, end_time, Lx);
        simulation.dt=min([simulation.dt,result.T_pe/10]);
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);
    case 1 % �ų��µ��Ӳ���˹̹��
        % ��t<T_ce/10, end_time>10*T_ce and end_time>T_pe, Lx>2*r_ce
        simulation.B0=1; %���Ⱥ㶨�ų�
        result=check_simulation_parameters( simulation, 1 );
        end_time=max([12*result.T_ce, result.T_pe]);
        Lx=max([simulation.Lx, 2*result.r_ce]);
        simulation=get_simulation('auto0',ne0, Te0, end_time, Lx);
        simulation.dt=min([simulation.dt,result.T_ce/20]);
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);
        simulation.B0=1;
end
% �޽��������
simulation.field_boundaries_type=3;%��λ�߽���������
simulation.particle_boundaries_type=0;
%--------���Ӳ���-----------------------------------------------------------------------------------
% TODO�����ִ�����particle_group
% simulation.num0_macro_e=20000;%��ʼ������Ŀ Particle e Num
simulation.num0_macro_e=2000;
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
% ������ֹ�����ӱ���
rho_positive_background=ne0*constants.e;  % ������ܶ�
%--------ʵʱ����������ز���---------------------------------------------------------------------------
switch str2double(wave_type(1))
    case 0 % �������壨���磩��
        % ��t<T_pe/10, end_time>10*T_pe
        display_timesteps=max(floor([simulation.dt,result.T_pe/10]/simulation.dt));
    case 1 % �ų��µ��Ӳ���˹̹��
        % ��t<T_ce/10, end_time>10*T_ce and end_time>T_pe
%         display_timesteps=max(floor([simulation.dt,result.T_ce/10]/simulation.dt));
        display_timesteps=max(floor([simulation.dt,result.T_ce/20]/simulation.dt));
end
average_timesteps=499;%���ٲ����ڵ�ƽ��ֵ
node_sample=floor(simulation.num_grid_point/2);

%�����ʼ��־
log_name = get_log_init( path_list.output, '' );
% diary(log_name) % �ض���ʽ�����־
disp(simulation)
result=check_simulation_parameters( simulation, 2 );
disp('')
fprintf('display_timesteps=%d ;ÿ���ٲ���������ʾ \n',display_timesteps)
fprintf('average_timesteps= %d ;���ٲ����ڵ�ƽ��ֵ \n',average_timesteps)
diary off
pause %һ��ֻ����Ʋ���ʱʹ�ã��Բ鿴��������Ϣ

%--------��ʼ��-------------------------------------------------------------------------------------
% TODO�����ӽṹ�壬�볡�ṹ��
%--------��������
%����Maxwell�ֲ����ɺ����ӵ��ٶȷֲ� x y z����(Particle e velocity)
ve=get_v_init( result.vth_e, 'Maxwellian velocity', [simulation.num0_macro_e,3]);




%%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������
switch str2double(wave_type(1))
    case 0 % �������壨���磩�� 1D1V
        ve(:,2:3)=0;
        Ek=@(v,q_m_ratio) 3*sum(v.*v/(2*abs(q_m_ratio))); %���ܣ���λeV
    case 1 % �ų��µ��Ӳ���˹̹�� 1D3V
        ve=get_v_init( result.vth_e, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
        Ek=@(v,q_m_ratio) sum(v.*v/(2*abs(q_m_ratio))); %���ܣ���λeV
end
Eke0=simulation.weight*sum(Ek(ve,constants.q_m_ratio_e)); %��ʼ1D�ܶ���
if 0==Eke0
    Eke0=1;
end

% ��������ֲ����Դ��Ŷ�
position_e=get_position_init(simulation, 'entire domain uniform random', [simulation.num0_macro_e,1] );

ae=zeros(simulation.num0_macro_e,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
%--------������
[ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
%--------�������
j_average_counter=average_timesteps;
usum=0;
record_flag=0; %��ʼΪ0
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
    %������+�ų�ʸ���������µ���λ��
    ve(:,1)=ve(:,1)+(ae(:,1)-constants.q_m_ratio_e*simulation.B0*ve(:,3))*simulation.dt;
    ve(:,3)=ve(:,3)+constants.q_m_ratio_e*simulation.B0*simulation.dt*ve(:,1);% By
    position_e=position_e+ve(:,1)*simulation.dt;
    %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
    
    %------����Խ�紦��------------------------------------------------------------------------------------------------
    switch simulation.particle_boundaries_type
        case 0 % ���ڱ߽�
            % -------���ڱ߽�------------------------------------------------------
            tempFlag=(position_e<0);
            position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
            tempFlag=(position_e>simulation.Lx);
            position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
        otherwise
            throw(get_exception('error','no such type of particle boundaries.'))
    end
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
    rho=rhoe+rho_positive_background;
    %-----------�����ɽ���------------------------------------------------------------
    
    %--------------���糡---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------���糡����---------------------------------------------------------------
    
    % TODO�������E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    
    % ��ѹ�ۼӣ�Ϊ  ʵʱ�����ȡƽ��  ��׼��
    if rem(ti+j_average_counter,display_timesteps)==1
        usum=usum+u;
        j_average_counter=j_average_counter-1;
        if j_average_counter==-1
            j_average_counter=average_timesteps;
        end
    end
    
    %-----------ʵʱ���----------------------------------------------------------------------------------
    if rem(ti-1,display_timesteps)==0%ÿ��dptime��ʾһ��
        fprintf('��ǰʱ������%d/%d \n',ti,simulation.all_timesteps)
        
        display_i=ceil(ti/display_timesteps);
        e_num(display_i)=size(position_e,1);%�����ڵĵ��������仯
        e_vAve(display_i)=mean(abs(ve(:,1)));%x�����ϵ�ƽ�����ʱ仯
        Eke(display_i)=simulation.weight*sum(Ek(ve,constants.q_m_ratio_e))/Eke0; %��һ��1D�ܶ���
        Ep(display_i)=0.5*constants.eps0*simulation.dx*sum(E(2:end-1).*E(2:end-1))/constants.e/Eke0; %��һ��1D�ܾ����ܣ���λeV
        % ���ӵ�����
        Phi_e=((simulation.dx-eassigndx).*u(enearx)+eassigndx.*u(enearx+1))/simulation.dx;
        Epe(display_i)=simulation.weight*sum(abs(Phi_e))/Eke0;
        
        % ����ƽ������
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
        subplot(3,2,5);%��ɷֲ�
        % ������test_energy.m�޸�
%         % ��ͼ��ʾ �����ܶ�����糡�ܴ���-ʱ�䲽
%         plot(ave_Eke,'-b')
%         hold on
%         plot(ave_Epe,'--r')
%         plot(ave_Eke+ave_Epe,'-.k')
%         xlabel(['t [' num2str(display_timesteps) 'ʱ��]']);
%         ylabel('��һ������');
%         %         axis([0,inf,0,1.2])
%         L1=legend('Ek','Ep','Ek+Ep');
%         set(L1,'location','east');
%         hold off
        subplot(3,2,6);
        %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%         plot_vineV_versus_position( ve, 0*ve, position_e,0*position_e,constants.q_m_ratio_e, constants.q_m_ratio_e, 'e','', simulation )
%         ylim([-200,200])
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

%% Ƶ�׷���
% �����ź���Ϣ
dt_sampling=display_timesteps*simulation.dt; %����ʱ����
f_sampling=1/dt_sampling; %����Ƶ��
L_sampling=length(u_node_sample); %��������
% �����һ������������ʵ�ʳ���Ϊ׼�������жϺ�������
if mod(L_sampling,2) % ��Ϊż��
    L_sampling=L_sampling-1;
    E_node_sample=E_node_sample(1:L_sampling);
    u_node_sample=u_node_sample(1:L_sampling);
end

display_t=dt_sampling*(1:L_sampling); %����ʱ���

figure
subplot(2,1,1) % �����ź�
title('�����ź�')
switch str2double(wave_type(1))
    case 0 % �������壨���磩��
        % ��t<T_pe/10, end_time>10*T_pe
        plot_t=display_t/result.T_pe;
        plot(plot_t,E_node_sample,'-r')
        xlabel('t [2\pi\omega_{pe}^{-1}]');
        ylabel('E_{node sample}');
        for i=1:5
            line([i,i],[min(E_node_sample),max(E_node_sample)],'linestyle',':','color','k');
        end
        hold off
    case 1 % �ų��µ��Ӳ���˹̹��
        % ��t<T_ce/10, end_time>10*T_ce and end_time>T_pe, Lx>2*r_ce
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
Y_fft=fft(E_node_sample); %����Ƶ��
P2 = abs(Y_fft/L_sampling); %˫�߷���Ƶ��
f_plot2=f_sampling*(1:L_sampling)/L_sampling;

P1 = P2(1:L_sampling/2+1);
P1(2:end-1) = 2*P1(2:end-1); %���߷���Ƶ��
f_plot1 = f_sampling*(0:(L_sampling/2))/L_sampling;

subplot(2,1,2) %��ֵƵ��
switch str2double(wave_type(1))
    case 0 % �������壨���磩��
        % ��t<T_pe/10, end_time>10*T_pe
        plot_omega=2*pi*f_plot1/result.omega_pe;
        h_fig4=plot(plot_omega,P1,'b-'); %���߷���Ƶ��
        xlabel('\omega/\omega_{pe}')
        ylabel('|P1(f)|')
        for i=1:floor(max(plot_omega))
            line([i,i],[0,max(P1)],'linestyle',':','color','k');
        end
        
        
    case 1 % �ų��µ��Ӳ���˹̹��
        % ��t<T_ce/10, end_time>10*T_ce and end_time>T_pe, Lx>2*r_ce
        yyaxis left;
        % plot(f_plot2,P2) %˫�߷���Ƶ��
        % plot(f_plot1,P1) %���߷���Ƶ��
        % xlabel('f [Hz]')
        % ylabel('|P1(f)|')
        plot_omega=2*pi*f_plot1/result.omega_ce;
        h_fig4=plot(plot_omega,P1,'b-'); %���߷���Ƶ��
        xlabel('\omega/\omega_{ce}')
        ylabel('|P1(f)|')
        for i=1:floor(max(plot_omega))
            line([i,i],[0,max(P1)],'linestyle',':','color','k');
        end
        yyaxis right;
        h_fig5=plot(plot_omega,log10(P1),'--','Color',[0.8500    0.3250    0.0980]); %���߷���Ƶ��ȡ����
        ylabel('lg|P1(f)|')
        legend([h_fig4,h_fig5],'��ֵ','lg(��ֵ)')
        % % �˹����Ʒ�Χ����ͻ����Ϣ
        % yyaxis left;
        % ylim([0,7]) 
end
saveas(gcf,[path_list.output '/�ź���Ƶ��' now_str '.png'])

%��ֹʱȫ�����ݴ洢
save([path_list.output '/final_data' now_str '.mat'])