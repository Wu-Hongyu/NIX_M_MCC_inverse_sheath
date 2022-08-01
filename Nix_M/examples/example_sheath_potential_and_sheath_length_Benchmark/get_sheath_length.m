function [ sheath_length] = get_sheath_length( n0 )
%GET_SHEATH_POTENTIAL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%1D PIC

%% ��ʼ��
%clear
%close all;

%�ļ�·��
addpath('../packages')

% ����˵������λΪ���ʵ�λ��
% ȫ�ֱ���
constants=get_constants();% ȫ�ֳ��� �ṹ��

%--------�������-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %������� �ṹ��
% ��ͨ��ȡ�����´���ע�����޸ķ������
% simulation.dt=10*10^-12;%ʱ�䲽��
% simulation.end_time=200000*10*10^-12;%����ʱ���ܳ�
% simulation.all_timesteps=floor(simulation.end_time/simulation.dt);%��ѭ������
% simulation.Lx=0.01;%�������򳤶�
simulation.num_grid_point=401;%��������Ŀ
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
simulation.n0=n0;%���������ܶ�
debye=sqrt(constants.eps0*simulation.Te/(constants.e*simulation.n0));
simulation.dx=debye;
simulation.Lx=debye*(simulation.num_grid_point-1);%�������򳤶�
simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];

simulation.dt=0.1*2*pi/(56.41*sqrt(simulation.n0));%ʱ�䲽��
simulation.end_time=100000*simulation.dt;%����ʱ���ܳ�
simulation.all_timesteps=floor(simulation.end_time/simulation.dt);%��ѭ������

simulation.num0_macro_e=40000;%��ʼ������Ŀ Particle e Num
simulation.num0_macro_H=simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.weight=simulation.n0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
% simulation.TH=1;%��λeV
% simulation.field_boundaries_type=0;%��λ�߽���������
% simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
veth=sqrt(constants.e*simulation.Te/constants.me);%�����¶ȶ�Ӧ�����ٶ�
vHth=sqrt(constants.e*simulation.TH/constants.mH);%�����¶ȶ�Ӧ�����ٶ�
%--------�������-----------------------------------------------------------------------------------

%--------��ʼ��-------------------------------------------------------------------------------------
% TODO�����ӽṹ�壬�볡�ṹ��
%--------��������
%����Maxwell�ֲ����ɺ����ӵ��ٶȷֲ� x y z����(Particle e velocity)
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
%%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������

%���ȷֲ��ڿռ�Lx�ڣ�����һ������Ŷ�
position_e=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_e,1] );
position_H=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_H,1] );

ae=zeros(simulation.num0_macro_e,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
aH=zeros(simulation.num0_macro_H,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
%--------������
[ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
%
%--------��ʼ��--------------------------------------------------------------------------------------

% ����������ز���
dptime=1000;%ÿ���ٲ���������ʾ
avgSteps=499;%���ٲ����ڵ�ƽ��ֵ
jj=499;
usum=0;
recordFlag=0;
u_record=zeros(simulation.num_grid_point,1);%��ʱ����

log_name = get_log_init( simulation, '' );%�����ʼ��־
h = figure('Unit','Normalized','position',...
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
    
    vH(:,1)=vH(:,1)+aH(:,1)*simulation.dt;%��������������λ��
    position_H=position_H+vH(:,1)*simulation.dt;
    %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
    
    %------����Խ�紦��------------------------------------------------------------------------------------------------
    % TODO�����޸ģ���Ϊר�ú���
    %     % -------���ڱ߽�------------------------------------------------------
    %     tempFlag=(position_e<0);
    %     position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
    %     tempFlag=(position_e>simulation.Lx);
    %     position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
    
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
    
    %--------thermalization �Ȼ�  ÿ5��һ��
    % ������������ڱ߽�ʹ��
    if rem(ti,5)==1
        flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
        ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
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
        fprintf('��ǰʱ������%d \n',ti)
        
        e_num(floor(ti/dptime)+1)=size(position_e,1);%�����ڵĵ��������仯
        H_num(floor(ti/dptime)+1)=size(position_H,1);
        e_vAve(floor(ti/dptime)+1)=mean(abs(ve(:,1)));%x�����ϵ�ƽ�����ʱ仯
        H_vAve(floor(ti/dptime)+1)=mean(abs(vH(:,1)));
        
        subplot(3,2,1);
        plot_no_versus_position( position_e,position_H,simulation );
        subplot(3,2,2);
        plot_num_versus_timestep( e_num,H_num,dptime )
        subplot(3,2,3);
        plot_2u1E_versus_x( u, ti, usum, avgSteps, E, simulation )
        subplot(3,2,4);
        plot_Ek_versus_timestep( e_vAve,H_vAve, constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', dptime )
        ylabel('ƽ��Ek [eV]'); %ʵ������3��x��ƽ������
        subplot(3,2,5);%��ɷֲ�
        plot_density_versus_x( rhoe, rhoH, rho, simulation )
        subplot(3,2,6);
        %         plot_v_versus_position( ve, vH, position_e,position_H,simulation )
        plot_vineV_versus_position( ve, vH, position_e,position_H,constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', simulation )
        drawnow;
        
        usum=0; %����
    end
    %-----------ʵʱ����----------------------------------------------------------------------------------
    
    if ti==70000 %�ȶ��о�
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

% saveas(gcf,'./others/real_time_display.png')

x_final=(0:1:(simulation.num_grid_point-1))*simulation.dx;
u_final=mean(u_record,2);
figure
plot(x_final,u_final,'-b','LineWidth',3)%���������ѹƽ��ֵ
axis([0,simulation.Lx,-inf,inf])
%         title('���ƿռ�ֲ�', 'FontSize', 18);

%ȡ�м��3/5��4/5���ȵĵ�ѹ��ƽ������Ϊ�ʲ��ѹ��ֵ��
Llimit=floor(simulation.num_grid_point*3/5);
Hlimit=floor(simulation.num_grid_point*4/5);
sheath_potential=mean(u_final(Llimit:Hlimit));

xlabel('x [m]')
ylabel('\phi [V]');
hold on;
% line([simulation.Lx-6.66e-4,simulation.Lx],[2.51,2.51],'linestyle','-.','color','r','LineWidth',3);
line([Llimit*simulation.dx,Llimit*simulation.dx],[0,max(u_final)],'linestyle','-.','color','r','LineWidth',3);
L1=legend('������','��ѹƽ��ȡֵ��');
set(L1,'location','south');
set(L1,'AutoUpdate','off');
% line([simulation.Lx-6.66e-4,simulation.Lx-6.66e-4],[0,2.51],'linestyle','-.','color','r','LineWidth',3);
line([Hlimit*simulation.dx,Hlimit*simulation.dx],[0,max(u_final)],'linestyle','-.','color','r','LineWidth',3);

for i=Hlimit:simulation.num_grid_point
    if u_final(i)<0.95*sheath_potential
        sheath_length=simulation.Lx-(i-1)*simulation.dx;
        break
    end
end
% saveas(gcf,'./others/final_display.png')

% save('./others/final_data.mat')

% TODO�����޸Ĺ�ʽ����eT�滻kT
% debye=sqrt(constants.eps0*kB/constants.e/constants.e/(1e16/Te+1e16/TH))     �°ݳ��ȼ���ʽ
%
% Vf=0.5*log((2*pi*constants.me/constants.mH)*(1+TH/Te))*constants.e*Te/constants.e    ���������������ֵ


end

