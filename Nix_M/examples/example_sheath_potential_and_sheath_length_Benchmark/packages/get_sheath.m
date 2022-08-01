function [ sheath_potential, sheath_length] = get_sheath( Te, n0 ,i_out)
% function [ u_final] = get_sheath_potential( Te )
%GET_SHEATH_POTENTIAL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%1D PIC

%% ��ʼ��
%clear
%close all;

%�ļ�·��
addpath('../../packages')

% ����˵������λΪ���ʵ�λ��
% ȫ�ֱ���
constants=get_constants();% ȫ�ֳ��� �ṹ��

%--------�������-----------------------------------------------------------------------------------
simulation=get_simulation('default'); %������� �ṹ��
% ��ͨ��ȡ�����´���ע�����޸ķ������

% simulation.Lx=0.01;%�������򳤶�
% simulation.num_grid_point=201;%��������Ŀ
% simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
simulation.Te=Te;%e��λeV
simulation.n0=n0;%���������ܶ�
% debye=sqrt(constants.eps0*simulation.Te/(constants.e*simulation.n0));
simulation.Lx=0.002;%�������򳤶�
simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
simulation.source_region=[0.74*simulation.Lx,0.76*simulation.Lx];
% simulation.n0=1E16;%���������ܶ�

simulation.dt=0.1*2*pi/(56.41*sqrt(simulation.n0));%ʱ�䲽��
simulation.end_time=2500000*simulation.dt;%����ʱ���ܳ�
simulation.all_timesteps=floor(simulation.end_time/simulation.dt);%��ѭ������
simulation.num0_macro_e=40000;%��ʼ������Ŀ Particle e Num
simulation.num0_macro_H=simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.weight=simulation.n0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
simulation.TH=0.8;%��λeV
% simulation.field_boundaries_type=0;%��λ�߽���������
% simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
veth=sqrt(constants.e*simulation.Te/constants.me);%�����¶ȶ�Ӧ�����ٶ�
vHth=sqrt(constants.e*simulation.TH/constants.mH);%�����¶ȶ�Ӧ�����ٶ�
% numH1npertime=5;%ÿ��ʱ�䲽�������ĸ���������
% v0H1n=10000;%�������ӷ����ٶ�
simulation.pressure=0.6;%��ѹ ��λPa
        simulation.Tgas=300;%����   K
n0H=5e19;%������ԭ���ܶ�
T0H=0.8;%������ԭ���¶ȣ���λeV
e0H1n=1;%�⸺���ӷ�������

vineV=@(v,q_m_ratio)v.*v*3/(2*abs(q_m_ratio));
eVinv=@(eV,q_m_ratio)sqrt(eV*2*abs(q_m_ratio)/3);
gamaH=0.25*eVinv(T0H,constants.q_m_ration_H)*n0H;
%--------�������-----------------------------------------------------------------------------------

%--------MCC����-----------------------------------------------------------------------------------
%     [mcc_para,Xsec]=MCC_init(simulation,constants,pressure);
   [mcc_para,~,~,~,Xsec_Hp_H]= MCC_init(simulation,constants,simulation.pressure, '2020LZS');

%--------��ʼ��-------------------------------------------------------------------------------------
% TODO�����ӽṹ�壬�볡�ṹ��
%--------��������
%����Maxwell�ֲ����ɺ����ӵ��ٶȷֲ� x y z����(Particle e velocity)
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
vH1n=[zeros(0,1),zeros(0,1),zeros(0,1)];
%%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������

%���ȷֲ��ڿռ�Lx�ڣ�����һ������Ŷ�
position_e=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_e,1] );
position_H=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_H,1] );
% position_H1n=simulation.Lx*0*(ones(numH1npertime,1)+rand(numH1npertime,1)*0.01);
 position_H1n=zeros(0,1);

ae=zeros(simulation.num0_macro_e,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
aH=zeros(simulation.num0_macro_H,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
aH1n=zeros(0,3);%���ٶȳ�ֵ ���Ӽ��ٶ�
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

% MOV=moviein(floor(simulation.all_timesteps/dptime));
% count_mov=0;
% video = VideoWriter(['sheath_vB',num2str(i_out),'.avi']);
% open(video);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vH_acc_flag=zeros(size(vH,1),1);
vH1n_acc_flag=zeros(size(vH1n,1),1);

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
    if(~isempty(position_e))
    ve(:,1)=ve(:,1)+ae(:,1)*simulation.dt;
    position_e=position_e+ve(:,1)*simulation.dt;
    end
    
    vH(:,1)=vH(:,1)+aH(:,1)*simulation.dt;%��������������λ��
    position_H=position_H+vH(:,1)*simulation.dt;
    
    vH1n(:,1)=vH1n(:,1)+aH1n(:,1)*simulation.dt;
    position_H1n=position_H1n+vH1n(:,1)*simulation.dt;
    %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
    
    %------����Խ�紦��------------------------------------------------------------------------------------------------
    % TODO�����޸ģ���Ϊר�ú���
    %     % -------���ڱ߽�------------------------------------------------------
    %     tempFlag=(position_e<0);
    %     position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
    %     tempFlag=(position_e>simulation.Lx);
    %     position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
    
    % -------�ɶ���ע��Դ�� pair re-injection into source region--------------
    if(~isempty(position_e))
    tempFlag=(position_e<0|position_e>simulation.Lx);
    ve(tempFlag,:)=[];%Ĩ����������ĵ���
    position_e(tempFlag)=[];
    end
    
    k1=(position_H>simulation.Lx);%ȡ������Χ�����ӵı��
    tempsum=sum(k1);
    
    if(tempsum>0)%����ǿ�
        vH(k1,:)=normrnd(0,vHth,[tempsum,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
        % ��ע������λ�þ�������ֲ���Դ����
        position_H(k1)=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
        temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
        position_e=reshape([position_e; temp_matrix],[],1);%�ϲ���ԭ��������
        ve=reshape([ve; normrnd(0,veth,[tempsum,3])],[],3);
    end
    
    k2=(position_H<0);%ȡ������Χ�����ӵı��
    tempsum=sum(k2);
    
    if(tempsum>0)%����ǿ�
        vH(k2,:)=normrnd(0,vHth,[tempsum,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
        % ��ע������λ�þ�������ֲ���Դ����
        position_H(k2)=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
        temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum,1);
        position_e=reshape([position_e; temp_matrix],[],1);%�ϲ���ԭ��������
        ve=reshape([ve; normrnd(0,veth,[tempsum,3])],[],3);
    end
    
    tempsum=sum(k2);%���¼����������Ӳ������⸺����
    positionH1nnew=zeros(tempsum,1);
    vH1nnew=zeros(tempsum,3);
    p_reflect=zeros(tempsum,1);
%     p_reflect=ones(tempsum,1)*0.7;
    j=1;
    for i=1:size(k2,1)
        if(j>tempsum)
            break;
        end
        if(k2(i,1)==1)
            vH1nnew(j,1)=eVinv(e0H1n,constants.q_m_ration_H1n);
            vH1nnew(j,2)=0;
            vH1nnew(j,3)=0;
            p_reflect(j,1)=0.3*(1-2/vineV(vH(i,1),constants.q_m_ration_H));
            j=j+1;            
        end
    end
    delete_flag=0;
    if(size(positionH1nnew,1)>1)
        for i=1:size(positionH1nnew,1)
            if rand(1)>p_reflect(i-delete_flag,1)
                positionH1nnew(i-delete_flag,:)=[];
                vH1nnew(i-delete_flag,:)=[];
                p_reflect(i-delete_flag,:)=[];
                delete_flag=delete_flag+1;
            end
        end
    end
    position_H1n=reshape([position_H1n; positionH1nnew],[],1);
    vH1n=reshape([vH1n; vH1nnew],[],3);
    
    %���¼�����ԭ�Ӳ������⸺����
    HtoH1n=round(gamaH*simulation.dt/simulation.weight);
    positionH1nnew2=zeros(HtoH1n,1);
    vH1nnew2=zeros(HtoH1n,3);
    vH1nnew2(:,1)=ones(HtoH1n,1)*eVinv(e0H1n,constants.q_m_ration_H1n);
    p_reflect=ones(HtoH1n,1)*0.42*exp(-1.05/T0H);
    delete_flag=0;
    if(size(positionH1nnew2,1)>1)
        for i=1:size(positionH1nnew2,1)
            if rand(1)>p_reflect(i-delete_flag,1)
                positionH1nnew2(i-delete_flag,:)=[];
                vH1nnew2(i-delete_flag,:)=[];
                p_reflect(i-delete_flag,:)=[];
                delete_flag=delete_flag+1;
            end
        end
    end
    position_H1n=reshape([position_H1n; positionH1nnew2],[],1);
    vH1n=reshape([vH1n; vH1nnew2],[],3);
    
    tempFlag2=(position_H1n<0|position_H1n>simulation.Lx);
    vH1n(tempFlag2,:)=[];%Ĩ����������ĸ�������
    position_H1n(tempFlag2)=[];
%     delete_flag=0;
%     for i=1:size(position_H1n,1)
%         if rand(1)<2e-5
%             position_H1n(i-delete_flag,:)=[];
%             vH1n(i-delete_flag,:)=[];
%             delete_flag=delete_flag+1;
%         end
%     end
    
%     while(length(position_H1n)>length(position_H))
%         delete=floor(rand()*min([length(position_H1n),length(vH1n),length(aH1n)]));
%         if(delete>0&&delete<length(position_H))
%             position_H1n(delete,:)=[];
%             vH1n(delete,:)=[];
%             aH1n(delete,:)=[];
%         end
%     end
    
    %--------thermalization �Ȼ�  ÿ5��һ��
    % ������������ڱ߽�ʹ��
    if(~isempty(position_e))
    if rem(ti,5)==1
        flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
        ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
    end
    end
    %------����Խ�紦�����------------------------------------------------------------------------------------------
    
    
     %------����MCC����------------------------------------------------------------------------------------------
        [collision_index,collision_type]=null_collision( vH,constants,mcc_para,Xsec_Hp_H,'H' );
        [colli_result]=collision_process(Xsec_Hp_H, mcc_para,vH,constants,collision_index,collision_type ,'H');
        v_cold=sqrt(constants.e*0.4/constants.mH);
        if colli_result.happen==1
            vH(colli_result.ion_index,:)=get_v_init( v_cold, 'Maxwellian velocity', [length(colli_result.ion_index),3]);%*eVinv(T0H,constants.q_m_ration_H);
            vH(colli_result.ela_index,:)=get_v_init( v_cold, 'Maxwellian velocity', [length(colli_result.ela_index),3]);%*eVinv(T0H,constants.q_m_ration_H);
            vH(colli_result.exc_index,:)=get_v_init( v_cold, 'Maxwellian velocity', [length(colli_result.exc_index),3]);%*eVinv(T0H,constants.q_m_ration_H);
%             vH(colli_result.att_index,:)=[];
%             position_H(colli_result.att_index,:)=[];
            vH=[vH; colli_result.ion_gen_v ];
            position_H=[position_H;position_H(colli_result.ion_index)];
            %         vH2=[vH2; normrnd(0,vHth,[length(colli_result.ion_gen_v) , 3]);  ];
        end
        
        %------����MCC�������------------------------------------------------------------------------------------------
        
    
    
    %-----------������----------------------------------------------------------------------------------
    rhoea=zeros(simulation.num_grid_point,1);%��������������������ܶ� ����
    rhoeb=zeros(simulation.num_grid_point,1);
    rhoaH=zeros(simulation.num_grid_point,1);%����
    rhobH=zeros(simulation.num_grid_point,1);
    rhoaH1n=zeros(simulation.num_grid_point,1);%��������
    rhobH1n=zeros(simulation.num_grid_point,1);
    
    if(~isempty(position_e))
    enearx=floor(position_e/simulation.dx)+1;%ÿ������������������� �������
    eassigndx=position_e-(enearx-1)*simulation.dx;%ÿ�����Ӿ���������������ľ���
    Eenearx=floor((position_e+0.5*simulation.dx)/simulation.dx)+1;%ÿ���������ڰ���������� ����
    Eedx=position_e-(Eenearx*simulation.dx-3/2*simulation.dx);%ÿ�����Ӿ����������İ�����ľ���
    end
    
    Hnearx=floor(position_H/simulation.dx)+1;
    Hassigndx=position_H-(Hnearx-1)*simulation.dx;
    EHnearx=floor((position_H+0.5*simulation.dx)/simulation.dx)+1;%ÿ���������ڰ���������� ����
    EHdx=position_H-(EHnearx*simulation.dx-3/2*simulation.dx);
    
    H1nnearx=floor(position_H1n/simulation.dx)+1;
    H1nassigndx=position_H1n-(H1nnearx-1)*simulation.dx;
    EH1nnearx=floor((position_H1n+0.5*simulation.dx)/simulation.dx)+1;%ÿ���������ڰ���������� ����
    EH1ndx=position_H1n-(EH1nnearx*simulation.dx-3/2*simulation.dx);
    
   if(~isempty(position_e))
    for j=1:size(position_e)
        rhoea(enearx(j))=rhoea(enearx(j))+(simulation.dx-eassigndx(j));
        rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)+(eassigndx(j));
        %         rhoeb(enearx(j)+1)=rhoeb(enearx(j)+1)-e/(dx^2)*(assigndx(j))*weight;
    end
    end
    rhoe=(rhoea+rhoeb)*(-constants.e)/(simulation.dx^2)*simulation.weight;
    %rhoP=rhoP*0;
    for j=1:size(position_H)
        rhoaH(Hnearx(j))=rhoaH(Hnearx(j))+(simulation.dx-Hassigndx(j));
        rhobH(Hnearx(j)+1)=rhobH(Hnearx(j)+1)+(Hassigndx(j));
        %         rhobH(Hnearx(j)+1)=rhobH(Hnearx(j)+1)+e/(dx^2)*(Hassigndx(j))*weight;
    end
    rhoH=(rhoaH+rhobH)*constants.e/(simulation.dx^2)*simulation.weight;
    
    for j=1:size(position_H1n)
        rhoaH1n(H1nnearx(j))=rhoaH1n(H1nnearx(j))+(simulation.dx-H1nassigndx(j));
        rhobH1n(H1nnearx(j)+1)=rhobH1n(H1nnearx(j)+1)+(H1nassigndx(j));
        %         rhobH(Hnearx(j)+1)=rhobH(Hnearx(j)+1)+e/(dx^2)*(Hassigndx(j))*weight;
    end
    rhoH1n=(rhoaH1n+rhobH1n)*-constants.e/(simulation.dx^2)*simulation.weight;
    
    rho=rhoe+rhoH+rhoH1n;
    %-----------�����ɽ���------------------------------------------------------------
    
    vH_grid=zeros(simulation.num_grid_point,1);
    count_vH=zeros(simulation.num_grid_point,1);
    
    for j=1:size(position_H)
        if vH_acc_flag(j,1)>=1
            vH_acc_flag(j,1)=1+vH_acc_flag(j,1);
        end
        if k1(j)==1
            vH_acc_flag(j,1)=1;
        end
        if vH_acc_flag(j,1)==1
            vH_acc_flag(j,1)=0;
        end
        if vH_acc_flag(j,1)==0
        vH_grid(Hnearx(j))=vH_grid(Hnearx(j))+abs(vH(j,1));
        count_vH(Hnearx(j))=count_vH(Hnearx(j))+1;
        end
    end
    vH_grid=vH_grid./count_vH;
    
    vH1n_grid=zeros(simulation.num_grid_point,1);
    count_vH1n=zeros(simulation.num_grid_point,1);
    
    vH1n_acc_flag=zeros(size(vH1n,1),1);
    for j=1:size(position_H1n)
        if vH1n_acc_flag(j,1)>=1
            vH1n_acc_flag(j,1)=1+vH1n_acc_flag(j,1);
        end
%         if k1(j)==1
%             vH1n_acc_flag(j,1)=1;
%         end
        if vH1n_acc_flag(j,1)==1
            vH1n_acc_flag(j,1)=0;
        end
        if vH1n_acc_flag(j,1)==0
        vH1n_grid(H1nnearx(j))=vH1n_grid(H1nnearx(j))+abs(vH1n(j,1));
        count_vH1n(H1nnearx(j))=count_vH1n(H1nnearx(j))+1;
        end
    end
    vH1n_grid=vH1n_grid./count_vH1n;
    
    
    
    
    %--------------���糡---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %--------------���糡����---------------------------------------------------------------
    
    % TODO�������E_particle
    ae=-((simulation.dx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/simulation.dx;
    aH=((simulation.dx-EHdx).*E(EHnearx)+EHdx.*E(EHnearx+1))*constants.e/(constants.mH)/simulation.dx;
    aH1n=-((simulation.dx-EH1ndx).*E(EH1nnearx)+EH1ndx.*E(EH1nnearx+1))*constants.e/(constants.mH1n)/simulation.dx;
    
    % ��ѹ�ۼӣ�Ϊ  ʵʱ������ȡƽ��  ��׼��
    if rem(ti+jj,dptime)==1
        usum=usum+u;
        jj=jj-1;
        if jj==-1
            jj=avgSteps;
        end
    end
    
%     count_mov=count_mov+1;
    
    
    %-----------ʵʱ���----------------------------------------------------------------------------------
    if rem(ti-1,dptime)==0%ÿ��dptime��ʾһ��
        fprintf('��ǰʱ������%d \n',ti)
        
        e_num(floor(ti/dptime)+1)=size(position_e,1);%�����ڵĵ��������仯
        H_num(floor(ti/dptime)+1)=size(position_H,1);
        H1n_num(floor(ti/dptime)+1)=size(position_H1n,1);
        e_vAve(floor(ti/dptime)+1)=mean(abs(ve(:,1)));%x�����ϵ�ƽ�����ʱ仯
        H_vAve(floor(ti/dptime)+1)=mean(abs(vH(:,1)));
        figure(3*i_out-2)
        subplot(3,2,1);
        plot_no_versus_position( position_e,position_H,simulation );
        subplot(3,2,2);
        plot_num_versus_timestep( e_num,H_num,H1n_num,dptime )
        subplot(3,2,3);
        plot_2u1E_versus_x( u, ti, usum, avgSteps, E, simulation )
        subplot(3,2,4);
%         plot_Ek_versus_timestep( e_vAve,H_vAve, constants.q_m_ration_e, constants.q_m_ration_H, 'e','H', dptime )
plot(vH1n_grid,'-b','LineWidth',3)%������������ٶ�ƽ��ֵ
        xlabel('������');
        ylabel('H-�ٶ�');
        axis([0,200,-inf,inf])
%         ylabel('ƽ��Ek [eV]'); %ʵ������3��x��ƽ������
        subplot(3,2,5);%��ɷֲ�
        plot_density_versus_x( rhoe, rhoH, rho, simulation )
        subplot(3,2,6);
%         plot_vH_acc_timestep(vH,k1,ti,vH_acc_flag)
%                 plot_v_versus_position( ve, vH, position_e,position_H,simulation )
%         plot_vineV_and_u( vH_grid, u, constants.q_m_ration_H, constants,Te )
        plot(rhoH1n/-constants.e,'-b','LineWidth',3)%�⸺���ӷֲ�
        hold on;
        plot(rhoH/constants.e,'-r','LineWidth',3)%�⸺���ӷֲ�
        hold on;
        plot(rhoe/-constants.e,'-k','LineWidth',3)%�⸺���ӷֲ�
        legend('H-','H+','e')
        xlabel('������');
        ylabel('�����ܶ�');
        axis([0,size(rhoH1n,1),-inf,inf])
        hold off
%         uB=sqrt(constants.e*Te/constants.mH);
%         line([0,size(vH_grid,1)],[uB,uB],'linestyle','-.','color','r','LineWidth',3)
        drawnow;
        
        usum=0; %����
        
        figure(3*i_out)        
%         plot(vH_grid,'-b','LineWidth',3)%������������ٶ�ƽ��ֵ
%         axis([0,size(vH_grid,1),-inf,inf])
%         uB=sqrt(constants.e*Te/constants.mH);
%         line([0,size(vH_grid,1)],[uB,uB],'linestyle','-.','color','r','LineWidth',3)
plot_vineV_and_u( vH_grid, u, constants.q_m_ration_H, constants,Te )        
% MOV(count_mov)=getframe;
%         writeVideo(video,MOV(count_mov) );
% H=1000;W=1000;
%         MOV=getframe;
%         MOV.cdata=imresize(MOV.cdata,[H W]);
%         writeVideo(video,MOV );
        
    end
    %-----------ʵʱ����----------------------------------------------------------------------------------
    
    if ti==50000 %�ȶ��о�
        recordFlag=1;
    end
    if recordFlag>0 %�ȶ���ȡ��ѹ������ʱ��ȡƽ���Խ���
        u_record(:,recordFlag)=u;
        recordFlag=recordFlag+1;
    end
    
%     if mod(ti,5000)==0
%         filename=strcat('����ʱ����',num2str(ti),'.mat');
%         save(filename);
%     end
    
end
%----��ѭ������-----------------------------------------------------------------
% movie(MOV,1);
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
diary(log_name)
disp(now_str)
diary off

% saveas(gcf,'./others/real_time_display.png')

x_final=(0:1:200)*simulation.dx;
u_final=mean(u_record,2);
figure(3*i_out-1)
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
% line([Llimit*simulation.dx,Llimit*simulation.dx],[0,max(u_final)],'linestyle','-.','color','r','LineWidth',3);
% faiw=-0.5*Te*log(2*pi*constants.me/constants.mH*(1+simulation.TH/simulation.Te));
% line([0,simulation.Lx],[faiw,faiw],'linestyle','-.','color','k','LineWidth',2);hold on;
% L1=legend('������','��ѹƽ��ȡֵ��','���ۼ����ʲ��λ');
L1=legend('������');
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

filename=strcat('��ԭ���ܶ�',num2str(n0H),'.mat');
        save(filename);



% vH_all=zeros(1,size(vH,1)/10);
% for i=1:size(vH,1)/10
%     for j=1:10
%     vH_all(1,i)=vH_all(1,i)+vH(j,1)*3/2;
%     end
% end
% vH_all=vH_all/10;

% figure(3*i_out)
% x_v=linspace(0,simulation.Lx,size(vH,2));
% plot(x_v,vH,'-b','LineWidth',3)%���������ѹƽ��ֵ
% hold on;
% axis([0,simulation.Lx,-inf,inf])
% uB=sqrt(constants.e*Te/constants.mH);
% line([0,simulation.Lx],[uB,uB],'linestyle','-.','color','r','LineWidth',3)

% saveas(gcf,'./others/final_display.png')

% save('./others/final_data.mat')

% TODO�����޸Ĺ�ʽ����eT�滻kT
% debye=sqrt(constants.eps0*kB/constants.e/constants.e/(1e16/Te+1e16/TH))     �°ݳ��ȼ���ʽ
%
% Vf=0.5*log((2*pi*constants.me/constants.mH)*(1+TH/Te))*constants.e*Te/constants.e    ���������������ֵ


end