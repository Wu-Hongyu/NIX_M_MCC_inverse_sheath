%1D PIC
% ʹ�ã���������ģ�ͣ��������ø��������������Ӧ��ʼ����
% ������ѭ����ʹ����Ӧģ�飨������pusher�����ӱ߽�������

%% ��ʼ��
clear
close all;
HpH2ela=0;
HpH2exc=0;
HpHCX=0;
HpHela=0;

H2pH2CX=0;
H2pH2PX=0;
H2pHela=0;


%�ų�1000A��-34.022*x^5+13.888*x^4-2.0218*x^3+0.1308*x^2-0.0001*x+0.0006
if exist('get_path_init','file')~=2
    addpath('../../packages') % ���./packages��·��
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
simulation.ne0=2E16;%���������ܶ�
simulation.Te0=6;%��λeV
% %         TH=1;%��λeV
% % ʱ�ճ߶�
simulation.dt=50*10^-12;%ʱ�䲽��
simulation.end_time=1E-3;%����ʱ���ܳ�
simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
simulation.Lx=0.1;%�������򳤶�
simulation.source_region=[0.1*simulation.Lx,0.2*simulation.Lx];
simulation.num_grid_point=2001;%��������Ŀ
simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��

temp_time_rhoP=zeros(simulation.num_grid_point,1);
% % ��ų�
% simulation.field_boundaries_type=0;%��λ�߽���������
simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
% simulation.B0=0.005; %�����ԴŸ�Ӧǿ�� T
simulation.B0=0.00;
Bx_posi=0.05;%λ��
Bdelta=0.01;%��״����

Bfun=@(position_e) (exp(-(position_e-Bx_posi).^2/(2*Bdelta^2)));
% Bfun=@(x) (-34.022*x.^5+13.888*x.^4-2.0218*x.^3+0.1308*x.^2-0.0001*x+0.0006)*1000.*(x/simulation.Lx);
% Bfun=@(x) (-34.022*x.^5+13.888*x.^4-2.0218*x.^3+0.1308*x.^2-0.0001*x+0.0006)*1000;
% Bfun=@(position_e) 1;


% ����
simulation.particle_boundaries_type=1;%���ӱ߽��������� 0���� 1�ɶ�+�Ȼ�  2�ɶԺ�ע��+�Ȼ�
simulation.field_boundaries_type=0;%��λ�߽��������� 0����һ��߽�  3���ڱ߽�
inject_num_pertime=20;%ÿ��ע��2�����Ӻ�1��H+��2��H2+ 1��H3+  �����4�ı���
source_density=0;
source_left_bdidx=ceil((simulation.num_grid_point-1)*simulation.source_region(1)/simulation.Lx);
source_right_bdidx=floor((simulation.num_grid_point-1)*simulation.source_region(2)/simulation.Lx);
%��ѹ
simulation.pressure=0.6;%��ѹ ��λPa

ddx=simulation.dx;
ddt=simulation.dt;
%--------���Ӳ���-----------------------------------------------------------------------------------
% TODO�����ִ�����particle_group
simulation.THp=0.1;%��λeV
simulation.TH2p=0.1;%��λeV
simulation.TH3p=0.1;%��λeV
simulation.THn=0.1;%��λeV

Hp_part=0.25;
H2p_part=0.5;
H3p_part=0.25;

simulation.num0_macro_e=40000;%��ʼ������Ŀ Particle e Num
simulation.num0_macro_Hp=Hp_part*simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.num0_macro_H2p=H2p_part*simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.num0_macro_H3p=H3p_part*simulation.num0_macro_e;%��ʼ������ĿParticle H Num
simulation.weight=simulation.ne0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
simulation.weight_Hp=simulation.weight;
simulation.weight_H2p=simulation.weight;
simulation.weight_H3p=simulation.weight;
simulation.weight_Hn=simulation.weight;



% TODO: ���ڿ���ͨ��result=check_simulation_parameters����ֵ���result.vth_e
veth=sqrt(-constants.q_m_ratio_e*simulation.Te0);%�����¶ȶ�Ӧ�����ٶ�
vHpth=sqrt(constants.q_m_ratio_Hp*simulation.THp);%�����¶ȶ�Ӧ�����ٶ�
vH2pth=sqrt(constants.q_m_ratio_H2p*simulation.TH2p);%�����¶ȶ�Ӧ�����ٶ�
vH3pth=sqrt(constants.q_m_ratio_H3p*simulation.TH3p);%�����¶ȶ�Ӧ�����ٶ�
vHnth=sqrt(constants.q_m_ratio_Hp*simulation.THn);%�����¶ȶ�Ӧ�����ٶ�

%--------MCC����-----------------------------------------------------------------------------------

[mcc_para,...
    Xsec,Xsec_H,Xsec_Hp_H2,Xsec_Hp_H,Xsec_H2p_H2,Xsec_H2p_H,Xsec_H3p_H2,Xsec_H3p_H,Xsec_Hn_H2,Xsec_Hn_H]...
    =MCC_init(simulation,constants,simulation.pressure, '2020LZS');

%--------��ϲ���---------------------------------------------------------------------------
display_timesteps=5000;%ÿ���ٲ���������ʾ
average_timesteps=499;%���ٲ����ڵ�ƽ��ֵ

OB_position1=(simulation.source_region(1)+simulation.source_region(2))/2;
OB_position2=0.5*simulation.Lx;%EEPF�Ĺ۲��
OB_position3=0.8*simulation.Lx;
ve_static1=[];
ve_static2=[];
ve_static3=[];

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
% if simulation.particle_boundaries_type==2
%     simulation.num0_macro_e=0.01* simulation.num0_macro_e;
%     simulation.num0_macro_Hp=0.01*simulation.num0_macro_Hp;
%     simulation.num0_macro_H2p=0.01*simulation.num0_macro_H2p;
% end
%����Maxwell�ֲ����ɺ����ӵ��ٶȷֲ� x y z����(Particle e velocity)
ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
vHp=get_v_init( vHpth, 'Maxwellian velocity', [simulation.num0_macro_Hp,3]);
vH2p=get_v_init( vH2pth, 'Maxwellian velocity', [simulation.num0_macro_H2p,3]);
vH3p=get_v_init( vH3pth, 'Maxwellian velocity', [simulation.num0_macro_H3p,3]);
vHn=get_v_init( vHnth, 'Maxwellian velocity', [1,3]);


% ve=[1E6*ones(size(ve,1),1) zeros(size(ve,1),1) zeros(size(ve,1),1)];

%%��ʼ�ٶȷֲ�ֱ�ӹ涨Ϊ-dt/2ʱ�̵�ֵ��ֱ�ӿ�ʼ������

%���ȷֲ��ڿռ�Lx�ڣ�����һ������Ŷ�
position_e=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_e,1] );%source region uniform random
position_Hp=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_Hp,1] );%entire domain uniform+noise
position_H2p=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_H2p,1] );
position_H3p=get_position_init(simulation, 'entire domain uniform+noise', [simulation.num0_macro_H3p,1] );%entire domain uniform+noise
position_Hn=0.5*simulation.Lx;

ae=zeros(simulation.num0_macro_e,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
aHp=zeros(simulation.num0_macro_Hp,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
aH2p=zeros(simulation.num0_macro_H2p,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
aH3p=zeros(simulation.num0_macro_H3p,1);%���ٶȳ�ֵ ���Ӽ��ٶ�
aHn=zeros(1,1);%���ٶȳ�ֵ ���Ӽ��ٶ�


% if simulation.particle_boundaries_type==2
%     simulation.num0_macro_e=100* simulation.num0_macro_e;
%     simulation.num0_macro_Hp=100*simulation.num0_macro_Hp;
%     simulation.num0_macro_H2p=100*simulation.num0_macro_H2p;
% end

% if simulation.particle_boundaries_type==2
%     simulation.weight=simulation.ne0*(simulation.source_region(2)-simulation.source_region(1))/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
%     simulation.weight_Hp=simulation.weight;
%     simulation.weight_H2p=simulation.weight;
% end


Be=simulation.B0*Bfun(position_e);
BHp=simulation.B0*Bfun(position_Hp);
BH2p=simulation.B0*Bfun(position_H2p);
BH3p=simulation.B0*Bfun(position_H3p);

%--------������
[ u, A, b_extra ] = get_field_init( simulation );%��ɢ���ɷ���
rho=zeros(simulation.num_grid_point,1);%����������ܶȳ�ֵ
E=zeros(simulation.num_grid_point+1,1);%��ڵ�糡��ֵ  NG-1+����ǽ�ڵİ���
B=zeros(simulation.num_grid_point+1,1);%��ڵ�ų���ֵ
%--------�������
j_average_counter=average_timesteps;
usum=0;
rhoesum=0;
rhoPsum=0;
rhosum=0;
record_flag=0;
u_record=zeros(simulation.num_grid_point,1);
%--------��ʼ��--------------------------------------------------------------------------------------

h_fig1 = figure('Unit','Normalized','position',...
    [0.02 0.1 0.8 0.8]); %ʵʱ���ʹ�õĴ󴰿�
h_fig2_EEPF=figure;

video = VideoWriter('MF_cool_down.avi');
open(video);
video_flag=0;%0��һֱ���ɶ��������1�ͰѶ����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------��ѭ����ʼ-----------------------------------------------------------------------------------
% tic
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
    clear vp vn
    
    %%%%%%%%%%%%%%%%%%
    position_e=position_e+ve(:,1)*ddt;
    
    
    %     vHp(:,1)=vHp(:,1)+aHp*ddt;%��������������λ��
    
    vHpn(:,1)=vHp(:,1)+aHp*ddt/2;
    vHpn(:,3)=vHp(:,3);
    tHp=constants.q_m_ratio_Hp*BHp*ddt/2;
    sHp=2*tHp./(1+tHp.^2);
    
    vpHp(:,1)=vHpn(:,1)-(vHpn(:,3)+vHpn(:,1).*tHp).*sHp;
    vpHp(:,3)=vHpn(:,3)+(vHpn(:,1)-vHpn(:,3).*tHp).*sHp;
    vHp(:,1)=vpHp(:,1)+aHp*ddt/2;
    vHp(:,3)=vpHp(:,3);
    clear vpHp vHpn
    
    position_Hp=position_Hp+vHp(:,1)*ddt;
    
    
    %     vH2p(:,1)=vH2p(:,1)+aH2p*ddt;%��������������λ��
    vH2pn(:,1)=vH2p(:,1)+aH2p*ddt/2;
    vH2pn(:,3)=vH2p(:,3);
    tH2p=constants.q_m_ratio_H2p*BH2p*ddt/2;
    sH2p=2*tH2p./(1+tH2p.^2);
    
    vpH2p(:,1)=vH2pn(:,1)-(vH2pn(:,3)+vH2pn(:,1).*tH2p).*sH2p;
    vpH2p(:,3)=vH2pn(:,3)+(vH2pn(:,1)-vH2pn(:,3).*tH2p).*sH2p;
    vH2p(:,1)=vpH2p(:,1)+aH2p*ddt/2;
    vH2p(:,3)=vpH2p(:,3);
    clear vpH2p vH2pn
    
    position_H2p=position_H2p+vH2p(:,1)*ddt;
    
    %     vH3p(:,1)=vH3p(:,1)+aH3p*ddt;%��������������λ��
    vH3pn(:,1)=vH3p(:,1)+aH3p*ddt/2;
    vH3pn(:,3)=vH3p(:,3);
    tH3p=constants.q_m_ratio_H3p*BH3p*ddt/2;
    sH3p=2*tH3p./(1+tH3p.^2);
    
    vpH3p(:,1)=vH3pn(:,1)-(vH3pn(:,3)+vH3pn(:,1).*tH3p).*sH3p;
    vpH3p(:,3)=vH3pn(:,3)+(vH3pn(:,1)-vH3pn(:,3).*tH3p).*sH3p;
    vH3p(:,1)=vpH3p(:,1)+aH3p*ddt/2;
    vH3p(:,3)=vpH3p(:,3);
    clear vpH3p vH3pn
    position_H3p=position_H3p+vH3p(:,1)*ddt;
    
    vHn(:,1)=vHn(:,1)+aHn*ddt;%��������������λ��
    position_Hn=position_Hn+vHn(:,1)*ddt;
    
    %------�ƶ����ӽ���-----------------------------------------------------------------------------------------------
    
    %------����Խ�紦��------------------------------------------------------------------------------------------------
    switch simulation.particle_boundaries_type
        case 0 % ���ڱ߽�
            % -------���ڱ߽�------------------------------------------------------
            tempFlag=(position_e<0);
            position_e(tempFlag)=position_e(tempFlag)+simulation.Lx;
            tempFlag=(position_e>simulation.Lx);
            position_e(tempFlag)=position_e(tempFlag)-simulation.Lx;
            tempFlag=(position_Hp<0);
            position_Hp(tempFlag)=position_Hp(tempFlag)+simulation.Lx;
            tempFlag=(position_Hp>simulation.Lx);
            position_Hp(tempFlag)=position_Hp(tempFlag)-simulation.Lx;
            tempFlag=(position_H2p<0);
            position_H2p(tempFlag)=position_H2p(tempFlag)+simulation.Lx;
            tempFlag=(position_H2p>simulation.Lx);
            position_H2p(tempFlag)=position_H2p(tempFlag)-simulation.Lx;
            
            tempFlag=(position_H3p<0);
            position_H3p(tempFlag)=position_H3p(tempFlag)+simulation.Lx;
            tempFlag=(position_H3p>simulation.Lx);
            position_H3p(tempFlag)=position_H3p(tempFlag)-simulation.Lx;
            
            tempFlag=(position_Hn<0);
            position_Hn(tempFlag)=position_Hn(tempFlag)+simulation.Lx;
            tempFlag=(position_Hn>simulation.Lx);
            position_Hn(tempFlag)=position_Hn(tempFlag)-simulation.Lx;
        case 1 % �ɶ���ע��+�Ȼ�
            
            % -------�ɶ���ע��Դ�� pair re-injection into source region--------------
            %if(source_density<(1.5*simulation.ne0))
            if(length(position_Hp)+length(position_H2p)+length(position_H3p)<simulation.num0_macro_e)
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%Ĩ����������ĵ���
                position_e(tempFlag)=[];
                
                k1=(position_Hp<0|position_Hp>simulation.Lx);%ȡ������Χ�����ӵı��
                tempsum1=sum(k1);
                
                if(tempsum1>0)%����ǿ�
                    vHp(k1,:)=normrnd(0,vHpth,[tempsum1,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
                    % ��ע������λ�þ�������ֲ���Դ����
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum1,1);
                    position_Hp(k1)=temp_matrix;
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
                
                k3=(position_H3p<0|position_H3p>simulation.Lx);%ȡ������Χ�����ӵı��
                tempsum3=sum(k3);
                
                if(tempsum3>0)%����ǿ�
                    vH3p(k3,:)=normrnd(0,vH3pth,[tempsum3,3]);%xyz ��ע�������ٶ�Maxwell�ֲ�
                    % ��ע������λ�þ�������ֲ���Դ����
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(tempsum3,1);
                    position_H3p(k3)=temp_matrix;
                    position_e=[position_e; temp_matrix];%�ϲ���ԭ��������
                    ve=[ve; normrnd(0,veth,[tempsum3,3])];
                end
                
                tempFlag5=(position_Hn<=0|position_Hn>=simulation.Lx);
                if length(vHn)>3
                    
                    vHn(tempFlag5,:)=[];%Ĩ�����������Hn
                    position_Hn(tempFlag5)=[];
                else
                    source_region_Lx=simulation.source_region(2)-simulation.source_region(1);
                    position_Hn(tempFlag5)=simulation.source_region(1)+source_region_Lx*rand(length(tempFlag5));
                end
            else
                tempFlag=(position_e<0|position_e>simulation.Lx);
                ve(tempFlag,:)=[];%Ĩ����������ĵ���
                position_e(tempFlag)=[];
                
                tempFlag2=(position_Hp<0|position_Hp>simulation.Lx);
                vHp(tempFlag2,:)=[];%Ĩ�����������H+
                position_Hp(tempFlag2)=[];
                
                tempFlag3=(position_H2p<0|position_H2p>simulation.Lx);
                vH2p(tempFlag3,:)=[];%Ĩ�����������H2+
                position_H2p(tempFlag3)=[];
                
                tempFlag4=(position_H3p<0|position_H3p>simulation.Lx);
                vH3p(tempFlag4,:)=[];%Ĩ�����������H2+
                position_H3p(tempFlag4)=[];
                
                tempFlag5=(position_Hn<=0|position_Hn>=simulation.Lx);
                if length(vHn)>3
                    
                    vHn(tempFlag5,:)=[];%Ĩ�����������Hn
                    position_Hn(tempFlag5)=[];
                else
                    source_region_Lx=simulation.source_region(2)-simulation.source_region(1);
                    position_Hn(tempFlag5)=simulation.source_region(1)+source_region_Lx*rand(length(tempFlag5));
                end
                
            end
            
            
            %--------thermalization �Ȼ�  ÿ10��һ��
            % ������������ڱ߽�ʹ��
            if rem(ti,100)==1
                flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
                ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
            end
            
        case 2 %�ɶԺ�ע��
            tempFlag=(position_e<0|position_e>simulation.Lx);
            ve(tempFlag,:)=[];%Ĩ����������ĵ���
            position_e(tempFlag)=[];
            
            tempFlag2=(position_Hp<0|position_Hp>simulation.Lx);
            vHp(tempFlag2,:)=[];%Ĩ�����������H+
            position_Hp(tempFlag2)=[];
            
            tempFlag3=(position_H2p<0|position_H2p>simulation.Lx);
            vH2p(tempFlag3,:)=[];%Ĩ�����������H2+
            position_H2p(tempFlag3)=[];
            
            tempFlag4=(position_H3p<0|position_H3p>simulation.Lx);
            vH3p(tempFlag4,:)=[];%Ĩ�����������H2+
            position_H3p(tempFlag4)=[];
            
            tempFlag5=(position_Hn<0|position_Hn>simulation.Lx);
            if length(vHn)>3
                
                vHn(tempFlag5,:)=[];%Ĩ�����������H2+
                position_Hn(tempFlag5)=[];
            else
                source_region_Lx=simulation.source_region(2)-simulation.source_region(1);
                position_Hn(tempFlag5)=simulation.source_region(1)+source_region_Lx*rand(length(tempFlag5));
            end
            %             source_density=0;
            
            if source_density<(1E17)
                if rem(ti,10)==1
                    temp_matrix=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime,1);
                    position_e=[position_e; temp_matrix];%�ϲ���ԭ��������
                    ve=[ve; normrnd(0,veth,[inject_num_pertime,3])];
                    
                    temp_matrix1=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime*Hp_part,1);
                    position_Hp=[position_Hp; temp_matrix1];%�ϲ���ԭ��������
                    vHp=[vHp; normrnd(0,vHpth,[inject_num_pertime*Hp_part,3])];
                    
                    temp_matrix2=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime*H2p_part,1);
                    position_H2p=[position_H2p; temp_matrix2];%�ϲ���ԭ��������
                    vH2p=[vH2p; normrnd(0,vH2pth,[inject_num_pertime*H2p_part,3])];
                    
                    temp_matrix3=simulation.source_region(1)+(simulation.source_region(2)-simulation.source_region(1))*rand(inject_num_pertime*H3p_part,1);
                    position_H3p=[position_H3p; temp_matrix3];%�ϲ���ԭ��������
                    vH3p=[vH3p; normrnd(0,vH3pth,[inject_num_pertime*H3p_part,3])];
                end
            end
            
            %             %--------thermalization �Ȼ�  ÿ10��һ��
            %             % ������������ڱ߽�ʹ��
            %             if rem(ti,100)==1
            %                 flag=(simulation.source_region(1)<position_e & position_e<simulation.source_region(2));
            %                 ve(flag,:)=normrnd(0,veth,[sum(flag),3]);
            %             end
            
        otherwise
            throw(get_exception('error','no such type of particle boundaries.'))
    end
    %------����Խ�紦�����------------------------------------------------------------------------------------------
    
    %------����MCC����------------------------------------------------------------------------------------------
    
    %     if rem(ti,4)==1
    %         disrec_idx_e=randperm(length(position_e),4);
    %         disrec_idx_Hp=randperm(length(position_Hp),1);
    %         disrec_idx_H2p=randperm(length(position_H2p),2);
    %         disrec_idx_H3p=randperm(length(position_H3p),1);
    %
    %         position_e(disrec_idx_e)=[];
    %         ve(disrec_idx_e,:)=[];
    %         position_Hp(disrec_idx_Hp)=[];
    %         vHp(disrec_idx_Hp,:)=[];
    %          position_H2p(disrec_idx_H2p)=[];
    %         vH2p(disrec_idx_H2p,:)=[];
    %         position_H3p(disrec_idx_H3p)=[];
    %         vH3p(disrec_idx_H3p,:)=[];
    %     end
    
    [collision_index,collision_type]=null_collision( ve,constants,mcc_para,Xsec,'H2' );
    [colli_result,ion_react_order]=collision_process(Xsec, mcc_para,ve,constants,collision_index,collision_type ,'H2');
    
    if colli_result.happen==1
        ve(colli_result.ion_index,:)=colli_result.ion_inc_v;
        ve(colli_result.ela_index,:)=colli_result.ela_v;
        ve(colli_result.exc_index,:)=colli_result.exc_v;
        %             ve(colli_result.att_index,:)=[];
        ve=[ve; colli_result.ion_gen_v ];
        position_e=[position_e; position_e(colli_result.ion_index)];
        %
        for i=1:length(ion_react_order)
            if ion_react_order(i)==1
                vH2p=[vH2p; normrnd(0,vH2pth,[1 , 3]);  ];
                position_H2p=[position_H2p; position_e(colli_result.ion_index(i))];
            else
                vHp=[vHp; normrnd(0,vHpth,[1 , 3]);  ];
                position_Hp=[position_Hp; position_e(colli_result.ion_index(i))];
            end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [collision_index_H,collision_type_H]=null_collision( ve,constants,mcc_para,Xsec_H,'H' );
    [colli_result,~]=collision_process(Xsec_H, mcc_para,ve,constants,collision_index_H,collision_type_H,'H' );
    
    if colli_result.happen==1
        ve(colli_result.ion_index,:)=colli_result.ion_inc_v;
        ve(colli_result.ela_index,:)=colli_result.ela_v;
        ve(colli_result.exc_index,:)=colli_result.exc_v;
        %             ve(colli_result.att_index,:)=[];
        ve=[ve; colli_result.ion_gen_v ];
        position_e=[position_e; position_e(colli_result.ion_index)];
        %
        vHp=[vHp; normrnd(0,vHpth,[size(colli_result.ion_gen_v,1) , 3]);  ];
        position_Hp=[position_Hp; position_e(colli_result.ion_index)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %Xsec_Hp_H2,Xsec_Hp_H,Xsec_H2p_H2,Xsec_H2p_H,Xsec_H3p_H2,Xsec_H3p_H,Xsec_Hn_H2,Xsec_Hn_H
    %%
    %Hp_H2
    [ C_index_heavy_Hp_H2,C_type_Hp_H2] = null_collision_heavy_par( vHp,constants,mcc_para,Xsec_Hp_H2,'Hp+H2');
    if ~isempty(C_index_heavy_Hp_H2)%1 ion; 2 att; 3 ela; 4~6exc
        C3=C_index_heavy_Hp_H2(C_type_Hp_H2==3);
        C4=C_index_heavy_Hp_H2(C_type_Hp_H2>3);
        vp3=vHp(C3,:);
        vp4=vHp(C4,:);
        
        HpH2ela=HpH2ela+length(C3);
        HpH2exc=HpH2exc+length(C4);
        
        
        vHp(C3,:)=Ion_ela_process(vp3);%ela
        vHp(C4,:) = Ion_exc_process(C_type_Hp_H2(C_type_Hp_H2>3)-mcc_para.ela_index_Hp_H2,vp4,constants.q_m_ratio_Hp,Xsec_Hp_H2.excThresh);%exc
        
        %         vHp(collision_index_heavy,:)=sqrt(2)/2*vecnorm(vHp(collision_index_heavy,:),2,2).*temp_p;
    end
    %%
    %Hp_H
    [ C_index_heavy_Hp_H,C_type_Hp_H] = null_collision_heavy_par( vHp,constants,mcc_para,Xsec_Hp_H,'Hp+H');
    if ~isempty(C_index_heavy_Hp_H)%1 ion; 2 CX; 3 ela; 4 exc
        C2=C_index_heavy_Hp_H(C_type_Hp_H==2);
        C3=C_index_heavy_Hp_H(C_type_Hp_H==3);
        
        %          C2=[];
        
        vp2=vHp(C2,:);
        vp3=vHp(C3,:);
        
        
        HpHCX=HpHCX+length(C2);
        HpHela=HpHela+length(C3);
        
        vHp(C2,:)=0*vp2;%CX
        vHp(C3,:)=Ion_ela_process(vp3);%ela
    end
    %%
    %H2p_H2
    [ C_index_heavy_H2p_H2,C_type_H2p_H2] = null_collision_heavy_par( vH2p,constants,mcc_para,Xsec_H2p_H2,'H2p+H2');
    if ~isempty(C_index_heavy_H2p_H2)%1 ion; 2 CX; 3 PX; 4 ela; 5 exc
        C2=C_index_heavy_H2p_H2(C_type_H2p_H2==2);
        C3=C_index_heavy_H2p_H2(C_type_H2p_H2==3);
        vp2=vH2p(C2,:);
        vp3=vH2p(C3,:);
        
        H2pH2CX=H2pH2CX+length(C2);
        H2pH2PX=H2pH2PX+length(C3);
        
        vH2p(C2,:)=0*vp2;%CX
        
        vH3p=[vH3p; normrnd(0,vH3pth,[length(C3) , 3]);  ];
        position_H3p=[position_H3p;position_H2p(C3,:);];
        vH2p(C3,:)=[];%PX
        position_H2p(C3,:)=[];
        
    end
    
    %%
    %H2p_H
    [ C_index_heavy_H2p_H,C_type_H2p_H] = null_collision_heavy_par( vH2p,constants,mcc_para,Xsec_H2p_H,'H2p+H');
    if ~isempty(C_index_heavy_H2p_H)%1 ion; 2  ; 3 ela; 4 exc
        
        C3=C_index_heavy_H2p_H(C_type_H2p_H==3);
        
        H2pHela=H2pHela+length(C3);
        
        vp3=vH2p(C3,:);
        
        vH2p(C3,:)=Ion_ela_process(vp3);%ela
    end
    %%
    %H3p_e
    
    %%
    %------����MCC�������------------------------------------------------------------------------------------------
    
    %-----------������----------------------------------------------------------------------------------
    %     rhoea=zeros(simulation.num_grid_point,1);%��������������������ܶ� ����
    %     rhoeb=zeros(simulation.num_grid_point,1);
    %     rhoaH=zeros(simulation.num_grid_point,1);%����
    %     rhobH=zeros(simulation.num_grid_point,1);
    %     rhoaH2=zeros(simulation.num_grid_point,1);%����
    %     rhobH2=zeros(simulation.num_grid_point,1);
    
    enearx=floor(position_e/ddx)+1;%ÿ������������������� �������
    eassigndx=position_e-(enearx-1)*ddx;%ÿ�����Ӿ���������������ľ���
    Eenearx=floor((position_e+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
    Eedx=position_e-(Eenearx*ddx-3/2*ddx);%ÿ�����Ӿ����������İ�����ľ���
    
    Hpnearx=floor(position_Hp/ddx)+1; %H+
    Hpassigndx=position_Hp-(Hpnearx-1)*ddx;
    EHpnearx=floor((position_Hp+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
    EHpdx=position_Hp-(EHpnearx*ddx-3/2*ddx);
    
    H2pnearx=floor(position_H2p/ddx)+1; %H2+
    H2passigndx=position_H2p-(H2pnearx-1)*ddx;
    EH2pnearx=floor((position_H2p+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
    EH2pdx=position_H2p-(EH2pnearx*ddx-3/2*ddx);
    
    H3pnearx=floor(position_H3p/ddx)+1; %H2+
    H3passigndx=position_H3p-(H3pnearx-1)*ddx;
    EH3pnearx=floor((position_H3p+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
    EH3pdx=position_H3p-(EH3pnearx*ddx-3/2*ddx);
    
    Hnnearx=floor(position_Hn/ddx)+1; %H2+
    Hnassigndx=position_Hn-(Hnnearx-1)*ddx;
    EHnnearx=floor((position_Hn+0.5*ddx)/ddx)+1;%ÿ���������ڰ���������� ����
    EHndx=position_Hn-(EHnnearx*ddx-3/2*ddx);
    
    rhoea=accumarray(enearx,ddx-eassigndx,[simulation.num_grid_point 1]);
    rhoeb=accumarray(enearx+1,eassigndx,[simulation.num_grid_point 1]);
    
    rhoe=(rhoea+rhoeb)*(-constants.e)/(ddx^2)*simulation.weight;
    %rhoP=rhoP*0;
    rhoaHp=accumarray(Hpnearx,ddx-Hpassigndx,[simulation.num_grid_point 1]);
    rhobHp=accumarray(Hpnearx+1,Hpassigndx,[simulation.num_grid_point 1]);
    
    rhoHp=(rhoaHp+rhobHp)*constants.e/(ddx^2)*simulation.weight_Hp;
    
    rhoaH2p=accumarray(H2pnearx,ddx-H2passigndx,[simulation.num_grid_point 1]);
    rhobH2p=accumarray(H2pnearx+1,H2passigndx,[simulation.num_grid_point 1]);
    
    rhoH2p=(rhoaH2p+rhobH2p)*constants.e/(ddx^2)*simulation.weight_H2p;
    
    rhoaH3p=accumarray(H3pnearx,ddx-H3passigndx,[simulation.num_grid_point 1]);
    rhobH3p=accumarray(H3pnearx+1,H3passigndx,[simulation.num_grid_point 1]);
    
    rhoH3p=(rhoaH3p+rhobH3p)*constants.e/(ddx^2)*simulation.weight_H3p;
    
    rhoaHn=accumarray(Hnnearx,ddx-Hnassigndx,[simulation.num_grid_point 1]);
    rhobHn=accumarray(Hnnearx+1,Hnassigndx,[simulation.num_grid_point 1]);
    
    rhoHn=(rhoaHn+rhobHn)*constants.e/(ddx^2)*simulation.weight_Hn;
    
    %     rhoHn=0;
    rho=rhoe+rhoHp+rhoH2p+rhoH3p+rhoHn;
    %-----------�����ɽ���------------------------------------------------------------
    
    source_density=-mean(rhoe(source_left_bdidx:source_right_bdidx))/constants.e;
    
    %--------------���糡---------------------------------------------------------------
    b=get_b( rho,b_extra,simulation );
    u=get_u( u,A,b,'direct inverse');
    E=get_E_at_half_node( u, simulation );
    %         E=1000*ones(size(E,1),1);
    %--------------���糡����---------------------------------------------------------------
    
    % TODO�������E_particle
    ae=-((ddx-Eedx).*E(Eenearx)+Eedx.*E(Eenearx+1))*constants.e/(constants.me)/ddx;
    aHp=((ddx-EHpdx).*E(EHpnearx)+EHpdx.*E(EHpnearx+1))*constants.e/(constants.mH)/ddx;
    aH2p=((ddx-EH2pdx).*E(EH2pnearx)+EH2pdx.*E(EH2pnearx+1))*constants.e/(constants.mH2)/ddx;
    aH3p=((ddx-EH3pdx).*E(EH3pnearx)+EH3pdx.*E(EH3pnearx+1))*constants.e/(constants.mH3)/ddx;
    aHn=((ddx-EHndx).*E(EHnnearx)+EHndx.*E(EHnnearx+1))*constants.e/(constants.mH)/ddx;
    
    %     Be=(ddx-Eedx).*B(Eenearx)+Eedx.*B(Eenearx+1);%Ӧ���ñ��ʽ��д
    Be=simulation.B0*Bfun(position_e);
    BHp=simulation.B0*Bfun(position_Hp);
    BH2p=simulation.B0*Bfun(position_H2p);
    BH3p=simulation.B0*Bfun(position_H3p);
    
    % ��ѹ�ۼӣ�Ϊ  ʵʱ�����ȡƽ��  ��׼��
    if rem(ti+j_average_counter,display_timesteps)==1 %ti����display��ǰaverage_timesteps��Χ��
        usum=usum+u;
        rhoesum=rhoesum+rhoe;
        rhoPsum=rhoPsum+rhoHp+rhoH2p+rhoH3p;
        rhosum=rhosum+rho;
        j_average_counter=j_average_counter-1;
        if j_average_counter==-1 %ti�ѳ���display��ǰaverage_timesteps��Χ
            j_average_counter=average_timesteps; %����j_average_counter
        end
    end
    if rem(ti,10)==1
        ve_static1=[ve_static1;ve(position_e>OB_position1-ddx & position_e<OB_position1+ddx,:)];
        ve_static2=[ve_static2;ve(position_e>OB_position2-ddx & position_e<OB_position2+ddx,:)];
        ve_static3=[ve_static3;ve(position_e>OB_position3-ddx & position_e<OB_position3+ddx,:)];
    end
    %-----------ʵʱ���----------------------------------------------------------------------------------
    if rem(ti-1,display_timesteps)==0%ÿ��dptime��ʾһ��
        fprintf('��ǰʱ������%d/%d \n',ti,simulation.all_timesteps)
        
        display_i=ceil(ti/display_timesteps);
        e_num(display_i)=size(position_e,1);%�����ڵĵ��������仯
        Hp_num(display_i)=size(position_Hp,1);
        H2p_num(display_i)=size(position_H2p,1);
        H3p_num(display_i)=size(position_H3p,1);
        
        e_v_mag=sqrt(ve(:,1).^2+ve(:,2).^2+ve(:,3).^2);
        Hp_v_mag=sqrt(vHp(:,1).^2+vHp(:,2).^2+vHp(:,3).^2);
        H2p_v_mag=sqrt(vH2p(:,1).^2+vH2p(:,2).^2+vH2p(:,3).^2);
        H3p_v_mag=sqrt(vH3p(:,1).^2+vH3p(:,2).^2+vH3p(:,3).^2);
        
        
        e_vAve(display_i)=mean(e_v_mag);%x�����ϵ�ƽ�����ʱ仯
        Hp_vAve(display_i)=mean(Hp_v_mag);
        H2p_vAve(display_i)=mean(H2p_v_mag);
        H3p_vAve(display_i)=mean(H3p_v_mag);
        
        figure(h_fig1)
        subplot(4,2,1);
        %         plot_no_versus_position( position_e(1:floor(size(position_e,1)/1000):end),position_H(1:floor(size(position_H,1)/1000):end),position_H2p(1:floor(size(position_H2p,1)/1000):end),simulation );
        plot_no_versus_position( position_e,position_Hp,position_H2p,position_H3p,simulation );
        subplot(4,2,2);
        plot_num_versus_timestep( e_num,Hp_num,H2p_num,H3p_num,display_timesteps )
        subplot(4,2,3);
        plot_2u1E_versus_x( u, ti, usum, average_timesteps, E, simulation )
        subplot(4,2,4);
        plot_Ek_versus_timestep( e_vAve,Hp_vAve,H2p_vAve,H3p_vAve, constants.q_m_ratio_e, constants.q_m_ratio_Hp, 'e','Hp','H2p','H3p', display_timesteps )
        ylabel('ƽ��Ek [eV]'); %ʵ������3��x��ƽ������
        subplot(4,2,5);%��ɷֲ�
        plot_density_versus_x( rhoesum/(average_timesteps+1), rhoPsum/(average_timesteps+1),  rhosum/(average_timesteps+1), simulation )
        temp_time_rhoP(:,(ti-1)/display_timesteps+1)=rhoPsum/(average_timesteps+1);
        subplot(4,2,6);
        plot_v_versus_position( ve(:,1), position_e,simulation,'e','b' )
        subplot(4,2,8);
        plot_v_versus_position(vHp(:,1),position_Hp,simulation,'Hp','r' )
        hold on
        plot_v_versus_position(vH2p(:,1),position_H2p,simulation,'H2p','g' )
        hold on
        plot_v_versus_position(vH3p(:,1),position_H3p,simulation,'H3p','m' )
        legend('Hp','H2p','H3p');
        
        subplot(4,2,7);
        plot_vineV_versus_position( e_v_mag, Hp_v_mag,H2p_v_mag,H3p_v_mag, position_e,position_Hp,position_H2p,position_H3p,constants.q_m_ratio_e, constants.q_m_ratio_Hp, 'e','Hp','H2p','H3p', simulation )
        
        figure(h_fig2_EEPF)%��ʾEEPF
        plot_position_EEPF( ve_static1,constants)
        plot_position_EEPF( ve_static2,constants)
        plot_position_EEPF( ve_static3,constants)
        ylabel('EEPF');
        hold off
        ve_static1=[];
        ve_static2=[];
        ve_static3=[];
        
        drawnow;
        
        if video_flag==0
            M(display_i) = getframe(h_fig1);
            writeVideo(video,M(display_i) );
            %         hold off
        else
            close(video);
        end
        
        usum=0; %����
        rhoesum=0;
        rhoPsum=0;
        rhosum=0;
        
        HpH2ela
        HpH2exc
        HpHCX
        HpHela
        
        [H2pH2CX  H2pH2PX H2pHela]
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
x_final=(0:1:(simulation.num_grid_point-1))*ddx;
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

% toc
%% ����
% ��ȡ�洢��������