function [ simulation ] = get_simulation( type, varargin )
% ����simulation�ṹ��
% type: 'default', 'auto1'
% varargin: �ɱ����������ͨ��name-value������
% type='auto1': ne0, Te0, end_time, Lx
% simulation: ����ģ�Ͳ��� �ṹ��
%% ȫ������
switch type
    case 'default' %A1DPIC_Ver1.0.m, 2019Montellano
        %%%%%%% ���ã������޸� %%%%%%%%
        % ��������������Ե���Ϊ����
        simulation.ne0=1E16;%���������ܶ�
        simulation.Te0=1;%��λeV
        %         TH=1;%��λeV
        % ʱ�ճ߶�
        simulation.dt=10*10^-12;%ʱ�䲽��
        simulation.end_time=6E-6;%����ʱ���ܳ�
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
        simulation.Lx=0.01;%�������򳤶�
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        simulation.num_grid_point=201;%��������Ŀ
        simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
        % ��ų�
        simulation.field_boundaries_type=0;%��λ�߽���������
        simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
        simulation.B0=0; %�����ԴŸ�Ӧǿ�� T
        % ����
        simulation.particle_boundaries_type=0;%���ӱ߽���������
        %��ѹ
        simulation.pressure=0.6;%��ѹ ��λPa
        simulation.Tgas=300;%����   K
        case 'BATMAN_old_MF'
            % TODO��������չǻģ����Ҫ�����޸�
        % ��������������Ե���Ϊ����
        simulation.ne0=1E18;%���������ܶ�
        simulation.Te0=10;%��λeV
        % ʱ�ճ߶�
        simulation.dt=10*10^-12;%ʱ�䲽��
        simulation.end_time=6E-6;%����ʱ���ܳ�
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
        simulation.Lx=0.01;%�������򳤶�
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        simulation.num_grid_point=150;%��������Ŀ
        simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
        % ��ų�
        simulation.field_boundaries_type=0;%��λ�߽���������
        simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
        simulation.B0=9e-3; %�����ԴŸ�Ӧǿ�� T
        % ����
        simulation.particle_boundaries_type=0;%���ӱ߽���������
         %��ѹ
        simulation.pressure=0.6;%��ѹ ��λPa
        simulation.Tgas=293;%����   K
        %     case 'BACON base'
        
    %% �Զ�����
    %TODO��'auto000'������type���ϲ�����
    case 'auto0' %����ne, Te, Lx, end_time���ϸ�������ɢ��������,��ͨ��type3��check
        % ��������
        p = inputParser;
        validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
        validScalarNotNegNum = @(x) isnumeric(x) && isscalar(x) && ~(x < 0);
        addRequired(p,'type');
        addOptional(p,'ne0',1E16,validScalarPosNum);
        addOptional(p,'Te0',1,validScalarNotNegNum);
        addOptional(p,'end_time',6E-6,validScalarPosNum);
        addOptional(p,'Lx',0.01,validScalarPosNum);
        parse(p,type,varargin{:});
        % p.Results
        % ��������������Ե���Ϊ����
        simulation.ne0=p.Results.ne0;%���������ܶ�
        simulation.Te0=p.Results.Te0;%��λeV
        simulation.end_time=p.Results.end_time;%����ʱ���ܳ�
        simulation.Lx=p.Results.Lx;%�������򳤶�
        % ʱ�ճ߶�
        %%%%%%%%%% ʱ�ղ���
        % ���� ����������г�˶��ľ���Ҫ�� ��p*dt/2<<1
        omega_pe=get_omega_pe(simulation.ne0);
        simulation.dt=(2/omega_pe)*(pi/32); % dt=0.098*(2/��p)=Tpe/32
        if simulation.Te0>0
            % ������ֵ�Լ��� dx<��D
            debye_length=get_debye_length(simulation.ne0,simulation.Te0);
            num_grid=ceil(simulation.Lx/(debye_length)); 
            simulation.num_grid_point=num_grid+1;%��������Ŀ
            simulation.dx=simulation.Lx/num_grid;
            % �����Լ�����������ʱ��CFL������������������Ҫ��һ��������㡣
        else
            num_grid=20;
            simulation.num_grid_point=num_grid+1;%��������Ŀ
            simulation.dx=simulation.Lx/num_grid;   
            % v_th=0������CFL����
        end
%         % ����CFL���� vth_e*��t<��x
%         constants=get_constants();
%         vth_e=sqrt(-constants.q_m_ration_e*simulation.Te0);
%         simulation.dt=min([simulation.dt simulation.dx/vth_e]);
        %%%%%%%%%% ����
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        % ��ų�
        simulation.field_boundaries_type=0;%��λ�߽���������
        simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
        simulation.B0=0; %�����ԴŸ�Ӧǿ�� T
        % ����
        simulation.particle_boundaries_type=0;%���ӱ߽���������

    case 'auto1' %����ne, Te, Lx, end_time������������ɢ��������,��ͨ��type2��check
        % ��������
        p = inputParser;
        validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
        validScalarNotNegNum = @(x) isnumeric(x) && isscalar(x) && ~(x < 0);
        addRequired(p,'type');
        addOptional(p,'ne0',1E16,validScalarPosNum);
        addOptional(p,'Te0',1,validScalarNotNegNum);
        addOptional(p,'end_time',6E-6,validScalarPosNum);
        addOptional(p,'Lx',0.01,validScalarPosNum);
        parse(p,type,varargin{:});
        % p.Results
        % ��������������Ե���Ϊ����
        simulation.ne0=p.Results.ne0;%���������ܶ�
        simulation.Te0=p.Results.Te0;%��λeV
        simulation.end_time=p.Results.end_time;%����ʱ���ܳ�
        simulation.Lx=p.Results.Lx;%�������򳤶�
        % ʱ�ճ߶�
        %%%%%%%%%% ʱ�ղ���
        % ���� ����������г�˶��ľ���Ҫ�� ��p*dt/2<<1
        omega_pe=get_omega_pe(simulation.ne0);
        simulation.dt=(2/omega_pe)*(pi/20); % dt=0.157*(2/��p)=Tpe/20
        if simulation.Te0>0
            % ������ֵ�Լ��� %��������dx<3*��D
            debye_length=get_debye_length(simulation.ne0,simulation.Te0);
            num_grid=ceil(simulation.Lx/(3*debye_length)); 
            simulation.num_grid_point=num_grid+1;%��������Ŀ
            simulation.dx=simulation.Lx/num_grid;
            % �����Լ�����������ʱ��CFL������������������Ҫ��һ��������㡣
        else
            num_grid=20;
            simulation.num_grid_point=num_grid+1;%��������Ŀ
            simulation.dx=simulation.Lx/num_grid;   
            % v_th=0������CFL����
        end
%         % ����CFL���� vth_e*��t<��x
%         constants=get_constants();
%         vth_e=sqrt(-constants.q_m_ration_e*simulation.Te0);
%         simulation.dt=min([simulation.dt simulation.dx/vth_e]);
        %%%%%%%%%% ����
        simulation.all_timesteps=ceil(simulation.end_time/simulation.dt);%��ѭ������
        simulation.source_region=[0.2*simulation.Lx,0.3*simulation.Lx];
        % ��ų�
        simulation.field_boundaries_type=0;%��λ�߽���������
        simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
        simulation.B0=0; %�����ԴŸ�Ӧǿ�� T
        % ����
        simulation.particle_boundaries_type=0;%���ӱ߽���������
    otherwise
        error('Not Done');
end

