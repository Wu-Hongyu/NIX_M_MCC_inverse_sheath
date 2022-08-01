function [result] = check_simulation_parameters( simulation, type )
% �����������Ч��
% ʹ��: % ���б�������������ʾ�޸�simulation�������ٴ����б�����
% type: ����ϸ�̶ȡ�1ȫ����Ϊ���棨����ֹ���򣩣�3ȫ��������ֹ����
% ʵ�֣�ÿ�ּ���׳�ָ�������쳣���ֱ�try/catch�Ա��ⲻ����ʱҲֻ����һ��check
constants=get_constants();
%% ���results��������Ļ
result.flag_check=true; %checkʧ��������
result.debye_length=get_debye_length(simulation.ne0,simulation.Te0);
result.omega_pe=get_omega_pe(simulation.ne0);
result.T_pe=2*pi/result.omega_pe;
result.N_D=get_N_D(simulation.ne0,simulation.Te0); %�°�����������Ŀ
result.vth_e=sqrt(-constants.q_m_ratio_e*simulation.Te0);
result.flag_magnetized=false;
if simulation.B0>0
    result.omega_ce=-constants.q_m_ratio_e*simulation.B0;
    result.T_ce=2*pi/result.omega_ce;
    result.r_ce=result.vth_e/result.omega_ce;
    if result.r_ce<simulation.Lx
        result.flag_magnetized=true;
    end
end
% result.ppc=

fprintf('[INFO] ʱ�������߶�\n')
fprintf('���ӵ�������Ƶ�� ��_pe = %.2e rad/s\n',result.omega_pe)
fprintf('ģ���򳤶� Lx = %.2e m;  ',simulation.Lx)
fprintf('�°ݳ��� ��_De = %.2e m;  ',result.debye_length)
fprintf('�ռ䲽�� ��x = %.2e m;  ',simulation.dx)
fprintf('����ڵ��� N_node = %d \n',simulation.num_grid_point)
if 3==simulation.field_boundaries_type
    fprintf('�޽�������壺ģ����ʹ�����ڱ߽�\n')
else
    fprintf('�н�������壺ģ����ʹ�÷����ڱ߽�\n')
end
fprintf('����ʱ�� t_end = %.2e s;  ',simulation.end_time)
fprintf('���Ӿ��������� T_pe = %.2e s;  ',result.T_pe)
fprintf('ʱ�䲽�� ��t = %.2e s;  ',simulation.dt)
fprintf('Ԥ����ʱ���� N_step = %d \n',simulation.all_timesteps)
fprintf('�������ٶ� vth_e = %.2e \n',result.vth_e)
fprintf('�°�����������Ŀ N_D = %.2e \n',result.N_D)
fprintf('�����ԴŸ�Ӧǿ�� B0 = %d T\n',simulation.B0)
if simulation.B0>0
    fprintf('��Ī��Ƶ�� ��_ce = %d rad/s;  ',result.omega_ce)
    fprintf('�������� T_ce = %d s;  ',result.T_ce)
    fprintf('��Ī���뾶 r_ce = %d m \n;  ',result.r_ce)
end
% fprintf('�����Ӵ��������� weighting = %.2e \n',result.N_D)
% fprintf('��ʼ���Ӻ������� N_macroe = %d \n',result.N_D)
% fprintf('ÿ��Ԫ�������� PPC = %d \n',result.N_D)

%% ����������
fprintf('[INFO] ����������\n')
ratio_Lx_debye_ratio=simulation.Lx/result.debye_length;
fprintf('Lx = %.2e ��_De\n',ratio_Lx_debye_ratio)
try
    if 3==simulation.field_boundaries_type
        if ratio_Lx_debye_ratio<1
            throw(get_exception('warn','Lx < ��_De���޷�ʶ��°���'))
        end
    else %�н��������
        if ratio_Lx_debye_ratio<1
            throw(get_exception('error','Lx < ��_De�����������о�Ҫ��Lx>>��_De'))
        elseif ratio_Lx_debye_ratio<10
            %     fprintf('[WARN]: Lx < 10*��_De\n')
            throw(get_exception('warn','��_De <= Lx < 10*��_De�����������о�Ҫ��Lx>>��_De'))
        end
    end
catch e
    treate_exception(e,type)
end

try
    if 0==simulation.Te0
        throw(get_exception('warn','Te0=0, ��N_D = 0'))
    elseif result.N_D<1
        throw(get_exception('error','N_D < 1�����������о�Ҫ��N_D>>>1'))
    elseif result.N_D<100
        throw(get_exception('warn','1 <= N_D < 100�����������о�Ҫ��N_D>>>1'))
    end
catch e
    treate_exception(e,type)
end
% TODO������������������ӵ�ƽ����ײ���ڣ����������оݣ�����MCCģ��ȷ������ʵ��

ratio_end_time_T_pe=simulation.end_time/result.T_pe;
fprintf('t_end = %.2e T_pe\n',ratio_end_time_T_pe)
try
    if ratio_end_time_T_pe<1
        throw(get_exception('warn','t_end < T_pe���޷���ԭ������������Ϣ'))
        % ����ʱ�䲽����error����˴˴�ֻwarn
    elseif ratio_end_time_T_pe<10
        throw(get_exception('warn','t_end < 10*T_pe�����Ի�ԭ������������Ϣ'))
    end
catch e
    treate_exception(e,type)
end

% �Ż�
if result.flag_magnetized
    disp('���ӱ��Ż�')
    ratio_end_time_T_ce=simulation.end_time/result.T_ce;
    fprintf('t_end = %.2e T_ce\n',ratio_end_time_T_ce)
    try
        if ratio_end_time_T_ce<1
            throw(get_exception('warn','t_end < T_ce���޷���ԭ��Ī���˶���Ϣ'))
            % ����ʱ�䲽����error����˴˴�ֻwarn
        elseif ratio_end_time_T_ce<10
            throw(get_exception('warn','t_end < 10*T_ce�����Ի�ԭ��Ī���˶���Ϣ'))
        end
    catch e
        treate_exception(e,type)
    end
    
    ratio_Lx_r_ce=simulation.Lx/result.r_ce;
    fprintf('Lx = %.2e r_ce\n',ratio_Lx_r_ce)
    if ratio_Lx_r_ce>10
        fprintf('Lx > 10*r_ce, \n')
    end
end

%% �����ֵ����
fprintf('[INFO] �����ֵ����-��ɢ��������\n')
% ������ֵ�Լ��ȵĿռ䲽���������� dx<�֦�D
ratio_dx_debye_length=simulation.dx/result.debye_length;
fprintf('��x = %.2e ��_De\n',ratio_dx_debye_length)
try
    if 0==simulation.Te0
        throw(get_exception('warn','Te0=0, ���_De = 0'))
    elseif ratio_dx_debye_length>3
        throw(get_exception('error','��x > 3*��_De���뿼����ֵ����'))
    elseif ratio_dx_debye_length>1
        throw(get_exception('warn','��x > ��_De'))
    end
catch e
    treate_exception(e,type)
end
% �����������Ӿ����񵴵�����������_pe*��t/2<1������Ҫ�� ��p*dt/2<<1
omega_pe_dt_2=result.omega_pe*simulation.dt/2;
fprintf('��_pe*��t/2 = %.2e rad\n',omega_pe_dt_2)
try
    if omega_pe_dt_2>1
        throw(get_exception('error','��_pe*��t/2>1 rad�������������������ڷ�ֵ���'))
    elseif omega_pe_dt_2>0.1
        throw(get_exception('warn','��_pe*��t/2>0.1 rad���뿼����������λ���'))
    end
catch e
    treate_exception(e,type)
end
% ��ʽʱ��FDM���ƫ΢�ַ��̵���������(CFL����)vth_e*��t<��x
ratio_v_dt_dx=result.vth_e*simulation.dt/simulation.dx;
fprintf('vth_e*��t = %.2e ��x\n',ratio_v_dt_dx)
try
    if ratio_v_dt_dx>1
        throw(get_exception('error','vth_e*��t>��x��������CFL����'))
    end
catch e
    treate_exception(e,type)
end
% TODO��result.ppc����paticle_groupȷ������ʵ��
% �ܶ�ͳ�Ƶ�ɢ������������1/sqrt(ppc)����˾���Ҫ��PPC>10
% TODO����ײ��ʱ�䲽����Ҫ�󡣴�MCCģ��ȷ������ʵ��

% TODO: ���Բ�Ӧ�ó��ֵĴ���
% assert particle_group[1]������simulation�е��Ӳ������  %Ҳ���Էŵ�test_particle_group��

if result.flag_check
    fprintf('OK! \n')
end
end

function treate_exception(exception,type)
switch type
    case 1 %warn��error��warning
        warning(exception.message)
    case 2 % ����warn��װ��warning������error��rethrow
        temp_e=get_exception('warn','');
        if exception.identifier==temp_e.identifier
            warning(exception.message)
        else
            result.flag_check=false;
            rethrow(exception)
        end
    case 3 %warn��error��rethrow
        result.flag_check=false;
        rethrow(exception)
end
end