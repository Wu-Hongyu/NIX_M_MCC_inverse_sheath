function [ Xsec ] = get_Xsec_e_H( type, flag_display)
% ����ImportLXCat�õ���e-Hԭʼ���ݣ���ȡMCC�������
% ������Χ 1e-2-1e3 eV�ڿɿ� ע����Ϊ���ܵ��Ա��������Բ�ֵ
% ʹ��
% type��e_Hɢ�����������������һ�㼴���ݿ���
% flag_display��0�򲻻�ͼ�������Ӧ���࣬1����ʾ
% ע�⣬���ļ�����.\others\choose_Xsec_e_n_of_hydrogen_plasma
% \compare_Xsec_e_H.mʵ��
% �����ԵıȽ���compare_Xsec_e_H
% �����Ե�test��test_get_Xsec

% % ������ʹ�þֲ�����������ʱ�ֶ�ʹ��
% test_compare_exc_data()
% test_combine_exc_data()
% % 20201217 these two test passed
path_list=get_path_init('test000');
recation_type={['e-H-' type]};  %�����cell
switch type
    case 'CCC-m' % �ϲ�CCCͬ��ֵ������Ӧ���Լ�����ֵ���
        % ����ʹ��CCC (��1)
        % ����ʹ��CCC (��n=2s/p, 3s/p/d, 4s/p/f/d) �ϲ���n=2, 3, 4
        % ����ʹ��CCC (��n=1)
        %%%%%%% ����CCC
        data_dir={[path_list.cross_sections '\e_H_CCC_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        %%% ׼��
        E_head=1e-3; % ����������Χ
        E_end=1e4;
        
        %%%%%%% ��������
        % �Ӻ���ǰcombine������Ҫ�����µ����
        Xsec=combine_exc_data(Xsec, 1, 6:9);
        Xsec=combine_exc_data(Xsec, 1, 3:5);
        Xsec=combine_exc_data(Xsec, 1, 1:2);
        % ���ڴ����񵴣����������Զ����㣬����ڴ�Ԥ���ֶ�����
        Xsec.exc{1}{1}=[[0 0];[Xsec.excThresh{1}{1} 0];Xsec.exc{1}{1}];
        Xsec.exc{1}{2}=[[0 0];[Xsec.excThresh{1}{2} 0];Xsec.exc{1}{2}];
        Xsec.exc{1}{3}=[[0 0];[Xsec.excThresh{1}{3} 0];Xsec.exc{1}{3}];
        
        %%%%%%% ��������
        for E_end_i=logspace(log10(100),log10(E_end),15)
            Xsec=Xsec.add_key_points2(E_head, E_end_i);
            Xsec=Xsec.add_key_points1(E_head, E_end_i);
        end
        % CCC����ֻ��1e3eV����ֻ��E_end����һ���㣬����total���Բ�ֵʱ����bug
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %��info�������ֵ��Ϣ
        
    case '2020LZS' % 2020����ɽ������CCC-m��2003Janev
        % ����ʹ��CCC (��1)
        % ����ʹ��CCC-m (��n=2, 3, 4)
        % ����ʹ��2003Janev (ѡȡn=1)
        %%% ׼��
        E_head=1e-3; % ����������Χ
        E_end=1e4;
        %%%%%%% ����CCC-m
        %         % �ݹ����
        %         Xsec=get_Xsec_e_H( 'CCC-m', path_list ,0);
        % Ϊ�˱����ظ�������Ϣ������CCC-m��������
        data_dir={[path_list.cross_sections '\e_H_CCC_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        Xsec=combine_exc_data(Xsec, 1, 6:9);
        Xsec=combine_exc_data(Xsec, 1, 3:5);
        Xsec=combine_exc_data(Xsec, 1, 1:2);
        Xsec.exc{1}{1}=[[0 0];[Xsec.excThresh{1}{1} 0];Xsec.exc{1}{1}];
        Xsec.exc{1}{2}=[[0 0];[Xsec.excThresh{1}{2} 0];Xsec.exc{1}{2}];
        Xsec.exc{1}{3}=[[0 0];[Xsec.excThresh{1}{3} 0];Xsec.exc{1}{3}];
        for E_end_i=logspace(log10(100),log10(E_end),15)
            Xsec=Xsec.add_key_points2(E_head, E_end_i);
            Xsec=Xsec.add_key_points1(E_head, E_end_i);
        end
        % CCC����ֻ��1e3eV����ֻ��E_end����һ���㣬����total���Բ�ֵʱ����bug
        %%%%%%% ���봦��
        % level 1
        Eth1=13.6; %��ֵ
        % ��ϱ��ʽϵ�� 2003Janev-tab. 3
        coeff1=[0.1845, -0.032226, -0.034539, 1.4003, -2.8115, 2.2986];
        com1=@(E, j) coeff1(j+1)*(1-Eth1/E)^j; % ��ϱ��ʽ����
        % ������ϱ��ʽ 2003Janev-ch 2.1.2-A-eq.14
        xsec1=@(E) 1e-17/Eth1/E...
            *(coeff1(0+1)*log(E/Eth1)+...
            com1(E,1)+com1(E,2)+com1(E,3)+com1(E,4)+com1(E,5));
        % ����Xsec.ion
        % log�����100������ ����
        x1=logspace(log10(Eth1),log10(E_end),100); % ��1000�����������
        y1=arrayfun(xsec1,x1);
        y1(1)=0; %��ֵ��Ӧ��������˸������������£�ֱ����Ϊ0
        Xsec.ion{1}{1}=[[0 0];[x1',y1']]; % Xsec.ion��energy��1��Ϊ0
        Xsec.ionThresh{1}{1}=Eth1;
        Xsec.info.ion{1}{1}='E+H(1S)->2E+H+(n=1),Ionization';
        %         % ��CCC����ion��energy�����ϸ���
        %         % CCCԭ����һ��ion��Ӧ������Ϊ0��13.6��14��1000eV��
        %         Xsec.ion{1}{1}(:,2)=arrayfun(xsec1,Xsec.ion{1}{1}(:,1));
        %         Xsec.ion{1}{1}(1, :)=[0 0]; % Xsec.ion��energy��1��Ϊ0
        %         Xsec.ion{1}{1}(2, :)=[Eth1 0]; % Xsec.ion��energy��2��Ϊ��ֵ
        %         Xsec.ionThresh{1}{1}=Eth1;
        %         Xsec.info.ion{1}{1}='E+H(1S)->2E+H+(n=1),Ionization';
        
        %%%%%%% ��������
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %��info�������ֵ��Ϣ
    case '1994Ness' % �ϲ�CCCͬ��ֵ������Ӧ���Լ�����ֵ���
        % ����ʹ��CCC (��1)
        % ����ʹ��CCC (��n=2s/p, 3s/p/d, 4s/p/f/d) �ϲ���n=2, 3, 4
        % ����ʹ��CCC (��n=1)
        %%%%%%% ����CCC
        data_dir={[path_list.root '\examples\example_ExB_drift']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
    otherwise
        error('[ERROR] No such type of Xsec_e_H')
end
% ���
disp(['[INFO] Get Xsec : ' Xsec.name{1} ' ok']);
if flag_display
    Xsec.output_info();
end
end

%% aid function
function compare_exc_data(Xsec, j, exc_pair)
% �Ƚ�exc i��exc j����ֵ���������ж��Ƿ���Ժϲ�
% exc_pair: 1*2 double�����Ƚϵ�����exc���
flag_Eth_equal=(roundn(Xsec.excThresh{j}{exc_pair(1)},-4)...
    ==roundn(Xsec.excThresh{j}{exc_pair(2)},-4));
flag_energy_equal=isequal(Xsec.exc{j}{exc_pair(1)}(:,1),Xsec.exc{j}{exc_pair(2)}(:,1));
if flag_Eth_equal && flag_energy_equal
    %     fprintf('%d��%d����ֵ��������ͬ���������ֱ�����\n',exc_pair(1),exc_pair(2))
else
    err_text=[num2str(exc_pair(1)) '��' num2str(exc_pair(2)) '����ֵ��������ͬ�����治��ֱ�����'];
    fprintf([err_text '\n'])
    throw(get_exception('error',err_text))
end
end

function Xsec=combine_exc_data(Xsec, j, k_used)
% combine exc in the array k_used into min(k_used) for the jth gas
num_exc=length(Xsec.exc{j});
k_used=sort(k_used);
assert(~isempty((k_used)))
assert(~(k_used(end)>num_exc)) %����k_used������������Ӧ��Ŀ
% �Ƚ���ֵ�ͽ��棬�ж��Ƿ���Ժϲ�
exc_pair=nchoosek(k_used,2);
size_pair=size(exc_pair);
try
    for k=1:size_pair(1)
        compare_exc_data(Xsec, j, exc_pair(k,:))
    end
catch e
    rethrow(e)
end
% ���ж���ֵ��ͬ����ֱ�������װ
% ��k_used�Ľ���ȫ����k_used(1)��
for i=k_used(2:end)
    Xsec.exc{1}{k_used(1)}(:,2)=Xsec.exc{1}{k_used(1)}(:,2)+Xsec.exc{1}{i}(:,2);
    Xsec.info.exc{1}{k_used(1)}=[Xsec.info.exc{1}{k_used(1)} ' & ' Xsec.info.exc{1}{i}];
end
%ɾ���ѽ�����͵ļ�����Ӧ
k_goal=sort(setdiff(1:num_exc,k_used(2:end)));
Xsec=Xsec.choose_exc(j, k_goal);
end

%% test aid function
% test �ֲ�������test�����ֶ�����
function test_compare_exc_data()
Xsec=ImportLXCat;
Xsec.excThresh{1}{1}=1; Xsec.exc{1}{1}=[2 3];
Xsec.excThresh{1}{2}=1; Xsec.exc{1}{2}=[2 3];
Xsec.excThresh{1}{3}=2; Xsec.exc{1}{3}=[2 3];
Xsec.excThresh{1}{4}=1; Xsec.exc{1}{4}=[1 3];
% test case 1&2: compare exc 1��exc 1/2 ��ͬ
try
    compare_exc_data(Xsec, 1, [1,1])
    compare_exc_data(Xsec, 1, [1,2])
catch e
    throw(get_exception('error','��Ӧ��catch������˲�Ӧ�����׳�'))
end
% test case 3: compare exc 1��exc 3 ��ֵ����ͬ
try
    compare_exc_data(Xsec, 1, [1,3])
catch e
    if e.identifier=='Nix_M:ERROR'
        output_exception(0, e);
    else
        throw(get_exception('error','catch�����쳣����Ԥ���쳣'))
    end
end
% test case 4: compare exc 1��exc 4 ��������ͬ
try
    compare_exc_data(Xsec, 1, [1,4])
catch e
    if e.identifier=='Nix_M:ERROR'
        output_exception(0, e);
    else
        throw(get_exception('error','catch�����쳣����Ԥ���쳣'))
    end
end
end

function test_combine_exc_data()
Xsec=ImportLXCat;
Xsec.excThresh{1}{1}=1; Xsec.exc{1}{1}=[2 3];
Xsec.info.exc{1}{1}='exc 1';
Xsec.excThresh{1}{2}=1; Xsec.exc{1}{2}=[2 4];
Xsec.info.exc{1}{2}='exc 2';
Xsec.excThresh{1}{3}=2; Xsec.exc{1}{3}=[3 4];
Xsec.info.exc{1}{3}='exc 3';
Xsec.excThresh{1}{4}=1; Xsec.exc{1}{4}=[5 6];
Xsec.info.exc{1}{4}='exc 4';
% test case 1&2: compare exc 1��exc 1/2 ��ͬ
e=get_exception('empty',''); %��չ�쳣�������÷�Χ
try
    Xsec=combine_exc_data(Xsec, 1, 1:2);
catch e
    rethrow(e)
end
assert(3==length(Xsec.exc{1}))
assert(3==length(Xsec.excThresh{1}))
assert(3==length(Xsec.info.exc{1}))
assert(isequal([2 7],Xsec.exc{1}{1}))
assert(isequal([3 4],Xsec.exc{1}{2}))
assert(2==Xsec.excThresh{1}{2})
Xsec.output_info();
% test case 3: compare exc 1��exc 3 ��ֵ����ͬ
try
    Xsec=combine_exc_data(Xsec, 1, 1:2);
catch e
    if e.identifier=='Nix_M:ERROR'
        output_exception(0, e);
    else
        throw(get_exception('error','catch�����쳣����Ԥ���쳣'))
    end
end
end