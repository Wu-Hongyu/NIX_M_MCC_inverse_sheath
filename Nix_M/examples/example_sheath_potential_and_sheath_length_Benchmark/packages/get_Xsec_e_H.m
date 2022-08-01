function [ Xsec ] = get_Xsec_e_H( type, flag_display)
% 处理ImportLXCat得到的e-H原始数据，获取MCC所需截面
% 能量范围 1e-2-1e3 eV内可靠 注意人为加密点以便于做线性插值
% 使用
% type：e_H散射截面数据种类名，一般即数据库名
% flag_display：0则不画图不输出反应种类，1则显示
% 注意，本文件根据.\others\choose_Xsec_e_n_of_hydrogen_plasma
% \compare_Xsec_e_H.m实现
% 物理性的比较在compare_Xsec_e_H
% 功能性的test在test_get_Xsec

% % 测试所使用局部函数，开发时手动使用
% test_compare_exc_data()
% test_combine_exc_data()
% % 20201217 these two test passed
path_list=get_path_init('test000');
recation_type={['e-H-' type]};  %处理成cell
switch type
    case 'CCC-m' % 合并CCC同阈值激发反应，以减少数值误差
        % 弹性使用CCC (有1)
        % 激发使用CCC (有n=2s/p, 3s/p/d, 4s/p/f/d) 合并成n=2, 3, 4
        % 电离使用CCC (有n=1)
        %%%%%%% 导入CCC
        data_dir={[path_list.cross_sections '\e_H_CCC_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        %%% 准备
        E_head=1e-3; % 可用能量范围
        E_end=1e4;
        
        %%%%%%% 激发处理
        % 从后往前combine，否则要计算新的序号
        Xsec=combine_exc_data(Xsec, 1, 6:9);
        Xsec=combine_exc_data(Xsec, 1, 3:5);
        Xsec=combine_exc_data(Xsec, 1, 1:2);
        % 由于存在振荡，后续难以自动补点，因此在此预先手动补点
        Xsec.exc{1}{1}=[[0 0];[Xsec.excThresh{1}{1} 0];Xsec.exc{1}{1}];
        Xsec.exc{1}{2}=[[0 0];[Xsec.excThresh{1}{2} 0];Xsec.exc{1}{2}];
        Xsec.exc{1}{3}=[[0 0];[Xsec.excThresh{1}{3} 0];Xsec.exc{1}{3}];
        
        %%%%%%% 其他处理
        for E_end_i=logspace(log10(100),log10(E_end),15)
            Xsec=Xsec.add_key_points2(E_head, E_end_i);
            Xsec=Xsec.add_key_points1(E_head, E_end_i);
        end
        % CCC数据只到1e3eV，若只在E_end处补一个点，后续total线性插值时可能bug
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %在info中添加阈值信息
        
    case '2020LZS' % 2020李增山，基于CCC-m和2003Janev
        % 弹性使用CCC (有1)
        % 激发使用CCC-m (有n=2, 3, 4)
        % 电离使用2003Janev (选取n=1)
        %%% 准备
        E_head=1e-3; % 可用能量范围
        E_end=1e4;
        %%%%%%% 导入CCC-m
        %         % 递归调用
        %         Xsec=get_Xsec_e_H( 'CCC-m', path_list ,0);
        % 为了避免重复附加信息，复制CCC-m操作如下
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
        % CCC数据只到1e3eV，若只在E_end处补一个点，后续total线性插值时可能bug
        %%%%%%% 电离处理
        % level 1
        Eth1=13.6; %阈值
        % 拟合表达式系数 2003Janev-tab. 3
        coeff1=[0.1845, -0.032226, -0.034539, 1.4003, -2.8115, 2.2986];
        com1=@(E, j) coeff1(j+1)*(1-Eth1/E)^j; % 拟合表达式因子
        % 截面拟合表达式 2003Janev-ch 2.1.2-A-eq.14
        xsec1=@(E) 1e-17/Eth1/E...
            *(coeff1(0+1)*log(E/Eth1)+...
            com1(E,1)+com1(E,2)+com1(E,3)+com1(E,4)+com1(E,5));
        % 更新Xsec.ion
        % log后均布100点能量 更新
        x1=logspace(log10(Eth1),log10(E_end),100); % 与1000点基本无区别
        y1=arrayfun(xsec1,x1);
        y1(1)=0; %阈值对应截面算出了负数，舍入误差导致，直接置为0
        Xsec.ion{1}{1}=[[0 0];[x1',y1']]; % Xsec.ion中energy第1个为0
        Xsec.ionThresh{1}{1}=Eth1;
        Xsec.info.ion{1}{1}='E+H(1S)->2E+H+(n=1),Ionization';
        %         % 在CCC现有ion的energy基础上更新
        %         % CCC原本有一个ion反应，能量为0，13.6，14到1000eV。
        %         Xsec.ion{1}{1}(:,2)=arrayfun(xsec1,Xsec.ion{1}{1}(:,1));
        %         Xsec.ion{1}{1}(1, :)=[0 0]; % Xsec.ion中energy第1个为0
        %         Xsec.ion{1}{1}(2, :)=[Eth1 0]; % Xsec.ion中energy第2个为阈值
        %         Xsec.ionThresh{1}{1}=Eth1;
        %         Xsec.info.ion{1}{1}='E+H(1S)->2E+H+(n=1),Ionization';
        
        %%%%%%% 其他处理
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %在info中添加阈值信息
    case '1994Ness' % 合并CCC同阈值激发反应，以减少数值误差
        % 弹性使用CCC (有1)
        % 激发使用CCC (有n=2s/p, 3s/p/d, 4s/p/f/d) 合并成n=2, 3, 4
        % 电离使用CCC (有n=1)
        %%%%%%% 导入CCC
        data_dir={[path_list.root '\examples\example_ExB_drift']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
    otherwise
        error('[ERROR] No such type of Xsec_e_H')
end
% 输出
disp(['[INFO] Get Xsec : ' Xsec.name{1} ' ok']);
if flag_display
    Xsec.output_info();
end
end

%% aid function
function compare_exc_data(Xsec, j, exc_pair)
% 比较exc i与exc j的阈值和能量，判断是否可以合并
% exc_pair: 1*2 double，待比较的两个exc序号
flag_Eth_equal=(roundn(Xsec.excThresh{j}{exc_pair(1)},-4)...
    ==roundn(Xsec.excThresh{j}{exc_pair(2)},-4));
flag_energy_equal=isequal(Xsec.exc{j}{exc_pair(1)}(:,1),Xsec.exc{j}{exc_pair(2)}(:,1));
if flag_Eth_equal && flag_energy_equal
    %     fprintf('%d和%d：阈值与能量相同，截面可以直接相加\n',exc_pair(1),exc_pair(2))
else
    err_text=[num2str(exc_pair(1)) '和' num2str(exc_pair(2)) '：阈值或能量不同，截面不能直接相加'];
    fprintf([err_text '\n'])
    throw(get_exception('error',err_text))
end
end

function Xsec=combine_exc_data(Xsec, j, k_used)
% combine exc in the array k_used into min(k_used) for the jth gas
num_exc=length(Xsec.exc{j});
k_used=sort(k_used);
assert(~isempty((k_used)))
assert(~(k_used(end)>num_exc)) %断言k_used不超出激发反应数目
% 比较阈值和截面，判断是否可以合并
exc_pair=nchoosek(k_used,2);
size_pair=size(exc_pair);
try
    for k=1:size_pair(1)
        compare_exc_data(Xsec, j, exc_pair(k,:))
    end
catch e
    rethrow(e)
end
% 已判断阈值相同，故直接相加组装
% 将k_used的截面全加在k_used(1)上
for i=k_used(2:end)
    Xsec.exc{1}{k_used(1)}(:,2)=Xsec.exc{1}{k_used(1)}(:,2)+Xsec.exc{1}{i}(:,2);
    Xsec.info.exc{1}{k_used(1)}=[Xsec.info.exc{1}{k_used(1)} ' & ' Xsec.info.exc{1}{i}];
end
%删除已进行求和的激发反应
k_goal=sort(setdiff(1:num_exc,k_used(2:end)));
Xsec=Xsec.choose_exc(j, k_goal);
end

%% test aid function
% test 局部函数的test，请手动运行
function test_compare_exc_data()
Xsec=ImportLXCat;
Xsec.excThresh{1}{1}=1; Xsec.exc{1}{1}=[2 3];
Xsec.excThresh{1}{2}=1; Xsec.exc{1}{2}=[2 3];
Xsec.excThresh{1}{3}=2; Xsec.exc{1}{3}=[2 3];
Xsec.excThresh{1}{4}=1; Xsec.exc{1}{4}=[1 3];
% test case 1&2: compare exc 1与exc 1/2 相同
try
    compare_exc_data(Xsec, 1, [1,1])
    compare_exc_data(Xsec, 1, [1,2])
catch e
    throw(get_exception('error','不应该catch到，因此不应该有抛出'))
end
% test case 3: compare exc 1与exc 3 阈值不相同
try
    compare_exc_data(Xsec, 1, [1,3])
catch e
    if e.identifier=='Nix_M:ERROR'
        output_exception(0, e);
    else
        throw(get_exception('error','catch到的异常不是预期异常'))
    end
end
% test case 4: compare exc 1与exc 4 能量不相同
try
    compare_exc_data(Xsec, 1, [1,4])
catch e
    if e.identifier=='Nix_M:ERROR'
        output_exception(0, e);
    else
        throw(get_exception('error','catch到的异常不是预期异常'))
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
% test case 1&2: compare exc 1与exc 1/2 相同
e=get_exception('empty',''); %扩展异常变量作用范围
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
% test case 3: compare exc 1与exc 3 阈值不相同
try
    Xsec=combine_exc_data(Xsec, 1, 1:2);
catch e
    if e.identifier=='Nix_M:ERROR'
        output_exception(0, e);
    else
        throw(get_exception('error','catch到的异常不是预期异常'))
    end
end
end