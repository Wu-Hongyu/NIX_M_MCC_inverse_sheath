%% Main function to generate tests
function tests = test_get_Xsec
% test get Xsec
tests = functiontests(localfunctions);
end

%% Test Functions
function test_use_get_Xsec(testCase)
% how to use different get_Xsec_reactantPair
% 主要是对get与使用各种反应物组合的Xsec做示例
%%%%%%%%%%%%%%%%%%%% get Xsec
% Xsec_reactantPair=get_Xsec_reactantPair(type, path_list ,flag_display)
% flag_display=1则自动画图和输出信息
Xsec_e_H=get_Xsec_e_H( '2020LZS', 0);

%%%%%%%%%%%%%%%%%%%% 使用Xsec
Xsec_e_H.plotXsections(14,'log','log'); %截面画图
Xsec_e_H.output_info() % 输出截面信息
%%%%%%%%%%%%%%% 对象Xsec的属性
Xsec_e_H.name{1} % char数组，反应物组合名+数据源类型名
%%%%%%%% exc
length(Xsec_e_H.exc{1}) %double, 若大于1，则为exc反应数目
k=2;
Xsec_e_H.info.exc{1}{k} % char数组，反应过程等信息
Xsec_e_H.excThresh{1}{k} % double，exc能量阈值，eV
Xsec_e_H.exc{1}{k} % 2*N double，第k个exc反应的截面数据
Xsec_e_H.exc{1}{k}(:,1) % 1*N double，exc截面对应能量，eV
Xsec_e_H.exc{1}{k}(:,2) % 1*N double，exc截面，m^2
%%%%%%%% ion
% 同exc
%%%%%%%% att
length(Xsec_e_H.att{1}) % 等于1，可能为att反应数目，也可能为
% 提示信息占位
k=1;
Xsec_e_H.info.att{1}{k} % 显示 no Attachment process
Xsec_e_H.attThresh{1}{k} %阈值为空, isempty
Xsec_e_H.att{1}{k} % 截面全为0
%%%%%%%% tot
Xsec_e_H.energy %1*M double，tot截面对应能量，eV
% 覆盖所有反应能量范围
Xsec_e_H.tot{1} %1*M double，tot截面，m^2
end

function test_check_get_Xsec(testCase)
% test use different get_Xsec_reactantPair
%%%%%%%%%%%%%%%%%%%% get Xsec e-n
Xsec1_e_H=get_Xsec_e_H( '2020LZS', 1);
Xsec2_e_H=get_Xsec_e_H( 'CCC-m', 1);
Xsec1_e_H2=get_Xsec_e_H2( '2020LZS' ,1);
Xsec2_e_H2=get_Xsec_e_H2( 'Phelps-m' ,1);
end


function test_combine_exc_data(testCase)
% 直接导入的CCC与combine_exc_data合并处理后的CCC对比total
flag_display=0; % 如想人工查看，请置1；否则置0
path_list=get_path_init('test000');
data_dir={[path_list.cross_sections '\e_H_CCC_20201214']};
Xsec = ImportLXCat;
Xsec=Xsec.init({'e-H-CCC_direct_imported'}, data_dir, flag_display);
Xsec=Xsec.add_threshold_to_info();
if flag_display
    Xsec.output_info();
end
Xsec_m=Xsec;
Xsec_m=combine_exc_data(Xsec_m, 1, 6:9);
Xsec_m=combine_exc_data(Xsec_m, 1, 3:5);
Xsec_m=combine_exc_data(Xsec_m, 1, 1:2);
Xsec_m=Xsec_m.totalXsection();
verifyEqual(testCase,Xsec_m.tot{1},Xsec.tot{1},'RelTol',5e-2); %总截面相同
if flag_display
    figure
    loglog(Xsec.energy,Xsec.tot{1},'-r','LineWidth',5)
    hold on
    loglog(Xsec_m.energy,Xsec_m.tot{1},'-.k')
    legend('before combine','after combine','location','best')
    grid on
end
end
%PengChen: 人工查看：二者一致 2020/12/17 20:32:43


function test_choose_exc_and_get_total(testCase)
% 直接导入的CCC与get_Xsec_e_H2处理后的CCC对比
% 在choose_exc和重新计算总截面后，应接近CCC的GTSC
flag_display=1; % 如想人工查看，请置1；否则置0
path_list=get_path_init('test000');
data_dir={[path_list.cross_sections '\e_H2_CCC_20201214']};
Xsec_CCC = ImportLXCat;
Xsec_CCC=Xsec_CCC.init({'e-H2-CCC_direct_imported'}, data_dir, flag_display);
Xsec_CCC=Xsec_CCC.add_threshold_to_info();
if flag_display
    Xsec_CCC.output_info();
end
Xsec_LZS=get_Xsec_e_H2( '2020LZS', flag_display);
% verifyEqual(testCase,Xsec_m.tot{1},Xsec.tot{1},'RelTol',5e-2); %总截面相同
idx1=Xsec_LZS.energy<300;
idx1(1)=false; % CCC-ela前侧水平补点，LZS-ela前侧loglog补点，因此忽略第一个
xsec1=[Xsec_LZS.energy(idx1) Xsec_LZS.tot{1}(idx1)];
idx2=Xsec_CCC.exc{1}{1}(:,1)<300;
idx2(1)=false;
xsec2=Xsec_CCC.exc{1}{1}(idx2,:);
verify_xsec_equal_roughly(testCase,xsec1,xsec2,0.05)
if flag_display
    figure
    legend_text={};
    legend_text{end+1}='GTSC CCC';
    loglog(Xsec_CCC.exc{1}{1}(:,1),Xsec_CCC.exc{1}{1}(:,2),'-g','LineWidth',5)
    hold on
    legend_text{end+1}='exc1 CCC-modified';
    loglog(Xsec_LZS.exc{1}{1}(:,1),Xsec_LZS.exc{1}{1}(:,2),'-r')
    legend_text{end+1}='total CCC';
    loglog(Xsec_CCC.energy,Xsec_CCC.tot{1},'--g','LineWidth',5)
    legend_text{end+1}='total CCC-modified';
    loglog(Xsec_LZS.energy,Xsec_LZS.tot{1},'--r')
    L1=legend(legend_text);
    set(L1,'location','best');
    grid on
end
end
%PengChen: 人工查看：ok 
% CCC-m保留所有反应的total与CCC的GTSC相同，
% 说明ImportLXCat的total功能ok，对CCC的理解ok

% .\others\choose_Xsec_e_n_of_hydrogen_plasma\Nix_M截面与文献对比.pptx
% 表明Nix_M使用公式与2003Janev绘图一致
% 因此以下弃用
% function test_e_H_2020LZS(testCase)
% % 比较'2020LZS' 和LZS的mat文件，验证'2020LZS'中ion公式
% flag_display=0; % 如想人工查看，请置1；否则置0
% path_list=get_path_init('test000');
% Xsec=get_Xsec_e_H( '2020LZS', flag_display);
% 
% mat_name=[path_list.others ...
%     '\choose_Xsec_e_n_of_hydrogen_plasma' ...
%     '\e_H_LiZengshan_20201214\ionCroSec1.mat'];
% in_struct=load(mat_name);
% in_cell=struct2cell(in_struct);
% xsec_LZS=cell2mat(in_cell);
% 
% % 通过插值到相同点的结果，自动比较不同离散曲线的重合度
% verify_xsec_equal_roughly(testCase,Xsec.ion{1}{1},xsec_LZS,0.05)
% 
% if flag_display
%     figure
%     loglog(Xsec.ion{1}{1}(:,1),Xsec.ion{1}{1}(:,2),'-r','LineWidth',5)
%     hold on
%     loglog(xsec_LZS(:,1), xsec_LZS(:,2), '-.k','LineWidth',4)
%     legend('NIX-LZS','origin-LZS','location','best')
%     grid on
% end
% end
% %PengChen: 人工查看：二者一致 2020/12/18 16:52:40

%% Optional file fixtures
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end

%% aid funciton
function verify_xsec_equal_roughly(testCase,xsec1,xsec2,rel_tol_lg)
% test中用于判断截面是否大致相同
% 在二者能量重叠段，线性插值100点，lg(截面)相对误差<rel_tol_lg
x1=xsec1(:,1);
y1=xsec1(:,2);
x2=xsec2(:,1);
y2=xsec2(:,2);
xmin=max([min(x1) min(x2)]);
xmax=min([max(x1) max(x2)]);
x=xmin:(xmax-xmin)/100:xmax;
y1_interp=interp1(x1,y1,x);
y2_interp=interp1(x2,y2,x);
verifyEqual(testCase, y1_interp, y2_interp, 'RelTol',rel_tol_lg);
verifyEqual(testCase, log10(y1_interp), log10(y2_interp), 'RelTol',rel_tol_lg);
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
