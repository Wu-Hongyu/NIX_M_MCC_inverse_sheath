%% Main function to generate tests
function tests = test_get_Xsec
% test get Xsec
tests = functiontests(localfunctions);
end

%% Test Functions
function test_use_get_Xsec(testCase)
% how to use different get_Xsec_reactantPair
% ��Ҫ�Ƕ�get��ʹ�ø��ַ�Ӧ����ϵ�Xsec��ʾ��
%%%%%%%%%%%%%%%%%%%% get Xsec
% Xsec_reactantPair=get_Xsec_reactantPair(type, path_list ,flag_display)
% flag_display=1���Զ���ͼ�������Ϣ
Xsec_e_H=get_Xsec_e_H( '2020LZS', 0);

%%%%%%%%%%%%%%%%%%%% ʹ��Xsec
Xsec_e_H.plotXsections(14,'log','log'); %���滭ͼ
Xsec_e_H.output_info() % ���������Ϣ
%%%%%%%%%%%%%%% ����Xsec������
Xsec_e_H.name{1} % char���飬��Ӧ�������+����Դ������
%%%%%%%% exc
length(Xsec_e_H.exc{1}) %double, ������1����Ϊexc��Ӧ��Ŀ
k=2;
Xsec_e_H.info.exc{1}{k} % char���飬��Ӧ���̵���Ϣ
Xsec_e_H.excThresh{1}{k} % double��exc������ֵ��eV
Xsec_e_H.exc{1}{k} % 2*N double����k��exc��Ӧ�Ľ�������
Xsec_e_H.exc{1}{k}(:,1) % 1*N double��exc�����Ӧ������eV
Xsec_e_H.exc{1}{k}(:,2) % 1*N double��exc���棬m^2
%%%%%%%% ion
% ͬexc
%%%%%%%% att
length(Xsec_e_H.att{1}) % ����1������Ϊatt��Ӧ��Ŀ��Ҳ����Ϊ
% ��ʾ��Ϣռλ
k=1;
Xsec_e_H.info.att{1}{k} % ��ʾ no Attachment process
Xsec_e_H.attThresh{1}{k} %��ֵΪ��, isempty
Xsec_e_H.att{1}{k} % ����ȫΪ0
%%%%%%%% tot
Xsec_e_H.energy %1*M double��tot�����Ӧ������eV
% �������з�Ӧ������Χ
Xsec_e_H.tot{1} %1*M double��tot���棬m^2
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
% ֱ�ӵ����CCC��combine_exc_data�ϲ�������CCC�Ա�total
flag_display=0; % �����˹��鿴������1��������0
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
verifyEqual(testCase,Xsec_m.tot{1},Xsec.tot{1},'RelTol',5e-2); %�ܽ�����ͬ
if flag_display
    figure
    loglog(Xsec.energy,Xsec.tot{1},'-r','LineWidth',5)
    hold on
    loglog(Xsec_m.energy,Xsec_m.tot{1},'-.k')
    legend('before combine','after combine','location','best')
    grid on
end
end
%PengChen: �˹��鿴������һ�� 2020/12/17 20:32:43


function test_choose_exc_and_get_total(testCase)
% ֱ�ӵ����CCC��get_Xsec_e_H2������CCC�Ա�
% ��choose_exc�����¼����ܽ����Ӧ�ӽ�CCC��GTSC
flag_display=1; % �����˹��鿴������1��������0
path_list=get_path_init('test000');
data_dir={[path_list.cross_sections '\e_H2_CCC_20201214']};
Xsec_CCC = ImportLXCat;
Xsec_CCC=Xsec_CCC.init({'e-H2-CCC_direct_imported'}, data_dir, flag_display);
Xsec_CCC=Xsec_CCC.add_threshold_to_info();
if flag_display
    Xsec_CCC.output_info();
end
Xsec_LZS=get_Xsec_e_H2( '2020LZS', flag_display);
% verifyEqual(testCase,Xsec_m.tot{1},Xsec.tot{1},'RelTol',5e-2); %�ܽ�����ͬ
idx1=Xsec_LZS.energy<300;
idx1(1)=false; % CCC-elaǰ��ˮƽ���㣬LZS-elaǰ��loglog���㣬��˺��Ե�һ��
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
%PengChen: �˹��鿴��ok 
% CCC-m�������з�Ӧ��total��CCC��GTSC��ͬ��
% ˵��ImportLXCat��total����ok����CCC�����ok

% .\others\choose_Xsec_e_n_of_hydrogen_plasma\Nix_M���������׶Ա�.pptx
% ����Nix_Mʹ�ù�ʽ��2003Janev��ͼһ��
% �����������
% function test_e_H_2020LZS(testCase)
% % �Ƚ�'2020LZS' ��LZS��mat�ļ�����֤'2020LZS'��ion��ʽ
% flag_display=0; % �����˹��鿴������1��������0
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
% % ͨ����ֵ����ͬ��Ľ�����Զ��Ƚϲ�ͬ��ɢ���ߵ��غ϶�
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
% %PengChen: �˹��鿴������һ�� 2020/12/18 16:52:40

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
% test�������жϽ����Ƿ������ͬ
% �ڶ��������ص��Σ����Բ�ֵ100�㣬lg(����)������<rel_tol_lg
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
