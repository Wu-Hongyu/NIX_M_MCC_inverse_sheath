%% Main function to generate tests
function tests = test_importLXCat
% test importLXCat
tests = functiontests(localfunctions);
end

%% Test Functions

% TODO
% ���¼ӵ�һЩ��Ա��������test

% TODO��testԭ�к���
% ����e-H2-CCC��exc GTSC��importLXcat�õ���tot��ȣ�
% importLXcat���ܽ��湦�� ok

function test_combine_exc(testCase)
% test combine_exc

end

% TODO���޸�����test��һ�����������test����󣬻�������

% function test_one_gas(testCase)
% % test one type of gas
% path_list=get_path_init('test000');
% % sumforumla gases
% gasName = {'Ar'};
% % directories of cross sections
% gasDir = {
%     [path_list.cross_sections '\Ar_Biagi'],...
%     };
% % plot data (1) or not(0)
% interactive  = 0;
% 
% Xsec = ImportLXCat;
% Xsec=Xsec.init(gasName, gasDir, interactive);
% 
% disp(Xsec)
% Xsec.output_info()
% % Xsec.ion 1��1 cell��Ԫ��Ϊ1��1 cell
% % Xsec.ion{1} cell��Ԫ��Ϊ[222��2 double]
% % Xsec.ion{1}{1} 222��2 double��������(eV)����� ��IONIZATION���ݺ���Ŀһ��
% % Xsec.ionThresh{1}{1}Ϊһ��double����Xsec.ion{1}{1}����ֵ
% % Xsec.exc{1} 1��44 cell ���飬��Ϊ����(eV)����� �� EXCITATION��Ŀһ��
% % Xsec.ionThresh{1} 1��44 cell ���飬��Ϊ������ֵ
% % Xsec.energy 11156��1 double�����з�Ӧ�����ȫ������ֵ��ȥ�ء�����
% % Xsec.tot{1} 11156��1 double��
% verifyEqual(testCase,class(Xsec),'ImportLXCat')
% verifyEqual(testCase,length(Xsec.ion),1)
% verifyEqual(testCase,size(Xsec.ion{1}{1}),[222,2])
% verifyEqual(testCase,Xsec.ionThresh{1}{1},15.76)
% verifyEqual(testCase,size(Xsec.exc{1}),[1,44])
% verifyEqual(testCase,size(Xsec.tot{1}),[11156,1])
% 
% % % % unittest���˹��ж�
% % Xsec=Xsec.init(gasName, gasDir, 1);
% % answer = questdlg('�����ͼ����Ԥ�ڣ�','�˹��ж�','Y','N','Y');
% % verifyEqual(testCase,answer,'Y')
% end
% 
% function test_two_gas(testCase)
% % test one type of gas
% path_list=get_path_init('test000');
% % sumforumla gases
% gasName = {'Ar','N2'};
% % directories of cross sections
% gasDir = {
%     [path_list.cross_sections '\Ar_Biagi'],...
%     [path_list.cross_sections '\N2_Biagi'],...
%     };
% % plot data (1) or not(0)
% interactive  = 0;
% 
% Xsec = ImportLXCat;
% Xsec=Xsec.init(gasName, gasDir, interactive);
% 
% disp(Xsec)
% Xsec.output_info()
% % Xsec.ion 1��2 cell ���� 
% % Xsec.ion{1} cell��Ԫ��Ϊ[222��2 double]
% % Xsec.ion{1}{1} 222��2 double��������(eV)�����
% % ��e:\GitRepos\CPIP_code\_Xsection\Ar_Biagi\xsections.txt��IONIZATION����һ��
% % Xsec.exc{1} 1��44 cell ���飬��Ϊ����(eV)����� �� EXCITATION������Ŀһ��
% verifyEqual(testCase,Xsec.name,gasName)
% verifyEqual(testCase,length(Xsec.ion),2)
% verifyEqual(testCase,size(Xsec.ion{1}{1}),[222,2])
% verifyEqual(testCase,Xsec.ionThresh{1}{1},15.76)
% verifyEqual(testCase,size(Xsec.exc{1}),[1,44])
% verifyEqual(testCase,size(Xsec.tot{1}),[18400,1])
% 
% % % unittest���˹��ж�
% % Xsec=Xsec.init(gasName, gasDir, 1);
% % answer = questdlg('�����ͼ����Ԥ�ڣ�','�˹��ж�','Y','N','Y');
% % verifyEqual(testCase,answer,'Y')
% end


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
