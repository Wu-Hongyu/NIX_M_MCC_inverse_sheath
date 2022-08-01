%% Main function to generate tests
function tests = test_importLXCat
% test importLXCat
tests = functiontests(localfunctions);
end

%% Test Functions

% TODO
% 对新加的一些成员函数做简单test

% TODO：test原有函数
% 根据e-H2-CCC的exc GTSC与importLXcat得到的tot相等，
% importLXcat求总截面功能 ok

function test_combine_exc(testCase)
% test combine_exc

end

% TODO：修改以下test，一个两种气体的test放最后，基本不用

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
% % Xsec.ion 1×1 cell，元素为1×1 cell
% % Xsec.ion{1} cell，元素为[222×2 double]
% % Xsec.ion{1}{1} 222×2 double，即能量(eV)与截面 与IONIZATION数据和数目一致
% % Xsec.ionThresh{1}{1}为一个double，即Xsec.ion{1}{1}的阈值
% % Xsec.exc{1} 1×44 cell 数组，各为能量(eV)与截面 与 EXCITATION数目一致
% % Xsec.ionThresh{1} 1×44 cell 数组，各为激发阈值
% % Xsec.energy 11156×1 double，所有反应种类的全部能量值，去重、排序
% % Xsec.tot{1} 11156×1 double，
% verifyEqual(testCase,class(Xsec),'ImportLXCat')
% verifyEqual(testCase,length(Xsec.ion),1)
% verifyEqual(testCase,size(Xsec.ion{1}{1}),[222,2])
% verifyEqual(testCase,Xsec.ionThresh{1}{1},15.76)
% verifyEqual(testCase,size(Xsec.exc{1}),[1,44])
% verifyEqual(testCase,size(Xsec.tot{1}),[11156,1])
% 
% % % % unittest中人工判断
% % Xsec=Xsec.init(gasName, gasDir, 1);
% % answer = questdlg('截面绘图符合预期？','人工判断','Y','N','Y');
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
% % Xsec.ion 1×2 cell 数组 
% % Xsec.ion{1} cell，元素为[222×2 double]
% % Xsec.ion{1}{1} 222×2 double，即能量(eV)与截面
% % 与e:\GitRepos\CPIP_code\_Xsection\Ar_Biagi\xsections.txt中IONIZATION数据一致
% % Xsec.exc{1} 1×44 cell 数组，各为能量(eV)与截面 与 EXCITATION数据数目一致
% verifyEqual(testCase,Xsec.name,gasName)
% verifyEqual(testCase,length(Xsec.ion),2)
% verifyEqual(testCase,size(Xsec.ion{1}{1}),[222,2])
% verifyEqual(testCase,Xsec.ionThresh{1}{1},15.76)
% verifyEqual(testCase,size(Xsec.exc{1}),[1,44])
% verifyEqual(testCase,size(Xsec.tot{1}),[18400,1])
% 
% % % unittest中人工判断
% % Xsec=Xsec.init(gasName, gasDir, 1);
% % answer = questdlg('截面绘图符合预期？','人工判断','Y','N','Y');
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
