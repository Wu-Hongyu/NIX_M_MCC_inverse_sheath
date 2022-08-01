clc
clear
addpath(genpath('../packages')) % ���./packages�������ļ��е�·��
results = runtests('test_get_path_init');%�Ȳ���·������
result_tables=table(results);

% ����Nix_M�����д�����test_ǰ׺���ļ�����ӵ�list
path_list=get_path_init('now');
cd(path_list.packages)
list=ls('**/test_*.m');
cd(path_list.root)
testfile_name_list=cellstr(list);
idx=zeros(2,1);
idx(1)=find(strcmp(testfile_name_list, 'test_all.m' ));
idx(2)=find(strcmp(testfile_name_list,  'test_get_path_init.m'));
testfile_name_list(idx)=[];
num_testfile=length(testfile_name_list);
%����
for i=1:num_testfile
    results = runtests(testfile_name_list{i}); %�ϲ����table����
    result_tables=[result_tables;table(results)];
end

%     results = runtests('MCC_react_rate_test'); %�ϲ����table����
%     result_tables=[result_tables;table(results)];

log_name=[path_list.output '\test_all.log'];
clc
diary(log_name) % �ض���ʽ�����־
disp(result_tables);
idx=find(~result_tables.Passed);
if isempty(idx)
    disp('All unit test passed.')
else
    num_failed=length(idx);
    for i=1:num_failed
        warning([result_tables.Name{idx(i)} ' failed!'])
        record=result_tables.Details{idx(i)}.DiagnosticRecord;
        disp(record.Report)
    end
end
diary off
disp('See the details in Nix_M/other/project_name/test_all.log.')