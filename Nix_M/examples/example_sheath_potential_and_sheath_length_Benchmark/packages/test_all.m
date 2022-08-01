clc
clear
testfile_name_list={...
    'test_get_v_init.m',...
    'test_get_position_init.m',...
    'test_get_field_init.m',...
    'test_get_u.m',...
    'test_solve_possion_equation.m',... %get_field_init.m和get_u两个文件组成的package
    'test_get_v_and_position.m',... %各种pusher组成的package
    'test_get_log_init.m',...
    };
num_testfile=length(testfile_name_list);

log_name = get_log_init( testfile_name_list', 'test' );
for i=1:num_testfile
    results = runtests(testfile_name_list{i});
    diary(log_name) % 重定向方式输出日志
    disp(table(results));
    diary off
    % pause
end

% results = runtests('test_get_v_init.m');
% disp(table(results));
% 
% results = runtests('test_get_position_init.m');
% disp(table(results));
% % pause
% results = runtests('test_get_field_init.m');
% disp(table(results));
% % pause
% results = runtests('test_get_u.m');
% disp(table(results));
% % pause
% results = runtests('test_solve_possion_equation.m');
% %solve_possion_equation是get_field_init.m和get_u两个文件组成的package
% disp(table(results));
% % pause
% results = runtests('test_get_v_and_position.m');
% disp(table(results));
% % pause
% results = runtests('test_get_log_init.m');
% disp(table(results));
% % pause