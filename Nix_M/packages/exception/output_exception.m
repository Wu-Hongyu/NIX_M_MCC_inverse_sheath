% TODO
% ?? exception(I)
function [ output ] = output_exception(output_type, exception)
% 输出异常
num=length(exception(:));
% 结构体数组难以预分配
for i=1:num
    output(i).out_type=['exception ' exception.identifier];
    output(i).message=exception.message;
    output(i).stack=struct2table(exception.stack);
end

switch output_type
    case 0 %输出报错信息到控制台
        for i=1:num
            disp(output(i).message)
            disp(output(i).stack)
            fprintf('\n')
        end
    case 1 %输出异常变量到控制台，用于调试
        disp(output(i).out_type)
            disp(output(i).message)
            disp(output(i).stack)
            fprintf('\n')
    otherwise
        if 'file'==class(output_type)
            % TODO: 输出到文件
            error('not done')
        else
            e=get_exception('not_supported','this out_type is not supported by output_exception()');
            output_exception(0, e)
            throw(e)
        end
end
end