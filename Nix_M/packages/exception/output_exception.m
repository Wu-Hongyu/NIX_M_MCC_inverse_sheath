% TODO
% ?? exception(I)
function [ output ] = output_exception(output_type, exception)
% ����쳣
num=length(exception(:));
% �ṹ����������Ԥ����
for i=1:num
    output(i).out_type=['exception ' exception.identifier];
    output(i).message=exception.message;
    output(i).stack=struct2table(exception.stack);
end

switch output_type
    case 0 %���������Ϣ������̨
        for i=1:num
            disp(output(i).message)
            disp(output(i).stack)
            fprintf('\n')
        end
    case 1 %����쳣����������̨�����ڵ���
        disp(output(i).out_type)
            disp(output(i).message)
            disp(output(i).stack)
            fprintf('\n')
    otherwise
        if 'file'==class(output_type)
            % TODO: ������ļ�
            error('not done')
        else
            e=get_exception('not_supported','this out_type is not supported by output_exception()');
            output_exception(0, e)
            throw(e)
        end
end
end