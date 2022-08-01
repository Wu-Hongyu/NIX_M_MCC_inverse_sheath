function write_file( filename, type, output )
% 写文件
% type: 参考 fopen-输入参数 permission - 文件访问类型
switch type
    case 'char2wt' % 写字符数组到文本文件
        out_file=fopen(filename,'wt');
        if out_file<0
            error(['[ERROR] can not write ' filename])
        end
        fprintf(out_file, '%s', output);
        % 注意，output中的\r\n是无法识别的
        fclose(out_file);
    otherwise
        error('[ERROR] not done')
end
end