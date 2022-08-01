function [ input ] = read_file( filename, type )
% 读文件
% type: 参考 fopen-输入参数 permission - 文件访问类型
switch type
    case 'rt2char' % 从文本文件读字符数组
        in_file=fopen(filename,'rt');
        if in_file<0
            error(['[ERROR] can not read ' filename])
        end
        line_no=0;
        input='';
        while ~feof(in_file) % 逐行读入
            line = fgetl(in_file);
            line_no=line_no+1;
            input=[input line];
        end
        fclose(in_file);
    otherwise
        error('[ERROR] not done')
end
end