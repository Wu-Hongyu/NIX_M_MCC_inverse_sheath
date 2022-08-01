function [ input ] = read_file( filename, type )
% ���ļ�
% type: �ο� fopen-������� permission - �ļ���������
switch type
    case 'rt2char' % ���ı��ļ����ַ�����
        in_file=fopen(filename,'rt');
        if in_file<0
            error(['[ERROR] can not read ' filename])
        end
        line_no=0;
        input='';
        while ~feof(in_file) % ���ж���
            line = fgetl(in_file);
            line_no=line_no+1;
            input=[input line];
        end
        fclose(in_file);
    otherwise
        error('[ERROR] not done')
end
end