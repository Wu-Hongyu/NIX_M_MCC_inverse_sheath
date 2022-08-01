function [new_y]=interp_linear_loglog(known_x,known_y,new_x)
% ��Ҫ��loglog����ϵ�������Ե��ڲ�����
% ���ܼ�test_plasma_physic_tools/test_interp_linear_loglog
% ����ֵ��ˮƽ�ߵĲ�ֵ���� ok
% ʹ��
% known_xӦ������ȥ��
% known_y���м䲻Ӧ����0

if 0==known_x(1)
    known_x(1)=1e-6;
end
new_y=interp1(log10(known_x),log10(known_y),...
    log10(new_x),'linear','extrap');
new_y=10.^new_y;

% ������ֵ����0�������ڲ�ֵ���õ�0.
% �������ĩ��0���0���紦���������Բ�ֵ
idx_p1=find(known_y>0,1);
if idx_p1>1 && 0==known_y(idx_p1-1)
        idx=find((new_x>known_x(idx_p1-1))&new_x<known_x(idx_p1));
        if ~isempty(idx)
            new_y(idx)=interp1(known_x,known_y,new_x(idx),'linear','extrap');
        end
end

idx_p2=find(known_y>0,1,'last');
if idx_p2<length(known_y) && 0==known_y(idx_p2+1)
        idx=find((new_x>known_x(idx_p2))&new_x<known_x(idx_p2+1));
        if ~isempty(idx)
            new_y(idx)=interp1(known_x,known_y,new_x(idx),'linear','extrap');
        end
end

assert(~any(isinf(new_y)))
assert(~any(isnan(new_y)))

end