function [new_y]=interp_linear_loglog(known_x,known_y,new_x)
% 主要在loglog坐标系中做线性的内插和外插
% 性能见test_plasma_physic_tools/test_interp_linear_loglog
% 对零值、水平线的插值能力 ok
% 使用
% known_x应已排序，去重
% known_y在中间不应该有0

if 0==known_x(1)
    known_x(1)=1e-6;
end
new_y=interp1(log10(known_x),log10(known_y),...
    log10(new_x),'linear','extrap');
new_y=10.^new_y;

% 上述插值，在0的邻域内插值均得到0.
% 因此在首末的0与非0交界处，改用线性插值
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