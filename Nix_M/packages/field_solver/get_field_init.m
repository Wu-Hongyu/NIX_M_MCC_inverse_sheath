function [ u, A, b_extra ] = get_field_init( simulation )
n=simulation.num_grid_point;
% 处理电位边界条件，组装系数矩阵
u=zeros(n,1);%电位初值
leading_diagonal=-2*ones(n,1);
lower_diagonal=ones(n,1);
upper_diagonal=ones(n,1);
% spdiags组装n × n 的稀疏系数矩阵
% 位于主对角线下方的对角线首先从列的顶部获取元素。
% 位于主对角线上方的对角线首先从列的底部获取元素。
b_extra=zeros(2,1);%对右端向量的额外操作
switch simulation.field_boundaries_type
    case 0 %00，首末第一类边界
        u(1)=simulation.field_boundaries(1);
        leading_diagonal(1)=1;
        upper_diagonal(2)=0;
        b_extra(1)=u(1);
        
        u(end)=simulation.field_boundaries(2);
        leading_diagonal(end)=1;
        lower_diagonal(end-1)=0;
        b_extra(end)=u(end);
    case 1 %01，首第一类，末第二类
        u(1)=simulation.field_boundaries(1);
        leading_diagonal(1)=1;
        upper_diagonal(2)=0;
        b_extra(1)=u(1);
        
        lower_diagonal(end-1)=2;
        b_extra(end)=-2*simulation.field_boundaries(2)*simulation.dx;
    case 2 %10，首第二类，末第一类
        upper_diagonal(2)=2;
        b_extra(1)=2*simulation.field_boundaries(1)*simulation.dx;
        
        u(end)=simulation.field_boundaries(2);
        leading_diagonal(end)=1;
        lower_diagonal(end-1)=0;
        b_extra(end)=u(end);
    case 3 %11，周期边界条件
        disp('使用电场周期边界条件，边界电位取为零参考电位，忽略simulation.field_boundaries')
        u(1)=0;
        leading_diagonal(1)=1;
        upper_diagonal(2)=0;
        b_extra(1)=u(1);
        
        u(end)=0;
        leading_diagonal(end)=1;
        lower_diagonal(end-1)=0;
        b_extra(end)=u(end);
    otherwise
        error('Not Done');
end
A=spdiags([lower_diagonal,leading_diagonal,upper_diagonal],[-1,0,1],n,n);
% 减小右端向量计算量
constants=get_constants();% 全局常数 结构体
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
A=coeff_temp*A;
b_extra=coeff_temp*b_extra;
end
