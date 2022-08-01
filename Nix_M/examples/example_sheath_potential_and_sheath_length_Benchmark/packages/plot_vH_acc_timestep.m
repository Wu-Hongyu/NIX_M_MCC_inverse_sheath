function [ ] = plot_vH_acc_timestep( vH,k1,ti ,vH_acc_flag)
% 画图显示 新注入的离子的速度变化
for i=1:size(k1)
    if k1(i)==1
        vH_acc_flag(i)=1;
    end
end
for i=1:size(vH_acc_flag(i)==1)
    plot(vH(i),ti,'--r')
    hold on
end
xlabel('时间步长')
ylabel('离子速度')
end