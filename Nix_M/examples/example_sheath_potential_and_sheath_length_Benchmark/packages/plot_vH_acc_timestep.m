function [ ] = plot_vH_acc_timestep( vH,k1,ti ,vH_acc_flag)
% ��ͼ��ʾ ��ע������ӵ��ٶȱ仯
for i=1:size(k1)
    if k1(i)==1
        vH_acc_flag(i)=1;
    end
end
for i=1:size(vH_acc_flag(i)==1)
    plot(vH(i),ti,'--r')
    hold on
end
xlabel('ʱ�䲽��')
ylabel('�����ٶ�')
end