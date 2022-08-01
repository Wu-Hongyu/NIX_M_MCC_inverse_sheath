function [ ] = plot_vH_grid( vH_grid, constants ,Te)
plot(vH_grid,'-b','LineWidth',3)%������������ٶ�ƽ��ֵ
        xlabel('������');
        ylabel('�����ٶ�[m/s]');
        axis([0,size(vH_grid,1),-inf,inf])
        uB=sqrt(constants.e*Te/constants.mH);
        line([0,size(vH_grid,1)],[uB,uB],'linestyle','-.','color','r','LineWidth',3)
        drawnow;
end