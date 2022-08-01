function [ simulation ] = get_simulation( type )
% simulation ������� �ṹ��
switch type
    case 'default' %A1DPIC_Ver1.0.m, 2019Montellano
        simulation.dt=10*10^-12;%ʱ�䲽��
        simulation.end_time=6E-6;%����ʱ���ܳ�
        simulation.all_timesteps=floor(simulation.end_time/simulation.dt);%��ѭ������
        simulation.Lx=0.01;%�������򳤶�
        simulation.source_region=[0.1*simulation.Lx,0.4*simulation.Lx];
        simulation.num_grid_point=201;%��������Ŀ
        simulation.dx=simulation.Lx/(simulation.num_grid_point-1);%�ռ䲽��
        simulation.n0=1E16;%���������ܶ�
        simulation.num0_macro_e=20000;%��ʼ������Ŀ Particle e Num
        simulation.num0_macro_H=simulation.num0_macro_e;%��ʼ������ĿParticle H Num
        simulation.weight=simulation.n0*simulation.Lx/simulation.num0_macro_e;%ÿ�������Ӵ���ʵ�����ӵ�Ȩ��
        simulation.Te=1;%��λeV
        simulation.TH=1;%��λeV
        simulation.field_boundaries_type=0;%��λ�߽���������
        simulation.field_boundaries=[0,0];%��λ�߽�����ֵ
        simulation.pressure=0.6;%��ѹ ��λPa
        simulation.Tgas=300;%����   K
    case 'BACON base'
        error('Not Done');
end

