function [ collision_result,ion_react_order ] = collision_process( Xsec,mcc_para,ve,constants,collision_index,collision_type ,react)
%COLLISION_PROCESS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
ion_react_order=[];
switch react
case 'H2'
if ~isempty(collision_index)
    ion_index=collision_index(collision_type<=mcc_para.ion_index );%��Ӧ���Ͷ�Ӧ������index
    ion_react_order=collision_type(collision_type<=mcc_para.ion_index );%���뷴Ӧ�ڵ��뷴Ӧ����Ķ�Ӧ���
    att_index=collision_index(collision_type<=mcc_para.att_index & collision_type>mcc_para.ion_index);
    ela_index=collision_index(collision_type<=mcc_para.ela_index & collision_type>mcc_para.att_index);
    exc_index=collision_index(collision_type<=mcc_para.exc_index & collision_type>mcc_para.ela_index);
    exc_react_order=collision_type(collision_type<=mcc_para.exc_index & collision_type>mcc_para.ela_index)-mcc_para.ela_index;%��Ӧ�ļ�����Ӧ�ڼ�����Ӧ����ı��
    
    [ ion_inc_post_v ,ion_gen_post_v] = ion_process( ion_index , ion_react_order, ve , constants ,  Xsec.ionThresh );%�����ɢ�䴦��
    
    collision_result.att_index=att_index;
    
    [  ela_post_v ]=ela_process( ela_index,ve,constants,react);%������ײ��ɢ�䴦��
    
    [ exc_post_v ] = exc_process( exc_index,exc_react_order,ve,constants  ,  Xsec.excThresh );%������ɢ�䴦��
    
    
    collision_result.ion_index=ion_index;
    collision_result.ion_inc_v=ion_inc_post_v;
    collision_result.ion_gen_v=ion_gen_post_v;
    
    collision_result.ela_index=ela_index;
    collision_result.ela_v=ela_post_v;
    
    collision_result.exc_index=exc_index;
    collision_result.exc_v=exc_post_v;
    
    collision_result.happen=1;
else
    collision_result.happen=0;
    
end

case 'H'
if ~isempty(collision_index)
    ion_index=collision_index(collision_type<=mcc_para.ion_index_H );%��Ӧ���Ͷ�Ӧ������index
    ion_react_order=collision_type(collision_type<=mcc_para.ion_index_H );%���뷴Ӧ�ڵ��뷴Ӧ����Ķ�Ӧ���
    att_index=collision_index(collision_type<=mcc_para.att_index_H & collision_type>mcc_para.ion_index_H);
    ela_index=collision_index(collision_type<=mcc_para.ela_index_H & collision_type>mcc_para.att_index_H);
    exc_index=collision_index(collision_type<=mcc_para.exc_index_H & collision_type>mcc_para.ela_index_H);
    exc_react_order=collision_type(collision_type<=mcc_para.exc_index_H & collision_type>mcc_para.ela_index_H)-mcc_para.ela_index_H;%��Ӧ�ļ�����Ӧ�ڼ�����Ӧ����ı��
    
    [ ion_inc_post_v ,ion_gen_post_v] = ion_process( ion_index , ion_react_order, ve , constants ,  Xsec.ionThresh );%�����ɢ�䴦��
    
    collision_result.att_index=att_index;
    
    [  ela_post_v ]=ela_process( ela_index,ve,constants,react);%������ײ��ɢ�䴦��
    
    [ exc_post_v ] = exc_process( exc_index,exc_react_order,ve,constants  ,  Xsec.excThresh );%������ɢ�䴦��
    
    
    collision_result.ion_index=ion_index;
    collision_result.ion_inc_v=ion_inc_post_v;
    collision_result.ion_gen_v=ion_gen_post_v;
    
    collision_result.ela_index=ela_index;
    collision_result.ela_v=ela_post_v;
    
    collision_result.exc_index=exc_index;
    collision_result.exc_v=exc_post_v;
    
    collision_result.happen=1;
else
    collision_result.happen=0;
    
end


end

