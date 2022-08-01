function [ collision_result,ion_react_order ] = collision_process( Xsec,mcc_para,ve,constants,collision_index,collision_type ,react)
%COLLISION_PROCESS 此处显示有关此函数的摘要
%   此处显示详细说明
ion_react_order=[];
switch react
case 'H2'
if ~isempty(collision_index)
    ion_index=collision_index(collision_type<=mcc_para.ion_index );%反应类型对应的粒子index
    ion_react_order=collision_type(collision_type<=mcc_para.ion_index );%电离反应在电离反应组里的对应编号
    att_index=collision_index(collision_type<=mcc_para.att_index & collision_type>mcc_para.ion_index);
    ela_index=collision_index(collision_type<=mcc_para.ela_index & collision_type>mcc_para.att_index);
    exc_index=collision_index(collision_type<=mcc_para.exc_index & collision_type>mcc_para.ela_index);
    exc_react_order=collision_type(collision_type<=mcc_para.exc_index & collision_type>mcc_para.ela_index)-mcc_para.ela_index;%对应的激发反应在激发反应组里的编号
    
    [ ion_inc_post_v ,ion_gen_post_v] = ion_process( ion_index , ion_react_order, ve , constants ,  Xsec.ionThresh );%电离的散射处理
    
    collision_result.att_index=att_index;
    
    [  ela_post_v ]=ela_process( ela_index,ve,constants,react);%弹性碰撞的散射处理
    
    [ exc_post_v ] = exc_process( exc_index,exc_react_order,ve,constants  ,  Xsec.excThresh );%激发的散射处理
    
    
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
    ion_index=collision_index(collision_type<=mcc_para.ion_index_H );%反应类型对应的粒子index
    ion_react_order=collision_type(collision_type<=mcc_para.ion_index_H );%电离反应在电离反应组里的对应编号
    att_index=collision_index(collision_type<=mcc_para.att_index_H & collision_type>mcc_para.ion_index_H);
    ela_index=collision_index(collision_type<=mcc_para.ela_index_H & collision_type>mcc_para.att_index_H);
    exc_index=collision_index(collision_type<=mcc_para.exc_index_H & collision_type>mcc_para.ela_index_H);
    exc_react_order=collision_type(collision_type<=mcc_para.exc_index_H & collision_type>mcc_para.ela_index_H)-mcc_para.ela_index_H;%对应的激发反应在激发反应组里的编号
    
    [ ion_inc_post_v ,ion_gen_post_v] = ion_process( ion_index , ion_react_order, ve , constants ,  Xsec.ionThresh );%电离的散射处理
    
    collision_result.att_index=att_index;
    
    [  ela_post_v ]=ela_process( ela_index,ve,constants,react);%弹性碰撞的散射处理
    
    [ exc_post_v ] = exc_process( exc_index,exc_react_order,ve,constants  ,  Xsec.excThresh );%激发的散射处理
    
    
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

