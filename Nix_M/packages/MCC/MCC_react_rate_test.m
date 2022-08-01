%% Main function to generate tests
function tests = MCC_react_rate_test
% test filename
tests = functiontests(localfunctions);
end

%% Test Functions
function test_reactionNumNorm2rateCoefficient(testCase)
% test
constants=get_constants();% 全局常数 结构体
simulation=get_simulation('default'); %仿真参数 结构体
simulation.num0_macro_e=200000;%初始电子数目 Particle e Num
EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
% Te=[ 10:10:50];
Te=logspace(log10(0.1),log10(50),50);
% Te=50;
for i=1:length(Te)
    simulation.Te0=Te(i);
    veth=sqrt(-constants.q_m_ration_e*simulation.Te0);%电子温度对应的热速度
    % vHth=sqrt(constants.q_m_ration_H*simulation.TH0);%离子温度对应的热速度
    [mcc_para,Xsec]=MCC_init(simulation,constants,0.3);
    ve=get_v_init( veth, 'Maxwellian velocity', [simulation.num0_macro_e,3]);
    % v_value=sqrt(ve(:,1).^2+ve(:,2).^2+ve(:,3).^2);
    % Te2v=sqrt(constants.kB*1/constants.me);
    for j=1:mcc_para.ion_index
        Eev_list=Xsec.ion{1,1}{1,j}(:,1);
        CS_list=Xsec.ion{1,1}{1,j}(:,2);
        reaction_rate(j,i)=trapz(Eev_list , 2/sqrt(pi)*sqrt(Eev_list)*Te(i)^(-3/2).*exp(-Eev_list/Te(i)).*sqrt(2*Eev_list*constants.e/constants.me).*CS_list);
        % vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
    end
    
    %暂时没写att的
    
    for j=mcc_para.att_index+1:mcc_para.ela_index
        Eev_list=Xsec.ela{1,1}{1,1}(:,1);
        CS_list=Xsec.ela{1,1}{1,1}(:,2);
        reaction_rate(j,i)=trapz(Eev_list , 2/sqrt(pi)*sqrt(Eev_list)*Te(i)^(-3/2).*exp(-Eev_list/Te(i)).*sqrt(2*Eev_list*constants.e/constants.me).*CS_list);
        % vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
    end
    
    for j=mcc_para.ela_index+1:mcc_para.exc_index
        Eev_list=Xsec.exc{1,1}{1,j-mcc_para.ela_index}(:,1);
        CS_list=Xsec.exc{1,1}{1,j-mcc_para.ela_index}(:,2);
        reaction_rate(j,i)=trapz(Eev_list , 2/sqrt(pi)*sqrt(Eev_list)*Te(i)^(-3/2).*exp(-Eev_list/Te(i)).*sqrt(2*Eev_list*constants.e/constants.me).*CS_list);
        % vH=get_v_init( vHth, 'Maxwellian velocity', [simulation.num0_macro_H,3]);
    end
    
    count_colli_type=0;
    count_colli_energy=0;
    
    fprintf('当前时步数：%d/%d \n',i,length(Te))
    for ti=1:20000
        
        %------粒子MCC处理------------------------------------------------------------------------------------------
        [collision_index,collision_type]=null_collision( ve,constants,mcc_para,Xsec );
        
        %      [colli_result]=collision_process(Xsec, mcc_para,ve,constants,collision_index,collision_type );
        %     if colli_result.happen==1
        %         ve(colli_result.ion_index,:)=colli_result.ion_inc_v;
        %         ve(colli_result.ela_index,:)=colli_result.ela_v;
        %         ve(colli_result.exc_index,:)=colli_result.exc_v;
        %         ve(colli_result.att_index,:)=[];
        %         ve=[ve; colli_result.ion_gen_v ];
        %         position_e=[position_e;position_e(colli_result.ion_index)];
        %         vH2=[vH2; normrnd(0,vHth,[length(colli_result.ion_gen_v) , 3]);  ];
        %     end
        
        %------粒子MCC处理结束------------------------------------------------------------------------------------------
        
        count_colli_type=[count_colli_type ; collision_type];
        E=EineV(ve(collision_index,:),constants.q_m_ration_e);
        %     count_colli_energy=[count_colli_energy; E];
        
    end
    %     exc_V1_energy=count_colli_energy(count_colli_type==6);
    %     ela_energy=count_colli_energy(count_colli_type==3);
    for reac_num=1:mcc_para.exc_index
        reaction_count(reac_num,i)=length(count_colli_type(count_colli_type==reac_num));
        %     ela_colli_count(i)=length(ela_energy);
    end
end

path_list=get_path_init('test000');
temp_result_file=[path_list.output '\MCC_rate_testResults.mat'];
save(temp_result_file);
% load(temp_result_file)

figure
for i=1:length(reaction_rate(:,1))
    loglog(Te,reaction_count(i,:)/max(max(reaction_count)),'-.r')
    hold on
    loglog(Te,reaction_rate(i,:)/max(max(reaction_rate)),'--k')
    if i<=mcc_para.ion_index
        text2function(Te,reaction_rate(i,:)/max(max(reaction_rate)),'ion',12)
    elseif i<=mcc_para.att_index
    elseif i<=mcc_para.ela_index
        text2function(Te,reaction_rate(i,:)/max(max(reaction_rate)),'ela',12)
    else
        text2function(Te,reaction_rate(i,:)/max(max(reaction_rate)),['exc',num2str(Xsec.excThresh{1}{i-5}),'eV'],12)
    end
    pause
end
axis([-inf,inf,1E-7,1E1])
grid on
xlabel('Te/eV') ; ylabel('normalized reaction rate');
title('MCC reaction rate compared to calculated rate')
hold off
legend('MCC react count(normalized)','react rate calculated from Maxwellian(normalized)')
end

%% aid function
function text2function(x,y,string,textsize)
% write text "string" on the maximum of function y(x)

y_max = max(y);
y_max = y_max(1);
x_max = x(y==max(y));
x_max = x_max(1);
text(x_max,y_max, string ,'fontsize',textsize)

end

%% Optional file fixtures
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end
