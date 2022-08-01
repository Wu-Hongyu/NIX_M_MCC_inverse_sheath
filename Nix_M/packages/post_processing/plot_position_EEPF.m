function [ ] = plot_position_EEPF( ve,constants)
%PLOT_POSITION_EEPF 此处显示有关此函数的摘要
%   此处显示详细说明
if ~isempty(ve)
    EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
    EV=EineV(ve,constants.q_m_ratio_e);
%     EDGES_vec=[0 0.5 1 1.5 2 3 4 5 6 7 8 9 10 12 14 16 18 20 23 26 29];
    EDGES_vec=[0 0.1 0.2 0.4 0.6 0.8 1 1.5 2 2.5 3 3.5 4 5 6 7 8 9 10 12 14 16];
    [N,EDGES]=histcounts(EV,EDGES_vec,'Normalization','pdf');
%     [N,EDGES]=histcounts(EV,'Normalization','pdf');
    if length(EDGES)>15
        if all(N(1:end))
            %  h=histogram(EV,1000,'Normalization','pdf');
            % figure
            % semilogy(EDGES(2:end),(N(1:end)./sqrt((EDGES(1:end-1)+EDGES(2:end))/2)));
            plot(EDGES(2:end),log10(N(1:end)./sqrt((EDGES(1:end-1)+EDGES(2:end))/2)));
            grid on
            hold on
            
            xxx=EDGES(2:end);
            yy=(EDGES(1:(end-1))+EDGES(2:end))/2;
            yyy=log10(N(1:end)./sqrt(yy(1:end)));
            % TE=9/(1.5);
            % EEPF=2/sqrt(pi)*TE^(-3/2)*exp(-xxx/TE);
            % myfunc = inline('log(2/sqrt(pi)*TE^(-3/2)*exp(-xxx/TE))','TE','xxx');
            EEPF_func=@(TE,xxx) (log10(2/sqrt(pi)*TE^(-3/2)*exp(-xxx/TE)));
            TE0 = 5;%待定系数的预估值
            opts = statset('nlinfit');
            opts.FunValCheck = 'off';
            TE_cal = nlinfit(xxx,yyy,EEPF_func,TE0,opts);%
%             TE_cal = nlinfit(xxx,yyy,EEPF_func,TE0);%
            EEPF_cal=EEPF_func(TE_cal,EDGES);
            plot(EDGES,EEPF_cal)
            % semilogy(xxx,log(EEPF))
            % semilogy(h.BinEdges(2:end),h.BinCounts(1:end)./sqrt(h.BinEdges(2:end)));
            text2function(EDGES,EEPF_cal,[num2str(TE_cal) 'eV'],12)
        end
%          if all(N(1:end/2))
%             %  h=histogram(EV,1000,'Normalization','pdf');
%             % figure
%             % semilogy(EDGES(2:end),(N(1:end)./sqrt((EDGES(1:end-1)+EDGES(2:end))/2)));
%             plot(EDGES(2:end),log(N(1:end)./sqrt((EDGES(1:end-1)+EDGES(2:end))/2)));
%             grid on
%             hold on
%             
%             xxx=EDGES(2:(end+1)/2);
%             yy=(EDGES(1:(end-1))+EDGES(2:end))/2;
%             yyy=log(N(1:end/2)./sqrt(yy(1:end/2)));
%             % TE=9/(1.5);
%             % EEPF=2/sqrt(pi)*TE^(-3/2)*exp(-xxx/TE);
%             % myfunc = inline('log(2/sqrt(pi)*TE^(-3/2)*exp(-xxx/TE))','TE','xxx');
%             EEPF_func=@(TE,xxx) (log(2/sqrt(pi)*TE^(-3/2)*exp(-xxx/TE)));
%             TE0 = 5;%待定系数的预估值
%             TE_cal = nlinfit(xxx,yyy,EEPF_func,TE0);%
%             EEPF_cal=EEPF_func(TE_cal,EDGES);
%             plot(EDGES,EEPF_cal)
%             % semilogy(xxx,log(EEPF))
%             % semilogy(h.BinEdges(2:end),h.BinCounts(1:end)./sqrt(h.BinEdges(2:end)));
%             text2function(EDGES,EEPF_cal,[num2str(TE_cal) 'eV'],12)
%         end
    end
    
end
end

function text2function(x,y,string,textsize)
% write text "string" on the maximum of function y(x)

Y = (max(y)+min(y))/2;
X = (max(x)+min(x))/2;

text(X,Y, string ,'fontsize',textsize)

end
