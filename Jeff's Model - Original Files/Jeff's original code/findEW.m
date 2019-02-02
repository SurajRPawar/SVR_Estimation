function [EW PVA] = findEW(V,P,V2,P2,LV_par,T_per)




E_es_lv = LV_par(1)
A_lv = LV_par(2)
B_lv = LV_par(3)
V0_lv = LV_par(4)
T_lv = LV_par(5)
shift = LV_par(6)



EW = 0; 
PVA = 0;




% Integrate for EW
for i = 2:length(V)+1
    if i == length(V)+1
        if V(1) < V(i-1)
            EW = EW + abs(V(1)-V(i-1)) * (P(1)+P(i-1))/2;
        else
            EW = EW - abs(V(1)-V(i-1)) * (P(1)+P(i-1))/2;
        end
    else
        if V(i) < V(i-1)
            EW = EW + abs(V(i)-V(i-1)) * (P(i)+P(i-1))/2;
        else
            EW = EW - abs(V(i)-V(i-1)) * (P(i)+P(i-1))/2;
        end
    end
end

imax = round(length(V)*(shift*T_lv)/T_per);
V_ED = 5:(V(1)-V0_lv)/500:V(1);
V_ES = 5:(V(imax)-V0_lv)/100:V(imax);

P_ED = A_lv * (exp(B_lv * (V_ED-V0_lv))-1);
P_ES = (V_ES-V0_lv)*(E_es_lv);

%Integrate for PVA
for i = 2:imax
    if V(i) < V(i-1)
        PVA = PVA + abs(V(i)-V(i-1)) * (P(i)+P(i-1))/2;
    else
        PVA = PVA - abs(V(i)-V(i-1)) * (P(i)+P(i-1))/2;
    end
end
for i = 2:length(V_ED)
    PVA = PVA - abs(V_ED(i)-V_ED(i-1)) * (P_ED(i)+P_ED(i-1))/2;
end
for i = 2:length(V_ES)
    PVA = PVA + abs(V_ES(i)-V_ES(i-1)) * (P_ES(i)+P_ES(i-1))/2;
end

PVA = PVA;% * 133.28 * 10^-6; %convert ml-mmHg to J or W-s
EW  = EW ;% * 133.28 * 10^-6;


Vp = 5:1:500;
P_EDp = A_lv * (exp(B_lv * (Vp-V0_lv))-1);
P_ESp = (Vp-V0_lv)*(E_es_lv);






% fig8 = figure(8);
% plot(Vp,P_EDp,'--','Color',[0.4,0.4,0.4],'LineWidth',1)
% hold on
% plot(Vp,P_ESp,'--','Color',[0.4,0.4,0.4],'LineWidth',1)
% % plot(V,P,'Color',[.4 .4 .4],'LineWidth',3)
% % plot(V,P,'Color',[.75 0 0],'LineWidth',2)
% % plot(V,P,'Color',[0 .75 0],'LineWidth',2)
% % plot(V,P,'Color',[0 0 1],'LineWidth',2)
% plot(V,P,'Color',[0 0 0],'LineWidth',2)
% % set(gca,'XTick',[0:50:200]);
% % set(gca,'XTickLabel',{'0';'';'100';'';'200'});
% % set(gca,'YTick',[0:25:125]);
% % set(gca,'YTickLabel',{'0';'';'50';'';'100';''});
% xlim([0 310])
% ylim([0 130])
% ylabel('LVP (mmHg)','FontSize',9)
% xlabel('LVV (ml)','FontSize',9)
% set(gca,'FontSize',8)
% set(fig8,'Units','inches') 
% set(fig8,'Position',[2.5 4.5 3.25 3])
% 
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);


