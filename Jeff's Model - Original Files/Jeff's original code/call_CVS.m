clear;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
% These parameters can be adjusted to simulate various patient conditions
% and VAD settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heart failure (0 - healthy, 1 - severe failure, 2 - average failure, 3 - recovery)
    failure = 1;
    
% Type of assistance (0 - no VAD, 1 - HMII, 2- TORVAD)
    flow_type = 2;
    
% Asycnrhonous setting (2 - 9krpm, 3 - 9.6krpm.  9-9.6 krpm is the typical clinical
% range of the pump)
    asynch_setting = 2;
    
% TORVAD stroke volume (ml)  The M3T5 has a 30 ml stroke volume
    SV = 30;
    
% TORVAD phasing - normally delayed at .24 seconds for counter pulse
% (should do a phasing study on this)
% Tvad(1) is the R-wave delay for the first stroke
% Tvad(2) is the time when VAD ejection ends
    Tvad = [.24 .24+.2 0 0];     
    
% Open (1) or closed (0) loop system
% when open, the preload and afterload can be set to a constant
    open = 0;


%vary phase delay
% parameter = 0:.1:.6;
parameter = 0;

parameter1 = 30;
parameter2 = 0;

x0 = 0;

for j = 1:length(parameter1)
    for k = 1:length(parameter2)
    
%assign value of parameter to phase delay
% Tvad(1) = parameter(j);
% Tvad(2) = parameter(j)+.3;




if open == 1
    pre = parameter(j);
%     pre = 18;
%     after = parameter(j);
%     after = 0.8;
    after = 1.1;
%     after = 1.2;
    Tvad = [.24 .24+.3 .7 .7+.3];
else
    pre = 0; after = 0; 
end

% PD = mod(parameter2(k)*60/90-.15,60/90)

% Tvad = [PD PD+.3 0 0];
SV = parameter1(j);


% set an array called parameters to all of the input parameters for the
% cardiovascular system model.
parameters = [failure...
              flow_type...
              asynch_setting...
              SV...
              Tvad(1)...
              Tvad(2)...
              Tvad(3)...
              Tvad(4)...
              open...
              pre...
              after];

          
[data1 data2 LV_par] = CVS(parameters);


t(j,:) = data1(:,1);
V_la(j,:) = data1(:,2);  %left atrial volume
V_lv(j,:) = data1(:,3);  %left ventricular volume
P_sa(j,:) = data1(:,4);  %systemic arterial pressure
Q_sa(j,:) = data1(:,5);
P_st(j,:) = data1(:,6);  %arterial tree pressure 
P_sv(j,:) = data1(:,7);  %right atrial pressure
V_ra(j,:) = data1(:,8);  %right ventricular pressure
V_rv(j,:) = data1(:,9);  
P_pa(j,:) = data1(:,10);
Q_pa(j,:) = data1(:,11);
P_pt(j,:) = data1(:,12);
P_pv(j,:) = data1(:,13);
e_t_la(j,:) = data1(:,14);
e_t_lv(j,:) = data1(:,15);
e_t_ra(j,:) = data1(:,16);
e_t_rv(j,:) = data1(:,17);
P_la(j,:) = data1(:,18);
P_lv(j,:) = data1(:,19);
Q_vad(j,:) = data1(:,20);
Q_m(j,:) = data1(:,21);
Q_a(j,:) = data1(:,22);
Q_st(j,:)= data1(:,23);
Q_sv(j,:)= data1(:,24);     
P_ra(j,:)= data1(:,25);
P_rv(j,:)= data1(:,26);
Q_t(j,:)= data1(:,27);
Q_p(j,:)= data1(:,28);
Q_pt(j,:)= data1(:,29);
Q_pv(j,:)= data1(:,30);
Q_cor(j,:)= data1(:,31);

output_data2(:,j,k) = data2';
% output_data1(:,j,k) = data1';
% output_LV_par(:,j,k) = LV_par';

%     output_data2(1) = failure;
%     output_data2(2) = flow_type;
%     if flow_type ~= 1
%         output_data2(3) = 0;
%     else
%         output_data2(3) = asynch_setting;
%     end
%     if flow_type ~= 2
%         output_data2(4) = 0;
%         output_data2(5) = 0;
%         output_data2(6) = 0;
%     else
%         output_data2(4) = SV;
%         output_data2(5) = Tvad(1);
%         output_data2(6) = Tvad(2)-Tvad(1);   
%     end
%     output_data2(7) = HR;
%     output_data2(8) = CO*60/1000; 
%     output_data2(9) = max(P_sa(i1:i2));  
%     output_data2(10) = min(P_sa(i1:i2));  
%     output_data2(11) = map;  
%     output_data2(12) = max(P_sa(i1:i2))-min(P_sa(i1:i2));
%     output_data2(13) = mean(Q_a(i1:i2))*60/1000;
%     output_data2(14) = Q_vad_avg*60/1000;
%     output_data2(15) = Qcor_avg;
%     output_data2(16) = SHE;
%     output_data2(17) = EEP;
%     output_data2(18) = VO2;
%     output_data2(19) = LVEDV;
%     output_data2(20) = LVESV;
%     output_data2(21) = (LVEDV - LVESV)/LVEDV*100;
%     output_data2(22) = LVEDP;
%     output_data2(23) = mean(P_la(i1:i2));
%     output_data2(24) = mean(P_ra(i1:i2));
%     output_data2(25) = max(P_pa(i1:i2));
%     output_data2(26) = min(P_pa(i1:i2));
%     output_data2(27) = mean(P_pa(i1:i2));
%     output_data2(28) = PVA;
%     output_data2(29) = max(P_lv(i1:i2));
%     output_data2(30) = min(P_lv(i1:i2));
%     output_data2(31) = mean(P_lv(i1:i2));

    end
end

%%

for i = 1:length(parameter2)
    CO(i) = output_data2(8,1,i)
end


figure(6)
plot(parameter2-.15,CO)


[m,n] = size(t);
n = n-1;

RR = max(t(j,1:n+1));
tc = cat(2,t(j,1:n),t(j,1:n)+RR,t(j,1:n)+2*RR,t(j,1:n)+3*RR,t(j,1:n)+4*RR);
Plv = cat(2,P_lv(j,1:n),P_lv(j,1:n),P_lv(j,1:n),P_lv(j,1:n),P_lv(j,1:n));
Psa = cat(2,P_sa(j,1:n),P_sa(j,1:n),P_sa(j,1:n),P_sa(j,1:n),P_sa(j,1:n));
Pla = cat(2,P_la(j,1:n),P_la(j,1:n),P_la(j,1:n),P_la(j,1:n),P_la(j,1:n));
Qa = cat(2,Q_a(j,1:n),Q_a(j,1:n),Q_a(j,1:n),Q_a(j,1:n),Q_a(j,1:n));
Qvad = cat(2,Q_vad(j,1:n),Q_vad(j,1:n),Q_vad(j,1:n),Q_vad(j,1:n),Q_vad(j,1:n));
Vlv = cat(2,V_lv(j,1:n),V_lv(j,1:n),V_lv(j,1:n),V_lv(j,1:n),V_lv(j,1:n));
Vla = cat(2,V_la(j,1:n),V_la(j,1:n),V_la(j,1:n),V_la(j,1:n),V_la(j,1:n));



%%

fig1 =figure(4);
subplot(2,1,1)
set(gca,'FontSize',8,'LineWidth',2)
plot(tc,Psa,'black-',tc,Plv,'black--',tc,Pla,'black:','LineWidth',2)
% legend('AOP','LVP','LAP','Location','NorthEastOutside','Fontsize',8)
hold on
plot([.1 .3],[115 115],'-','Color',[0 0 0],'LineWidth',2)
text(.3,115,'AOP','FontSize',8)
plot([.7 .9],[115 115],'--','Color',[0 0 0],'LineWidth',2)
text(.9,115,'LVP','FontSize',8)
plot([1.3 1.5],[115 115],':','Color',[0 0 0],'LineWidth',2)
text(1.5,115,'LAP','FontSize',8)
hold off
% legend boxoff
xlim([0 2.25])
ylim([0 130])
title('Model LV failure','FontSize',9)
% title('Model CF support','FontSize',9)
xlabel('time (sec)','FontSize',8)
% title('Model TORVAD support','FontSize',10)
% title('Model Baseline','FontSize',10)



subplot(2,1,2)
% subaxis(2,2,4,'PaddingRight', -.03)
set(gca,'FontSize',8)
% plot(t(j,:),Q_a(j,:),'black-',t(j,:),Q_vad(j,:),'black--')
% plot(tc,Qa,'black-')
plot(tc,Qa,'black-',tc,Qvad,'black--','LineWidth',2)
hold on
plot([.1 .3],[300 300],'-','Color',[0 0 0],'LineWidth',2)
text(.3,300,'QAO','FontSize',8)
plot([.7 .9],[300 300],'--','Color',[0 0 0],'LineWidth',2)
text(.9,300,'QVAD','FontSize',8)
hold off
% legend('QAO','Location','SouthEast')
% legend('QAO','QVAD','Location','NorthEastOutside','Fontsize',8)
% legend boxoff
xlim([0 2.25])
ylim([-75 750])
xlabel('time (sec)','FontSize',8)







[t_ecg,ecg] = complete();
t_ecg = t_ecg*.8;
t_off = .55;

fig7 = figure(7);
subplot(4,1,[1 2])
set(gca,'FontSize',8,'LineWidth',2)
plot(tc-t_off,Psa,'black-',tc-t_off,Plv,'black--','LineWidth',2)
ylabel('')
xlim([0 2])
box off
set(gca,'XTickLabel','','YTickLabel','')
set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);




subplot(4,1,3)
set(gca,'FontSize',8)
plot(t_ecg-t_off+.02,ecg,'black','LineWidth',2)
ylabel('')
xlim([0 2])
box off
set(gca,'XTickLabel','','YTickLabel','')
set(gca,'FontSize',8)
set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);




subplot(4,1,4)
set(gca,'FontSize',8)
plot(tc-t_off,Qvad,'black','LineWidth',2)
ylim([min(Qvad-10) max(Qvad+10)])
xlim([0 2])
box off
ylabel('')
set(gca,'XTickLabel','','YTickLabel','')
set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);


set(fig7,'Units','inches') 
set(fig7,'Position',[2.5 4.5 6 5])



% subplot(3,1,3)
% plot(tc,Vla,tc,Vlv)
% xlim([0 2.25])

for i = 1:length(data2)
    fprintf('%f\n',data2(i));
end






figure(5)
subplot(2,1,1)
plot(tc,Qvad,':')
xlim([0 1])
subplot(2,1,2)
plot(tc,Psa-Plv,':')
xlim([0 1])
hold on




