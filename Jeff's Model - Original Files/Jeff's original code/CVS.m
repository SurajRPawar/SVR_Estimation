function [output_data1 output_data2 LV_par] = CVS(parameters)


% Input parameters
failure         = parameters(1);
flow_type       = parameters(2);
asynch_setting  = parameters(3);
SV              = parameters(4);
Tvad(1)         = parameters(5);
Tvad(2)         = parameters(6);
Tvad(3)         = parameters(7);
Tvad(4)         = parameters(8);
open            = parameters(9);
pre             = parameters(10);
after           = parameters(11);

% Set initial conditions and HR based on healthy or failure
if failure == 0
    HR = 75;
    x0 = [71.0 124.4 76.6 59.1 64.8 4.9 37.9 107.2 15.1 45.3 11.8 9.7 0 0 0 0 0];
elseif failure == 1
    HR = 90;
    x0 = [93.3 198.1 67.1 16.1 64.7 20.8 87.5 130.5 36.9 39.4 29.2 28.5 0 0 0 0 0];
    if flow_type ~= 0
        x0 = [93.3 98.1 67.1 16.1 64.7 20.8 87.5 32.4 36.9 39.4 29.2 28.5 0 0 0 0 0];
    end
elseif failure == 2
    HR = 90;
    x0 = [93.6 234.6694   58.2305   25.2608   54.5053    2.3311   17.4917   74.4681   34.1537   51.6578   30.8598   28.9859 0 0 0 0 0];
    if flow_type ~= 0
        x0 = [93.6496  234.6694-00   58.2305   25.2608   54.5053    2.3311   17.4917   74.4681   34.1537   51.6578   30.8598   28.9859 0 0 0 0 0];
    elseif flow_type == 1
        x0 = [75.5001  220.5770   77.1101   84.0615   64.4526   15.5583   78.1569  125.5997   18.3786   70.7016   14.8956   11.7362 0 0 0 0 0];
    elseif flow_type == 2
        x0 = [78.3812  231.5137   64.3567   31.3159   59.7964   13.6982   73.5115  124.0444   19.6268   66.7168   16.3929   13.4734 0 0 0 0 0];
    end
elseif failure == 3
    HR = 90;
    x0 = [93.6 234.6694   58.2305   25.2608   54.5053    2.3311   17.4917   74.4681   34.1537   51.6578   30.8598   28.9859 0 0 0 0 0];
    if flow_type ~= 0
        x0 = [93.6496  234.6694-00   58.2305   25.2608   54.5053    2.3311   17.4917   74.4681   34.1537   51.6578   30.8598   28.9859 0 0 0 0 0];
    elseif flow_type == 1
        x0 = [75.5001  220.5770   77.1101   84.0615   64.4526   15.5583   78.1569  125.5997   18.3786   70.7016   14.8956   11.7362 0 0 0 0 0];
    elseif flow_type == 2
        x0 = [78.3812  231.5137   64.3567   31.3159   59.7964   13.6982   73.5115  124.0444   19.6268   66.7168   16.3929   13.4734 0 0 0 0 0];
    end
end



RR = 60/HR;
tspan = [0 RR]; %ODE timespan

% number of points in output array (per heartbeat)
num = 500;
xp = ones(num+1,12);

%solve differential equations for a single cardiac cycle until solution
%stabilizes (for synchronous mode where each cycle is the same)
cycles = 1;
fprintf('Cycle\t Convergence(%%)\n')
while cycles < 45

    % solve ODE for one cardiac cycle (0->RR) w/ initial conditions x0
    [tode,xode] = ode23(@(t,x) CVS_dot_nv(t,x,parameters,HR),tspan,x0);

    %interpolate ODE results
    t = 0:RR/num:RR;
    x = interp1(tode,xode,t);                                                                                                           

    % caculate percent error from the last run
    if cycles == 1
        error = 100;
        errorp(cycles,:) = ones(1,12)*100;
    else
        error = 100 * max(abs(x(:,1:12)-xp)) ./ mean(xp);
        errorp(cycles,:) = error;
    end
    xp = x(:,1:12);

    % print error to screen to track convergence
    fprintf('%2.0f\t%2.1f\n',cycles,max(error));

    if min(1 - max(error)) > 0 %convergence
        break; %exit while loop - solution has converged
    else
        %update initial conditions to previous solution to solve again
        x0 = x(length(t),:);
    end  
    cycles = cycles+1;
end
fprintf('Cycles to converge: %g\n', cycles);

x0_out = x0;

% Separate all of the variables for post-processing
    V_la = x(:,1); 
    V_lv = x(:,2); 
    P_sa = x(:,3);
    Q_sa = x(:,4);
    P_st = x(:,5);
    P_sv = x(:,6);  
    V_ra = x(:,7);  
    V_rv = x(:,8);
    P_pa = x(:,9);  
    Q_pa = x(:,10);
    P_pt = x(:,11);
    P_pv = x(:,12);

% Geta all oif the secondary variables (not solved by ODE) by running the
% ODE file (CVS_dot) at the interpolated rate and pulling off the values
    x2 = zeros(length(t),18);
    for i = 1:length(t)
        [xd x2(i,:) LV_par] = CVS_dot_nv(t(i),x(i,:),parameters,HR);
    end
    e_t_la = x2(:,1);
    e_t_lv = x2(:,2);
    e_t_ra = x2(:,3);
    e_t_rv = x2(:,4);
    P_la = x2(:,5);
    P_lv = x2(:,6);
    Q_vad = x2(:,7);
    Q_m = x2(:,8);
    Q_a = x2(:,9);
    Q_st= x2(:,10);
    Q_sv= x2(:,11);     
    P_ra= x2(:,12);
    P_rv= x2(:,13);
    Q_t = x2(:,14);
    Q_p = x2(:,15);
    Q_pt= x2(:,16);
    Q_pv= x2(:,17);
    Q_cor=x2(:,18);
    Q = Q_a + Q_vad;
    QP = Q .* P_sa;

%generate output_data1 array with time and all variable arrays
    [m,n] = size(x);
    n = 12;
    [m2,n2] = size(x2);
    output_data1 = t';
    output_data1(:,2:1+n) = x(:,1:n);
    output_data1(:,n+2:n+1+n2) = x2;

% i1 and i2 are used to average the arrays
% num is used for i2, which is one less than the length of the array
% because the data is periodic
    i1 = 1;
    i2 = num; 
% End-systolic LV volume (LVESV)
    LVESV = min(V_lv(i1:i2)); 
% Average coronary flow
    Qcor_avg = mean(Q_cor(i1:i2));
% Mean arterial pressure (MAP)
    map = mean(P_sa(i1:i2)); 
% Left ventricular stroke volume
    LVSV = max(V_lv(i1:i2)) - min(V_lv(i1:i2));  
% Average VAD flow
    Q_vad_avg = mean(Q_vad(i1:i2)); 
% Cardiac Output
    CO = mean(Q_a(i1:i2))+Q_vad_avg; 
% Mean right atrial pressure
    P_ra_avg = mean(P_ra(i1:i2)); 
% Mean left atrial pressure
    P_la_avg = mean(P_la(i1:i2)); 
% Mean left ventricular pressure
    P_lv_avg = mean(P_lv(i1:i2)); 
% Energy equivalent pressure  
    EEP = trapz(t(i1:i2),QP(i1:i2))/trapz(t(i1:i2),Q(i1:i2)); 
% Surplus hemodynamic energy
    SHE = EEP - trapz(t(i1:i2),P_sa(i1:i2))/(t(i2)-t(i1));  
% Left ventriculaqr pressure volume area (PVA) and external work (EW)
    [EW PVA] = findEW(V_lv(i1:i2+1),P_lv(i1:i2+1),V_rv(i1:i2+1),P_rv(i1:i2+1),LV_par,RR); 
% Left ventricular end diatolic volume and pressure
    [LVEDV imax] = max(V_lv(i1:i2+1));
    LVEDP = P_lv(imax);
% Cardiac oxygen supply delivered to coronaries (CO2) and used by the
% ventricle (VO2)
    CO2 = 0.85 * 0.14 * Qcor_avg*(t(i2+1)-t(i1));
    E_es_lv = LV_par(1);
    VO2 = 0.000018*PVA + 0.0024*E_es_lv + 0.014; 

% Critical parameters that were averaged are output to a separate array
% called output_data2 that can be put into a spreadsheet
    output_data2(1) = failure;
    output_data2(2) = flow_type;
    if flow_type ~= 1
        output_data2(3) = 0;
    else
        output_data2(3) = asynch_setting;
    end
    if flow_type ~= 2
        output_data2(4) = 0;
        output_data2(5) = 0;
        output_data2(6) = 0;
    else
        output_data2(4) = SV;
        output_data2(5) = Tvad(1);
        output_data2(6) = Tvad(2)-Tvad(1);   
    end
    output_data2(7) = HR;
    output_data2(8) = CO*60/1000; 
    output_data2(9) = max(P_sa(i1:i2));  
    output_data2(10) = min(P_sa(i1:i2));  
    output_data2(11) = map;  
    output_data2(12) = max(P_sa(i1:i2))-min(P_sa(i1:i2));
    output_data2(13) = mean(Q_a(i1:i2))*60/1000;
    output_data2(14) = Q_vad_avg*60/1000;
    output_data2(15) = Qcor_avg;
    output_data2(16) = SHE;
    output_data2(17) = EEP;
    output_data2(18) = VO2;
    output_data2(19) = LVEDV;
    output_data2(20) = LVESV;
    output_data2(21) = (LVEDV - LVESV)/LVEDV*100;
    output_data2(22) = LVEDP;
    output_data2(23) = mean(P_la(i1:i2));
    output_data2(24) = mean(P_ra(i1:i2));
    output_data2(25) = max(P_pa(i1:i2));
    output_data2(26) = min(P_pa(i1:i2));
    output_data2(27) = mean(P_pa(i1:i2));
    output_data2(28) = PVA;
    output_data2(29) = max(P_lv(i1:i2));
    output_data2(30) = min(P_lv(i1:i2));
    output_data2(31) = mean(P_lv(i1:i2));