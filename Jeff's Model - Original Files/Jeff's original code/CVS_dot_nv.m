%Function file with the differential equations of the computational
%cardiovascular system model
%Jeff Gohean 7.26.13

function [xdot x2 LV_par] = CVS_dot_nv(t,x,parameters,HR)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Parameters
%   ventricle, atria, valves, systemic circulation, pulmonary circulation
%   values found in journal paper "Verification of a Computational
%   Cardiovascular System Model..." ASAIO, Gohean et al 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ventricules and Atria
    % Ventricular and atrial elastance (mmHg / ml)
        E_es_la = 0.25;
        E_es_lv = 3.25;     
        E_es_ra = 0.25;
        E_es_rv = 0.75;
    % Unstressed ventricule and atria volume (ml)
        V0_la = 10;
        V0_lv = 5;          
        V0_ra = 10;
        V0_rv = 5;
    % Exponential constant for passive elastance (mmHg)
        A_la = 0.45;
        A_lv = 0.03;        
        A_ra = 0.45;
        A_rv = 0.04;
    % Exponential constant for passive elastance (1 / ml)
        B_la = 0.05;
        B_lv = 0.05;        
        B_ra = 0.05;
        B_rv = 0.05;  
    % Equation for systolic time as function of HR
        T_v = (550 - 1.75*HR)/1000;     
        shift = 0.666;                  
        T_a = 0.09;

% Valve resistance (mmHg s/ml)
    R_p = 0.0025;
    R_t = 0.0025;
    R_m = 0.0025;
    R_a = 0.0025;
    L_v = 0.0000003; %hack used to make solver unstiff (valve diodes create a very stiff system)

% Systemic parameters
    R_sa =  0.15;       %arterial resistance (mmHg s / ml)
    L_sa =  0.0022;     %arterial inertia (mmHg s^2 / ml)
    C_sa =  1.25;       %arterial compliance (ml / mmHg)
    R_st =  0.8;        %tree resistance (mmHg s / ml)
    C_st =  2.0;        %tree compliance (ml / mmHg)
    R_sv =  0.025;      %vein resistance (mmmHg s / ml)
    C_sv = 20.0;        %vein compliance (ml / mmHg)

% Pulmonary parameters
    R_pa =  0.07;       %arterial resistance (mmHg s / ml)
    L_pa =  0.0018;     %arterial inertia (mmHg s^2 / ml)
    C_pa =  7.5;        %arterial compliance (ml / mmHg)
    R_pt =  0.04;       %tree resistance (mmHg s / ml)
    C_pt =  0.5;        %tree compliance (ml / mmHg)
    R_pv =  0.003;      %vein resistance (mmmHg s / ml)
    C_pv = 20.0;        %vein compliance (ml / mmHg)

% Heart Failure, change some model parameters to simulate heart failure 
    if failure == 1 %these are the values for "failure" in the journal paper
        E_es_lv = .3;   % reduced LV elastance (weak heart)
        A_lv = 0.18;    % modify A/B exp coef.s to enlarge LV
        B_lv = 0.016;
        E_es_rv = .45;  % reduced RV elastance
        C_sa = 0.65;    % reduced arterial compliance
        R_st = 0.9;     % increased SVR
        C_st = 1.5;     % reduced arterial compliance
        R_pa = 0.05;    % reduced PVR   !not in the journal paper but should have been
        R_pt = 0.195;   % increased PVR !not in the journal paper but should have been
        if flow_type ~= 0  
            R_st = 0.45;   % reduce SVR when VAD in place !not in the journal paper
            R_pt = 0.01;   % reduce PVR when VAD in place !not in the journal paper
        end
    elseif failure == 2
        E_es_lv = 0.45; % reduce LV elastance (active)
        A_lv = 0.25;    % reduce LV elastance (passive)
        B_lv = 0.018;  
        E_es_rv = 0.45; % reduce RV elastance
        C_sa = 0.65;    % reduce arterial compliance
        C_st = 1.5;     % reduce arterial compliance
        R_st = 1.1;     % increased SVR
        R_pt = 0.15;    % increase PVR
        if flow_type ~= 0
            R_st = 0.5;  % use when VAD in place
            R_pt = 0.04; % use when VAD in place
        end    
    elseif failure == 3
        E_es_lv = 0.48; % reduce LV elastance (active)
        A_lv = 0.25;    % reduce LV elastance (passive)
        B_lv = 0.018;  
        E_es_rv = 0.45; % reduce RV elastance
        C_sa = 0.65;    % reduce arterial compliance
        C_st = 1.5;     % reduce arterial compliance
        R_st = .85;     % increased SVR
        R_pt = 0.15;    % increase PVR
        if flow_type ~= 0
            R_st = 0.5;  % use when VAD in place
            R_pt = 0.04; % use when VAD in place
        end   
    end
% Output array for calculating PVA and EW (findEW.m)
    LV_par(1) = E_es_lv; 
    LV_par(2) = A_lv;
    LV_par(3) = B_lv;
    LV_par(4) = V0_lv;
    LV_par(5) = T_v;
    LV_par(6) = shift;
    
% Inflow / Outflow cannulas
    % Blood properties (density and viscosity)
        rho = 1050;     %kg/m^3
        mu  = 0.0035;   %Pa s
    % Inflow/Outflow cannula lengths (inches)
        Lin  = 4;       % Inflow
        Lout = 8;       % Outflow
    % Inflow/Outflow cannula radius (inches)
        Rin  = 5/8/2;   % Inflow
        Rout = 5/8/2;   % Outflow
    % Convert lengths and radii to meters
        Lin=Lin*.0254; Lout=Lout*.0254; Rin=Rin*.0254; Rout=Rout*.0254;
    % Cannula Resistance and Inertial 
        R_in  = 8*pi*Lin*mu/(pi*Rin^2)^2;
        L_in  = rho*Lin/(pi*Rin^2);
        R_out = 8*pi*Lout*mu/(pi*Rout^2)^2;
        L_out = rho*Lout/(pi*Rout^2);
    % Lump Inflow/Outflow Resistance and Inertia
        R_tot = (R_in + R_out)/(100^3*133.28);
        L_tot = (L_in + L_out)/(100^3*133.28);
    
% VAD Parameters
    % Heartmate II parameters for a limited number of rpms 
        % 8.4k, 9k, 9.6k rpm, 8, 7.2, 7.1, 7.3
        hm1_v = -1*[.41 .37 .33, .47, .570, .583, .557];
        hm2_v = [4.6 4.5 4.3, 4.8, 4.86, 4.86, 4.86];
        hm3_v = -1*[23 24 25, 22, 20, 20, 20];
        hm4_v = [96 108 123, 85, 62.38, 60, 65.19];
        hm_c1 = 0; hm_c2 = 0; hm_c3 = 0; hm_c4 = 0; %init to zero
    % TORVAD parameters      
        T_vad_1 = Tvad(1);
        T_vad_2 = Tvad(2);
        T_vad_3 = Tvad(3);
        T_vad_4 = Tvad(4);
        acc_vad_A = 0.5;
        dec_vad_A = 0.5;
        acc_vad_B = 0.5;
        dec_vad_B = 0.5;
        plat_vad_A = 1.0 - acc_vad_A - dec_vad_A;
        plat_vad_B = 1.0 - acc_vad_B - dec_vad_B;
        diff_vad_A = T_vad_2 - T_vad_1;
        diff_vad_B = T_vad_4 - T_vad_3;
        Q_max_vad_A = SV / (plat_vad_A + 0.5*acc_vad_A + 0.5*dec_vad_A) / diff_vad_A;
        Q_max_vad_B = SV / (plat_vad_B + 0.5*acc_vad_B + 0.5*dec_vad_B) / diff_vad_B;
        
    % Set VAD parameters based on type of assistance
        if flow_type == 0
            % No VAD
            SV = 0;
        elseif flow_type == 1   
            % HM II support
            % set HM II coefficients based on asynch speed setting
            % set TORVAD SV to zero
            if asynch_setting <= length(hm1_v)
                hm_c1 = hm1_v(asynch_setting);
                hm_c2 = hm2_v(asynch_setting);
                hm_c3 = hm3_v(asynch_setting);
                hm_c4 = hm4_v(asynch_setting);
            end
            SV = 0;
        else
            % TORVAD suppport - SV already defined
        end
    
% Heart rate and RR interval
    RR = 60/HR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the arrays used for ODE solver
    xdot = zeros(size(x));
    x2 = zeros(16,1);

% Normalize time for the cardiac cycle (t=0 at R-wave)
    tn = mod(t,RR);

    V_la = x(1);    %Left atrial volume
    V_lv = x(2);    %Left ventricular volume
    P_sa = x(3);    %Systemic artery pressure
    Q_sa = x(4);    %Systemic artery flow
    P_st = x(5);    %Systemic artery tree pressure
    P_sv = x(6);    %Systemic vein pressure
    V_ra = x(7);    %Right arterial volume
    V_rv = x(8);    %Right ventricular volume
    P_pa = x(9);    %Pulmonary artery pressure
    Q_pa = x(10);   %Pulmonary artery flow
    P_pt = x(11);   %Pulmonary artery tree pressure
    P_pv = x(12);   %Pulmonary vein pressure
    Q_m  = x(13);   %Mitral valve flow
    Q_a  = x(14);   %Aortic valve flow
    Q_t  = x(15);   %Tricuspid valve flow
    Q_p  = x(16);   %Pulmonary valve flow
    Q_vad= x(17);   %VAD flow
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Systemic Circulation equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Left atrial pressure

% Left atria pressure
    % Normalized elastance (min 0, max 1) 
        if tn >= RR-T_a && tn <= RR
            e_t_la = 1/2 * (1 - cos(2*pi*(tn-(RR-T_a))/T_a));
        else
            e_t_la = 0;
        end
    % Passive elastance pressure
        P_la_passive = (1-e_t_la) * A_la * (exp(B_la*(V_la-V0_la)) - 1);
    % Active elastance pressure
        P_la_active =  + (e_t_la)*E_es_la*(V_la-V0_la);
    % Left atria pressure
        P_la = P_la_passive + P_la_active;
    
% Left ventricular pressure
    % Normalized elastance (min 0, max 1)
        if tn <= T_v*shift
            e_t_lv = 1/2 * (1 - cos(pi*tn/(T_v*shift)));
        elseif tn > T_v*shift && tn <= T_v
            e_t_lv = 1/2 * (1 + cos(pi*(tn-T_v * shift)/(T_v*(1-shift))));
        else
            e_t_lv = 0;
        end
    % Passive elastance pressure
        P_lv_passive = (1-e_t_lv) * A_lv * (exp(B_lv*(V_lv - V0_lv)) - 1);
    % Active elastance pressure
        P_lv_active = (e_t_lv)*E_es_lv*(V_lv-V0_lv);
    % Left ventricular pressure
        P_lv = P_lv_active + P_lv_passive;
    
% VAD flow
    
    % No VAD
    if flow_type == 0
        Q_vad = 0;
        Q_vad_dot = 0;
        
    % HMII assistance
    elseif flow_type == 1
        
        % Pressure as a funciton of flow (Q) - from HMII PQ curves
            dpfq =   hm_c1*(Q_vad*60/1000)^3 ...
                   + hm_c2*(Q_vad*60/1000)^2 ...
                   + hm_c3*(Q_vad*60/1000) ...
                   + hm_c4;

        % Increase the resistance at volumles less than 30ml, when the
        % heart diameter is the diameter of the inflow cannula and suction
        % starts to occur
            if V_lv < 30 && V_lv > 20
                R_tot = R_tot + (30 - V_lv)/10 * 20;
            elseif V_lv < 20
                R_tot = R_tot + 20;
            end
   
        % ODE for CF VAD - uses PQ curve and cannula resistance and inertia
            Q_vad_dot = (dpfq - (P_sa-P_lv) - R_tot * Q_vad)/L_tot;
            
    % TORVAD assistance
    elseif flow_type == 2
        
        % Correct normalized time if the stroke pushes through next R-wave
            if T_vad_2 > RR && tn < (T_vad_2 - RR)
                tn = tn + RR;
            end
            if T_vad_4 > RR && tn < (T_vad_4 - RR)
                tn = tn + RR;
            end
            
        % First VAD ejection   
        if tn <= T_vad_2 && tn > T_vad_1
            if tn <= T_vad_1 + acc_vad_A * diff_vad_A
                Q_vad = Q_max_vad_A* 1/2 * (1 - cos(pi*(tn-T_vad_1)/(acc_vad_A * diff_vad_A)));
            elseif tn >= T_vad_1 + acc_vad_A * diff_vad_A && tn <= T_vad_1 + (plat_vad_A + acc_vad_A) * diff_vad_A
                Q_vad = Q_max_vad_A;
            else
                Q_vad = Q_max_vad_A * 1/2 * (1 + cos(pi*(tn-T_vad_1-diff_vad_A*(acc_vad_A+plat_vad_A))/(dec_vad_A * diff_vad_A)));
            end

        % Between strokes
        elseif tn <= T_vad_4 && tn > T_vad_3
            if tn <= T_vad_3 + acc_vad_A * diff_vad_B
                Q_vad = Q_max_vad_B* 1/2 * (1 - cos(pi*(tn-T_vad_3)/(acc_vad_A * diff_vad_B)));
            elseif tn >= T_vad_3 + acc_vad_A * diff_vad_B && tn <= T_vad_3 + (plat_vad_A + acc_vad_A) * diff_vad_B
                Q_vad = Q_max_vad_B;
            else
                Q_vad = Q_max_vad_B * 1/2 * (1 + cos(pi*(tn-T_vad_3-diff_vad_B*(acc_vad_A+plat_vad_A))/(dec_vad_A * diff_vad_B)));
            end       
%             Q_vad = Q_vad*-1;
        
        else
            Q_vad = 0;
        end

        
        % positive displacement, flow is known directly, Q_dot = 0
        Q_vad_dot = 0;
    end

% Mitral valve flow   
    Q_m_dot = (P_la - P_lv - R_m^2*Q_m^2)/L_v;
    if Q_m <= 0
        Q_m = 0;
        if P_la < P_lv
            Q_m_dot = 0;
        end
    end

% Aortic valve flow
    Q_a_dot = (P_lv - P_sa - R_a^2*Q_a^2)/L_v;
    if Q_a < 0
        Q_a = 0;
        if P_lv < P_sa
            Q_a_dot = 0;
        end
    end  
         
% Left atrial and ventricular volume
    % if simplified model, set the pulmonary vein pressure to preload
        if open == 1 
            P_pv = pre;
        end
    % Pulmonary vein flow    
        Q_pv = (P_pv - P_la) / R_pv; %right side return, flow into LA
    % Left atrial volume ODE
        V_la_dot = Q_pv - Q_m;
    % Left ventricular volume ODE
        V_lv_dot = Q_m - Q_a - Q_vad; 
        
% Systemic artery pressure and flow ODEs
    P_sa_dot = (Q_a + Q_vad - Q_sa)/C_sa;
    Q_sa_dot = (P_sa - P_st - Q_sa*R_sa)/L_sa;
    
% Systemic Arteries (sa) Pressure and Flow
    % if simplified model, set the afterload resistance
        if open == 1
            R_st = after - R_sa;
        end
    % Systemic artery tree flow and pressure ODE
    Q_st = (P_st - P_sv) / R_st;
    P_st_dot = (Q_sa  - Q_st) / C_st;   
    
% Coronary flow
    if tn > T_v %diastole
        R_cor = 37;
    else %systole
        R_cor = 80;
    end
    Q_cor = (P_sa - P_sv) / R_cor;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulmonary Circulation equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Right atria pressure
    % Normalized elastance (min 0, max 1) 
        e_t_ra = e_t_la;
    % Passive elastance pressure
        P_ra_passive = (1-e_t_ra) * A_ra * (exp(B_ra*(V_ra-V0_ra)) - 1);
    % Active elastance pressure
        P_ra_active =  + (e_t_ra)*E_es_ra*(V_ra-V0_ra);
    % Left atria pressure
        P_ra = P_ra_passive + P_ra_active;

% Right ventricle pressure
    % Normalized elastance (min 0, max 1) 
        e_t_rv = e_t_lv;
    % Passive elastance pressure
        P_rv_passive = (1-e_t_rv) * A_rv * (exp(B_rv*(V_rv-V0_rv)) - 1);
    % Active elastance pressure
        P_rv_active =  + (e_t_rv)*E_es_rv*(V_rv-V0_rv);
    % Left atria pressure
        P_rv = P_rv_passive + P_rv_active;
        
% Tricuspid valve flow
    Q_t_dot = (P_ra - P_rv - R_t^2*Q_t^2)/L_v;
    if Q_t <= 0
        Q_t = 0;
        if P_ra < P_rv
            Q_t_dot = 0;
        end
    end

% Pulmonary valve flow   
    Q_p_dot = (P_rv - P_pa - R_p^2*Q_p^2)/L_v;
    if Q_p <= 0
        Q_p = 0;
        if P_rv < P_pa
            Q_p_dot = 0;
        end
    end  
    
% Right atrial and ventricular volume
    % Systemic vein flow    
        Q_sv = (P_sv - P_ra) / R_sv; %right side return, flow into RA
    % Right atrial volume ODE
        V_ra_dot = Q_sv - Q_t;     
    % Right ventricular volume ODE
        V_rv_dot = Q_t - Q_p;
     
% Pulmonary artery pressure and flow ODEs
    P_pa_dot = (Q_p - Q_pa)/C_pa;
    Q_pa_dot = (P_pa - P_pt - Q_pa*R_pa)/L_pa;
    
% Pulmonary artery tree flow and pressure ODE
    Q_pt = (P_pt - P_pv) / R_pt;
    P_pt_dot = (Q_pa - Q_pt) / C_pt;    
   
% Pulmonary vein pressure ODE
    P_pv_dot = (Q_pt - Q_pv)/C_pv;
    
% Systemic vein pressure ODE - had to place after Q_sv is defined
    P_sv_dot = (Q_st - Q_sv)/C_sv;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
xdot(1) = V_la_dot;    
xdot(2) = V_lv_dot;    
xdot(3) = P_sa_dot;  
xdot(4) = Q_sa_dot;
xdot(5) = P_st_dot; 
xdot(6) = P_sv_dot;    
xdot(7) = V_ra_dot;        
xdot(8) = V_rv_dot;        
xdot(9)= P_pa_dot;  
xdot(10)= Q_pa_dot;
xdot(11)= P_pt_dot; 
xdot(12)= P_pv_dot;
xdot(13)= Q_m_dot;
xdot(14)= Q_a_dot;
xdot(15)= Q_t_dot;
xdot(16)= Q_p_dot;
xdot(17)= Q_vad_dot;


%data array x2 (pass out variables not computed using ODE)
x2(1)=e_t_la;
x2(2)=e_t_lv;
x2(3)=e_t_ra;
x2(4)=e_t_rv;
x2(5)=P_la;
x2(6)=P_lv;
x2(7)=Q_vad;
x2(8)=Q_m;
x2(9)=Q_a;
x2(10)=Q_st;
x2(11)=Q_sv;
x2(12)=P_ra;
x2(13)=P_rv;
x2(14)=Q_t;
x2(15)=Q_p;
x2(16)=Q_pt;
x2(17)=Q_pv;
x2(18)=Q_cor;





