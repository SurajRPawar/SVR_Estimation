function xout = est_state(x0,version, mode)
% x0 : The 4 states of the scaled system
% version of estimator being used
% mode = 1 : Populate the xhat0 and send back
% mode = 2 : Populate the estimate hat (Rsvrhat, Caohat, Crhat) and send
% back
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;
if version == 1
    % v1
    % x5 = 1/Rsvr*Cr;
    % x6 = 1/Rsvr*Cao;
    % x7 = 1/Cao;
    if mode == 1; 
        xout = [x0; 1/(Rsvr*Cr); 1/(Rsvr*Cao); 1/Cao];
    elseif mode == 2;
        Rsvrhat = x0(7)/x0(6);
        Caohat = 1/x0(7);
        Crhat = x0(6)/(x0(7)*x0(5));
        xout = [Rsvrhat; Caohat; Crhat];
    end
    
elseif version == 2
    % v2
    % x5 = 1/Rsvr*Cr
    % x6 = 1/Cao;
    if mode == 1
        xout = [x0; 1/(Cr*Rsvr); 1/Cao];
    elseif mode == 2
        Rsvrhat = 1/(Cr*x0(5));
        Caohat = 1/x0(6);
        xout = [Rsvrhat; Caohat];
    end
    
elseif version == 3
    % v3
    % x5 = 1/Rsvr;
    % x6 = 1/Cao;
    % x7 = 1/Cr;
    if mode == 1
        xout = [x0; 1/Rsvr; 1/Cao; 1/Cr];
    elseif mode == 2
        Rsvrhat = 1/x0(5);
        Caohat  = 1/x0(6);
        Crhat   = 1/x0(7);
        xout = [Rsvrhat; Caohat; Crhat];
    end
elseif version == 5;
    % v5
    % x5 = Rsvr;
    if mode == 1;
        xout = [x0; Rsvr];
    elseif mode == 2;
        Rsvrhat = x0(5);
        xout = [Rsvrhat];
    end
end
end