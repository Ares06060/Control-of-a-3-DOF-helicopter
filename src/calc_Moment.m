function Moments = calc_Moment(len, wei, var, Inertia)

%% MOMENT
%% alpha
Moments.Malpha = -(var.Ff + var.Fb)*sin(var.gamma)*len.l1*cos(var.beta);

%% beta
M_all = (wei.mw*wei.g*(len.l2*cos(wei.theta)+len.l4) + wei.marm2*wei.g*(len.l3*cos(wei.theta)+len.l4)..., 
        - (wei.mf + wei.mb + wei.mj2)*wei.g*len.l1 - wei.marm1*wei.g*len.l5 -wei.mmag*wei.g*len.l7)*cos(0);

mass_position = M_all / ((wei.mw + wei.marm2 + wei.marm1 + wei.mmag + wei.mj2 + wei.mf + wei.mb)*wei.g);

Moments.Mbeta = (var.Ff + var.Fb)*len.l1*cos(var.gamma) + (wei.mw + wei.marm2 + wei.marm1 + wei.mmag + wei.mj2...,
        + wei.mf + wei.mb)*wei.g*(mass_position*cos(var.beta) - len.h*sin(var.beta));

%% gamma
Moments.Mgamma = (var.Ff - var.Fb)*len.l6 + wei.mb*wei.g*cos(wei.phi + var.gamma)*len.l6/cos(wei.phi)...,
                 - wei.mf*wei.g*cos(wei.phi - var.gamma)*len.l6/cos(wei.phi);

%% SOLVE
Moments.alphadot_2 = Moments.Malpha/Inertia.Jalpha;
Moments.betadot_2 = Moments.Mbeta/Inertia.Jbeta;       
Moments.gammadot_2 = Moments.Mgamma/Inertia.Jgamma;

end
