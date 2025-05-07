function Inertia = calc_Moment_of_Inertia(len, wei, var)


%% alpha
aRb = [cos(var.beta) 0 sin(var.beta); 
       0 1 0; 
       -sin(var.beta) 0 cos(var.beta)];

bRc = [cos(var.gamma) -sin(var.gamma) 0; 
       sin(var.gamma) cos(var.gamma) 0; 
       0 0 1];

arb = [len.l1*sin(var.beta); 0; len.l1*cos(var.beta)];

cP_f = [-0.065; -len.l6; 0];
cP_b = [-0.065; len.l6; 0];

aP_f = aRb*bRc*cP_f + arb;
aP_b = aRb*bRc*cP_b + arb;

D_Pf = aP_f(2)^2 + aP_f(3)^2;
D_Pb = aP_b(2)^2 + aP_b(3)^2;


J_mw_alpha = (1/12)*wei.mw*((0.070^2 + 0.057^2)*(sin((pi/2)+var.beta-wei.theta))^2+(0.057^2 + 0.057^2)*(cos((pi/2)+var.beta-wei.theta))^2);
J_marm1_alpha = (1/12)*wei.marm1*((0.915^2+0.019^2)*(cos(var.beta))^2+(0.019^2+0.019^2)*(sin(var.beta))^2);
J_marm2_alpha = (1/12)*wei.marm2*((0.305^2+0.025^2)*(cos(var.beta-wei.theta))^2 +(0.025^2+0.025^2)*(sin(var.beta-wei.theta))^2);
J_j1_alpha = 1/12*0.029*((0.019^2 + 0.025^2)*(cos(var.beta))^2 + (0.019^2 + 0.019^2)*(sin(var.beta))^2);
J_bar_alpha =  1/12*0.322*((0.495^2 + 0.006^2)*((cos(var.beta))^2+(cos(var.gamma))^2) + (0.006^2 + 0.006^2)*(sin(var.beta))^2);
J_joint_alpha =   1/12*0.026*((0.019^2 + 0.025^2)*((cos(var.beta))^2+(cos(var.gamma))^2) + (0.019^2 + 0.019^2)*(sin(var.beta))^2);
J_body_alpha = 1/12*0.2*(3*0.114^2 + 0.027^2)*((sin(var.beta))^2+(cos(var.gamma))^2) + 1/2*0.2*0.114^2*(cos(var.beta))^2; 
J_moter_alpha =   1/12*0.287*(3*0.0195^2 + 0.069^2)*((sin(var.beta))^2+(cos(var.gamma))^2) + 1/2*0.287*0.0195^2*(cos(var.beta))^2;


Inertia.Jalpha = J_mw_alpha + wei.mw*((len.l4*cos(var.beta)+len.l2*cos(wei.theta-var.beta))^2) ...,                 % Moment of CounterWeight
    + J_marm1_alpha + wei.marm1*(len.l5*cos(var.beta))^2 ...,                                               % Moment of main arm
    + J_marm2_alpha + wei.marm2*((len.l3*cos(wei.theta-var.beta)+len.l4*cos(var.beta))^2)...,               % Moment of secondary arm
    + wei.mj2*(len.l1*cos(var.beta))^2 ...,                                                                 % Moment of joint2
    + wei.mmag*(len.l7*cos(var.beta))^2 ...,                                                                % Moment of magnet
    + J_j1_alpha + 0.029*((0.042-0.019/2)*sin(var.beta))^2 ...,                                             % Moment of joint1
    + J_bar_alpha + 0.322*(len.l1*cos(var.beta)+(0.042+0.019+0.028+0.041+0.006/2)*sin(var.beta))^2 ...,     % Moment of bar
    + J_joint_alpha + 0.026*(len.l1*cos(var.beta)+(0.042+0.019+0.028+0.041-0.019/2)*sin(var.beta))^2 ...,   % Moment of joint3
    + 2*J_body_alpha + 2*J_moter_alpha + (0.2+0.287)*(D_Pf+D_Pb);                                           % Moment of Motors and Bodies

Inertia.Jalpha = 1.1794;


%% beta
J_mw = (1/12)*wei.mw*(0.070^2 + 0.057^2);
J_marm1 = (1/12)*wei.marm1*(0.915^2 + 0.019^2);
J_marm2 = (1/12)*wei.marm2*(0.305^2 + 0.025^2);

Inertia.Jbeta = J_mw + wei.mw * ((len.l4+len.l2*cos(wei.theta))^2+(len.l2*sin(wei.theta)-len.h)^2) ...,             % Moment of CounterWeight
        + J_marm1 + wei.marm1 * (((len.l5)/2)^2+len.h^2) ...,                               % Moment of main arm
        + J_marm2 + wei.marm2 * (((len.l3*cos(wei.theta))+len.l4)^2+(len.l3*sin(wei.theta)-len.h)^2) ...,   % Moment of secondary arm
        + (wei.mf+wei.mb)*(len.l1^2+(0.044+len.h)^2) ...,                                       % Moment of heli(simplify)
        + wei.mj2*(len.l1^2+(len.h+0.019+0.028)^2) ...,                                     % Moment of joint 2
        + wei.mmag*(len.l7^2+(len.h+0.062*sin(pi/3))^2);                                    % Moment of Magnet
% 1.1435        


%% gamma
J_body  =   1/12*0.2*(3*0.114^2 + 0.027^2);           % Body = Cylinder
J_moter =   1/12*0.287*(3*0.0195^2 + 0.069^2);        % Motor = Cylinder
J_bar   =   1/12*0.322*(0.495^2 + 0.006^2);           % Bar = Cubic
J_joint =   1/12*0.026*(0.019^2 + 0.025^2);           % joint = Cubic

Inertia.Jgamma = 2*(J_body+0.2*((0.027/2+0.069-0.041)^2+(0.355/2)^2))+ ...,           % Moment of bodies
    2*(J_moter+0.287*((0.041-0.069/2)^2+(0.355/2)^2))+ ...,                     % Moment of Motors
    J_bar+0.322*(0.041+0.006/2)^2+ ...,                                         % Moment of bar
    J_joint+0.026*(0.041-0.019/2)^2;                                            % Moment of joints3
% 0.0402    


end






