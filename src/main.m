clear
clc

%% Parameters of the helicopter
% length
len.l1 = 0.655;         % Distance between travel axis to helicopter body  
len.l2 = 0.262;         % Length of the secondary arm
len.l3 = 0.1355;        % Force arm of the gravity of secondary arm
len.l4 = 0.2285;        % Length from Joint 1 to crossing point
len.l5 = 0.21;          % Length from Joint 1 to midpoint of main arm 
len.l6 = 0.355/2;       % Distance between pitch to each motor 
len.l7 = 0.534;
len.h = 0.042+0.019/2;  % Distance between axis Y and main arm

% Weight
wei.mw = 1.918;             % Mass of the counterweight
wei.marm2 = 0.138;          % Mass of secondary arm
wei.marm1 = 0.377;          % Mass of main arm
wei.mf = 0.661;             % Mass of front motor and bodies
wei.mb = 0.661;             % Mass of back motor and bodies
wei.mj2 = 0.106;            % Mass of the joint2
wei.mmag = 0.07;            % Mass of the magnet
wei.g = 9.81;               % Erdbeschleunigung
wei.theta = deg2rad(15);    % Crossing angle of main and secondary arm
wei.phi = 0.1427;           % Deviation of the helicopter relative to joint2
wei.d = 0.102;              % Distance between joint3 to the rail


%% Calculate the Moment of inertia of 3 axis and their Moments and Angular acceleration 
% Definition of the movements and force in 3 axis
syms alpha beta gamma alphadot betadot gammadot Ff Fb;

var.alpha = alpha;
var.beta = beta;
var.gamma = gamma;
var.alphadot = alphadot;
var.betadot = betadot;
var.gammadot = gammadot;
var.Ff = Ff;
var.Fb = Fb;

% Calculate the moment of inertia Jalpha, Jbeta, Jgamma
Inertia = calc_Moment_of_Inertia(len, wei, var);

% Calculate the Moments and Angular acceleration in 3 axis
Moments = calc_Moment(len, wei, var, Inertia);


%% System modelling and LQR controller, Luenberger observer design
% system modelling (continuous system)
[sys_Linear, m, n] = sys_Modelling(var, Moments);

% Luenberger observer and LQR controller design(针对连续性系统设计的)
[L, K_new, V_new] = LQR_Luenberger(sys_Linear);


%% System modelling and MPC controller, EKF design
% check the system characteristic frequency
poles = pole(sys_Linear);
natural_freq = abs(poles);
max_freq = max(natural_freq);       % 系统带宽fc
Dt = 2*pi/(20*max_freq);            % 采样频率为系统带宽的10-20倍, 计算得到的采样间隔为0.1096

dt = 0.1;                           % 为了简便，取0.1s

% discretize the linear system modell for MPC
sys_MPC = c2d(sys_Linear, dt, 'zoh');       % 对于线性系统，使用c2d精度高且计算效率好

% discretize and then linearize the nonlinear system for EKF
sys_EKF = sys_dis_linear(var, Moments, dt);






