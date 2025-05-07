function [sys, m, n] = sys_Modelling(var, Moments)

% Linearisierung
x = [var.alpha; var.beta; var.gamma; var.alphadot; var.betadot; var.gammadot];   % State quantity
u = [var.Ff; var.Fb];                       % input vector of the system
y = [var.alpha; var.beta; var.gamma];       % output vector of the system
x_s = [0; 0*pi/180; 0; 0; 0; 0];            % stable point of the system

% Get the derivative of the state quantity(left side of the State space eq.)
x_dot = [var.alphadot; var.betadot; var.gammadot; Moments.alphadot_2; Moments.betadot_2; Moments.gammadot_2];

% Find the derivative of the states at the stable point
x_dot_s = subs(x_dot, x, x_s);

% Solve the state space equation at the stable point to get the input value
u_stable = solve(x_dot_s==zeros(6,1));

u_s = zeros(2,1);

if isempty(u_stable.Ff)    
    u_s(1)=0; 
else
    u_s(1)=u_stable.Ff;
end    
if isempty(u_stable.Fb)    
    u_s(2)=0; 
else
    u_s(2)=u_stable.Fb;
end

% Jacobian derivatives of the system equations with respect to the state quantities 
J1 = jacobian(x_dot, x);
A = subs(J1, [var.alpha, var.beta, var.gamma, var.Ff, var.Fb], [0, 0*pi/180, 0, u_s(1), u_s(2)]);

% Jacobian derivatives of the system equations with respect to the input quantities
J2 = jacobian(x_dot, u);
B = subs(J2, [var.alpha, var.beta, var.gamma], [0,0*pi/180,0]);

C = jacobian(y, x);
D = jacobian(y, u);

A = double(A);
B = double(B);
C = double(C);
D = double(D);
sys = ss(A,B,C,D);

CONT=ctrb(A,B); % controllability
m = rank(CONT);
OBSER=obsv(A,C);% Observability 
n = rank(OBSER);

end

