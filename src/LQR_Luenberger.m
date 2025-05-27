function [L, K_new, V_new] = LQR_Luenberger(sys)
%% Observer
% The observer should have a fast convergence speed and have a significant influence to the measurement noise.
% Assigning observer poles
ew_L= [-50 -51 -52 -53 -59 -55];

L = place(sys.A', sys.C',ew_L)'; 


%% Parameters for LQR
gw = [100; 70; 1; 70; 1; 1];    % Weight coefficient for state vector
Q = diag(gw);

R = [400 0;0 400];              % Weight coefficient for input vector

[K, P, e] = lqr(sys.A, sys.B, Q, R);

% Feedforward controller
V = -pinv(sys.C*inv(sys.A - sys.B*K)*sys.B);


%% Parameters for LQI
A_new = [sys.A, zeros(6,2); -sys.C(1:2,:), zeros(2,2)];
B_new = [sys.B; zeros(2,2)];

gw_new = [100000; 100000; 500000; 100000; 500; 100000;100000;50000];

Q_new = diag(gw_new);

[K_new, P_new, e_new] = lqr(A_new, B_new, Q_new, R);

% New feedforward controller for LQI
V_new = -pinv(sys.C*inv(sys.A - sys.B*K_new(:,1:6))*sys.B);

end


