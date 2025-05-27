function sys_EKF = sys_dis_linear(var, Moments, dt)
% define the system state and input/ouput parameters
u = [var.Ff; var.Fb];
x = [var.alpha; var.beta; var.gamma; var.alphadot; var.betadot; var.gammadot];
y = [var.alpha; var.beta; var.gamma];

% define the system equation
f = @(x, u) [
    x(4);
    x(5);
    x(6);
    Moments.alphadot_2;
    Moments.betadot_2;
    Moments.gammadot_2
];

h = @(x) [
    x(1);
    x(2);
    x(3);
];

% RK4方法实现状态转移
x_k+ = @(x, u, dt) rk4_integration(x, u, dt, f);

% 在当前工作点进行线性化
x_k = x;  % 当前状态
u_k = u;  % 当前输入

%% 计算雅可比矩阵 - 状态转移矩阵 F
% 数值求导计算 F = ∂f/∂x
epsilon = 1e-8;
n_states = length(x_k);
F = zeros(n_states, n_states);

% 计算离散化后的状态转移矩阵
x_k_plus_1 = rk4_step(x_k, u_k, dt);  % 参考轨迹

for i = 1:n_states
    x_pert = x_k;
    x_pert(i) = x_pert(i) + epsilon;
    x_pert_plus_1 = rk4_step(x_pert, u_k, dt);
    F(:, i) = (x_pert_plus_1 - x_k_plus_1) / epsilon;
end

%% 计算控制输入矩阵 G
% G = ∂f/∂u (离散化后)
n_inputs = length(u_k);
G = zeros(n_states, n_inputs);

for i = 1:n_inputs
    u_pert = u_k;
    u_pert(i) = u_pert(i) + epsilon;
    x_pert_plus_1 = rk4_step(x_k, u_pert, dt);
    G(:, i) = (x_pert_plus_1 - x_k_plus_1) / epsilon;
end

%% 计算观测矩阵 H
% H = ∂h/∂x
n_outputs = length(y);
H = zeros(n_outputs, n_states);

y_k = h(x_k);  % 参考观测
for i = 1:n_states
    x_pert = x_k;
    x_pert(i) = x_pert(i) + epsilon;
    y_pert = h(x_pert);
    H(:, i) = (y_pert - y_k) / epsilon;
end

%% 构建EKF系统结构体
sys_EKF.F = F;              % 离散化状态转移矩阵
sys_EKF.G = G;              % 离散化控制输入矩阵
sys_EKF.H = H;              % 观测矩阵
sys_EKF.f_nonlinear = f;    % 非线性状态方程
sys_EKF.h_nonlinear = h;    % 非线性观测方程
sys_EKF.rk4_step = rk4_step; % RK4离散化函数
sys_EKF.dt = dt;            % 采样时间
sys_EKF.n_states = n_states;
sys_EKF.n_inputs = n_inputs;
sys_EKF.n_outputs = n_outputs;
sys_EKF.x_ref = x_k;        % 参考状态点
sys_EKF.u_ref = u_k;        % 参考输入点

end

%% RK4积分函数
function x_next = rk4_integration(x, u, dt, f_func)
    % 四阶龙格-库塔方法
    k1 = f_func(x, u);
    k2 = f_func(x + dt/2 * k1, u);
    k3 = f_func(x + dt/2 * k2, u);
    k4 = f_func(x + dt * k3, u);
    
    x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% 使用示例和说明
% 
% 使用方法:
% 1. 调用函数获取线性化系统
%    sys_EKF = sys_dis_linear(var, Moments, dt);
%
% 2. 在EKF中使用:
%    - 预测步骤: x_pred = sys_EKF.rk4_step(x_est, u, dt);
%    - 协方差预测: P_pred = sys_EKF.F * P * sys_EKF.F' + Q;
%    - 卡尔曼增益: K = P_pred * sys_EKF.H' / (sys_EKF.H * P_pred * sys_EKF.H' + R);
%    - 状态更新: x_est = x_pred + K * (y_meas - sys_EKF.h_nonlinear(x_pred));
%    - 协方差更新: P = (eye(6) - K * sys_EKF.H) * P_pred;
%
% 注意事项:
% 1. 线性化是在当前工作点进行的，需要在每个时间步重新计算
% 2. Moments结构体中的力矩项可能依赖于状态和输入，需要相应更新


