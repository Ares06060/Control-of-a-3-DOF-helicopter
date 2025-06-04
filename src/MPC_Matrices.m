function [E, H] = MPC_Matrices(A, B, Q, R, F, N)

n = size(A, 1);  % A is an n x n matrix, get n
p = size(B, 2);  % B is an n x p matrix, get p

%%
M = [eye(n); zeros(N*n, n)]; % Initialize M matrix. M is a (N+1)n x n matrix,
                             % with an n x n "I" on the top. This step writes
                             % the lower half as 0

C = zeros((N+1)*n, N*p); % Initialize C matrix, this step makes it (N+1)n x NP zeros

% Define M and C
tmp = eye(n);  % Define an n x n I matrix

% Update M and C
for i = 1:N % Loop, i from 1 to N
    rows = i*n + (1:n); % Define current row number, starting from i x n, for n rows
    C(rows, :) = [tmp*B, C(rows-n, 1:end-p)]; % Fill the C matrix
    tmp = A*tmp; % Multiply tmp by A each time
    M(rows, :) = tmp; % Fill the M matrix
end

% Define Q_bar and R_bar
Q_bar = kron(eye(N), Q);
Q_bar = blkdiag(Q_bar, F);
R_bar = kron(eye(N), R);

% Compute G, E, H
G = M' * Q_bar * M; % G: n x n
E = C' * Q_bar * M; % E: NP x n
H = C' * Q_bar * C + R_bar; % H: NP x NP

end