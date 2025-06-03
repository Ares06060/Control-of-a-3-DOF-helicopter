function ref = ref_trajectory_planning(time)
% Trajectory planning for alpha, beta, gamma and U_magnet over time t
N = length(time);
ref.gamma = zeros(1, N);
x = zeros(1, N);
y = zeros(1, N);

for i = 1:N

    t = time(i);

    if t <= 5
        ts = 0; te = 5;
        a0 = five_times_plan(ts, -27, 0, 0, te, 0, 0, 0);
        beta = polyval(flip(a0'), t);
        alpha = 0;

    elseif t > 5 && t <= 10
        ts = 5; te = 10;
        a0 = five_times_plan(ts, 0, 0, 0, te, 90, 0, 0);
        alpha = polyval(flip(a0'), t);
        beta = 0;

    elseif t > 10 && t <= 15
        ts = 10; te = 15;
        a0 = five_times_plan(ts, 0, 0, 0, te, -22, 0, 0);
        beta = polyval(flip(a0'), t);
        alpha = 90;

    elseif t > 15 && t <= 17
        alpha = 90;
        beta = -22;

    elseif t > 17 && t <= 22
        ts = 17; te = 22;
        a0 = five_times_plan(ts, -22, 0, 0, te, 0, 0, 0);
        beta = polyval(flip(a0'), t);
        alpha = 90;

    elseif t > 22 && t <= 42
        ts = 22; te = 42;
        a0 = five_times_plan(ts, 90, 0, 0, te, 450, 0, 0);
        alpha = polyval(flip(a0'), t);
        beta = 0;

    elseif t > 42 && t <= 47
        ts = 42; te = 47;
        a0 = five_times_plan(ts, 0, 0, 0, te, -22, 0, 0);
        beta = polyval(flip(a0'), t);
        alpha = 450;

    elseif t > 47 && t <= 49
        alpha = 450;
        beta = -22;

    elseif t > 49 && t <= 54
        ts = 49; te = 54;
        a0 = five_times_plan(ts, -22, 0, 0, te, 0, 0, 0);
        beta = polyval(flip(a0'), t);
        alpha = 450;

    elseif t > 54 && t <= 79
        ts = 54; te = 79;
        a0 = five_times_plan(ts, 450, 0, 0, te, 0, 0, 0);
        alpha = polyval(flip(a0'), t);
        beta = 0;

    elseif t > 79 && t < 84
        ts = 79; te = 84;
        a0 = five_times_plan(ts, 0, 0, 0, te, -27, 0, 0);
        beta = polyval(flip(a0'), t);
        alpha = 0;

    else
        alpha = 0;
        beta = -27;
    end

    % Output
    x(i) = alpha;
    y(i) = beta;
    
end

ref.alpha = x;
ref.beta = y;

% Five-times polynomial planning
function a = five_times_plan(ts, xs, vs, as, te, xe, ve, ae)
    para = [xs; vs; as; xe; ve; ae];
    Tran = [1,   ts,  ts^2,  ts^3,    ts^4,     ts^5;
            0,   1,   2*ts,  3*ts^2,  4*ts^3,   5*ts^4;
            0,   0,   2,     6*ts,    12*ts^2,  20*ts^3;
            1,   te,  te^2,  te^3,    te^4,     te^5;
            0,   1,   2*te,  3*te^2,  4*te^3,   5*te^4;
            0,   0,   2,     6*te,    12*te^2,  20*te^3];
    a = Tran \ para;  
end

end
