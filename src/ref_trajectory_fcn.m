function ref = ref_trajectory_fcn(t)
beta = -27;
if t <= 10
    a0 = five_times_plan(0, -27, 0, 0, 10, 0, 0, 0);
    alpha = 0;
    beta = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    U_magnet = 0;
elseif t > 10 && t <= 18
    a0 = five_times_plan(10, 0, 0, 0, 18, 90, 0, 0);
    alpha = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    beta = 0;
    U_magnet = 0;
elseif t > 18 &&t <= 26
    a0 = five_times_plan(18, 0, 0, 0, 26, -22, 0, 0);
    alpha = 90;
    beta = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    U_magnet = 1;
elseif t > 26 && t <= 30
    alpha = 90;
    beta = -22;
    U_magnet = 0;
elseif t > 30 && t <= 38
    a0 = five_times_plan(30, -22, 0, 0, 38, 0, 0, 0);
    alpha = 90;
    beta = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    U_magnet = 0;
elseif t > 38 && t <= 62
    a0 = five_times_plan(38, 90, 0, 0,62, 450, 0, 0);
    alpha = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    beta = 0;
    U_magnet = 0;
elseif t > 62 && t <= 70
    a0 = five_times_plan(62, 0, 0, 0, 70, -22, 0, 0); 
    alpha = 450;
    beta = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    U_magnet = 1;
elseif t > 70 && t <= 74
    alpha = 450;
    beta = -22;
    U_magnet = 0;
elseif t > 74 && t <= 82
    a0 = five_times_plan(74, -22, 0, 0, 82, 0, 0, 0);
    alpha = 450;
    beta = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    U_magnet = 0;
elseif t > 82 && t <= 110
    a0 = five_times_plan(82, 450, 0, 0, 110, 0, 0, 0);
    alpha = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    beta = 0;
    U_magnet = 0;
elseif t > 110 && t < 118
    a0 = five_times_plan(110, 0, 0, 0, 118, -27, 0, 0);
    alpha = 0;
    beta = a0(1)+a0(2)*t+a0(3)*t^2+a0(4)*t^3+a0(5)*t^4+a0(6)*t^5;
    U_magnet = 1;
else 
    alpha = 0; beta = -27;
    U_magnet = 0;
end
r = [alpha; beta;0];
function a = five_times_plan(ts, xs, vs, as, te, xe, ve, ae)
    para = [xs;vs;as;xe;ve;ae];
    Tran = [1,   ts,  ts^2,  ts^3,    ts^4,     ts^5;
            0,   1,   2*ts,  3*ts^2,  4*ts^3,   5*ts^4;
            0,   0,   2,     6*ts,    12*ts^2,  20*ts^3;
            1,   te,  te^2,  te^3,    te^4,     te^5;
            0,   1,   2*te,  3*te^2,  4*te^3,   5*te^4;
            0,   0,   2,     6*te,    12*te^2,  20*te^3];
    a = inv(Tran)*para;
end

end