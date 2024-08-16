%% Inter-arrival time error
function [error] = calc_IAT_err(x, y, N_r)
    % x: occurrence times set 1: Un
    % y: occurrence times set 2: DS

    error = NaN(5,1);

    % inter-arrival times of reaction 1 (r1) of set 1 (x)
    IAT_r1_x = x(2:N_r(1), 1) - x(1:N_r(1)-1, 1);
    IAT_r1_y = y(2:N_r(1), 1) - y(1:N_r(1)-1, 1);
    error(1)  = norm(IAT_r1_x - IAT_r1_y);

    IAT_r2_x = x(2:N_r(2), 2) - x(1:N_r(2)-1, 2);
    IAT_r2_y = y(2:N_r(2), 2) - y(1:N_r(2)-1, 2);
    error(2)  = norm(IAT_r2_x - IAT_r2_y);

    IAT_r7_x = x(2:N_r(3), 3) - x(1:N_r(3)-1, 3);
    IAT_r7_y = y(2:N_r(3), 3) - y(1:N_r(3)-1, 3);
    error(3)  = norm(IAT_r7_x - IAT_r7_y);

    IAT_r8_x = x(2:N_r(4), 4) - x(1:N_r(4)-1, 4);
    IAT_r8_y = y(2:N_r(4), 4) - y(1:N_r(4)-1, 4);
    error(4)  = norm(IAT_r8_x - IAT_r8_y);

    IAT_r9_x = x(2:N_r(5), 5) - x(1:N_r(5)-1, 5);
    IAT_r9_y = y(2:N_r(5), 5) - y(1:N_r(5)-1, 5);
    error(5)  = norm(IAT_r9_x - IAT_r9_y);

end