T = T_f;
P = P_f;

a(1) = popul(1)                   * rates(1);          % O      * k0
a(2) = popul(4)                   * rates(2);          % OX2    * k1
a(3) = popul(2) * abs(popul(2)-1) * rates(4);          % X(X-1) * chi
a(4) = popul(3)                   * rates(7)*rates(4); % X2     * beta*chi
a(5) = popul(1) * popul(3)        * rates(5);          % O * X2 * phi
a(6) = popul(4)                   * rates(6)*rates(5); % OX2    * alfa*phi
a(7) = popul(5) * popul(3)        * rates(3);          % Yi*X2  * k2
a(8) = popul(2) * popul(6)        * rates(8);          % X * Y  * lamda1
a(9) = popul(6)                   * rates(9);          % Y      * lamda2

Dt = ( P - T ) ./ a;
[D, m] = min(Dt);
