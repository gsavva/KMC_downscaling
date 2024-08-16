%% Modified Next Reaction Method
function [D,m,T,P,a,Dt] = Mod_NRM(D,Dt,T,P,a,m,popul,rates,map,GEN)
    % Random Number Generators, GEN, are automatically updated
    % Therefore, I don't have to pass them out (which slows down execution)

    % T, a are vectors. D is scalar
    % same performance as explict loop below
    T = T + a * D; 
%     for i=1:9
%         T(i) = T(i) + a(i)*D;
%     end
    if m>2 && m<7 % m= 3, 4, 5, 6
        P(m) = P(m) - log(rand());
    else
        P(m) = P(m) - log(rand(GEN{map(m)}));
    end

    % recalculate propensities
    a(1) = popul(1)                   * rates(1);          % O      * k0
    a(2) = popul(4)                   * rates(2);          % OX2    * k1
    a(3) = popul(2) * abs(popul(2)-1) * rates(4);          % X(X-1) * chi
    a(4) = popul(3)                   * rates(7)*rates(4); % X2     * beta*chi
    a(5) = popul(1) * popul(3)        * rates(5);          % O * X2 * phi
    a(6) = popul(4)                   * rates(6)*rates(5); % OX2    * alfa*phi
    a(7) = popul(5) * popul(3)        * rates(3);          % Yi*X2  * k2
    a(8) = popul(2) * popul(6)        * rates(8);          % X * Y  * lamda1
    a(9) = popul(6)                   * rates(9);          % Y      * lamda2

    Dt(1) = ( P(1) - T(1) ) / a(1);
    D = Dt(1);
    m = 1;
    % recalculate Dt and find the minimum, D.
    for i=2:9
        Dt(i) = ( P(i) - T(i) ) / a(i);
        if Dt(i) < D
            D = Dt(i);
            m = i;
        end
    end

end