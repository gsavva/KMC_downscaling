Rejected = Rejected + 1;
Decision(Accepted+Rejected+Upscaled) = 3;
% t = f_unscaled(end, 1);
% popul = f_unscaled(end, 3:8);
% t_f is larger than the last sampling time, and population is different
% because 1 more reaction occured after the last sampling point.
t = t_f; 
popul = popul_fun;
sampling_t = sampling_t_f;
% append the last UNscaled trajectory to
% the array of the "accepted" data-points
data(samples+1:samples+samples_f, :) = f_unscaled;
samples = samples + samples_f;
% save occurrence time for slow reactions
for ii=[1 2 7 8 9] % loop over the slow reactions: 1 2 7 8 9
    % retrieve the Non-NaN occurrence times
%     xx = ocT_Un(1:counters_f_UN(ii,N_UN), map(ii));
    xx = ocT_Un(~isnan(ocT_Un(:, map(ii))),map(ii));
    xi = counters(ii,1) + 1;
    xf = counters(ii,1) + size(xx,1); %counters_f_UN(ii,N_UN);
    ocT(xi:xf, map(ii)) = xx;
end
% the KMC steps of the last accepted trajectory are added to the reference counters
counters(:,1) = counters(:,1) + counters_f_UN(1:9,N_UN);

update_T_P
