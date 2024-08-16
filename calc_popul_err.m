%% population error
function [error, error_norm] = calc_popul_err(set1, set2)

    error = NaN(4,1);
    for i=1:4 % loop over species: O X X2 [OX2 Yi] Y
        error(i) = norm( set1(:,i) - set2(:,i) );
    end
    % normalize error by number of samples
    % make it independent of sampling_dt
    error_norm = error ./ size(set1,1);

end