function threshold(omega, epsilon0, mode_num, exp_factor, t)

epsilon = epsilon0;
for k = 1 : length(t)
    epsilon(3:end-2,:) = 1 * (epsilon0(3:end-2,:) < t(k)) + ...
        12.25 * (epsilon0(3:end-2,:) >= t(k));
    % imagesc(epsilon)
    calc_eff(omega, epsilon, mode_num, exp_factor);
end
