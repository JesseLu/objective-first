function [var, fval, ss_hist] = opt(f, g, c, var, max_iters, varargin)


    %
    % Resolve the step_sizes.
    %

if isempty(varargin)
    step_sizes = 2.^[0 : -1 : -40]; % None given, use default.
else
    step_sizes = sort(varargin{1}(:), 1, 'descend'); % Put in descending order.
end

% Make sure all values are positive.
if any(step_sizes <= 0)
    error('Only positive values of STEP_SIZES are allowed.');
end


    %
    % Perform the optimization.
    %

% Initial conditions.
fval(1) = f(var); % Starting function value.
ss_ind{1}(1) = 1; % Start with the largest step size.
ss_ind{2}(1) = 1; % Start with the largest step size.
k = 1; % Index variable.

% Iterate.
for k = 1 : max_iters

    % Evaluate the gradient at current location.
    grad_full = g(var);

    for l = 1 : 2
        grad = grad_full;
        switch l
            case 1
                grad.phi = 0 * grad.phi;
            case 2
                grad.x = 0 * grad.x;
        end

        % Try step sizes in descending order. 
        % If we can decrease the objective function, 
        % then commit to taking the step.
        while true

            % Compute the step and the its function value. 
            var1 = c(var, grad, step_sizes(ss_ind{l}(k)));
            fval1 = f(var1);

            % If function value decreases, take the step.
            if (fval1 < fval(k))
                % Update variables to prepare for next step.
                fval(k+1) = fval1; 
                var = var1; 
                if (ss_ind{l}(k) >= 2) % Increase the next iteration's step size.
                    ss_ind{l}(k+1) = ss_ind{l}(k) - 1;
                else
                    ss_ind{l}(k+1) = ss_ind{l}(k);
                end
                break % Move on to next iteration.

            % If function value does not decrease, take next smaller step.
            else
                ss_ind{l}(k) = ss_ind{l}(k) + 1;

                % If we just tried smallest step, we are done! 
                if (ss_ind{l}(k) > length(step_sizes))
                    ss_hist{l}(k) = 0;
                    return
                end
            end
        end
        % Record the size of the step taken.
        ss_hist{l}(k) = step_sizes(ss_ind{l}(k));

    end
%     % Done with a step, print status.
%     fprintf('%e | %e\n', fval(k), step_sizes(ss_ind(k)));

end


        
        


