function [var, fval, ss_hist] = opt(f, g, c, var, max_iters, varargin)
% [VAR, FVAL, SS_HIST] = OPT(F, G, C, VAR, MAX_ITERS, STEP_SIZES)
%
% Description
%     Perform constrained, gradient-optimization. Basically, follow the gradient
%     to minimize the objective function while satisfying the constraints.
% 
% Inputs
%     F: Function handle.
%         The objective function to minimize. Must take VAR as input, and return
%         a real scalar as output.
% 
%     G: Function handle.
%         The gradient of F. Must take VAR as input, and return a gradient in 
%         the same "form" (same structure and formatting) as VAR.
% 
%     C: Function handle.
%         Constraining function for the gradient. Takes VAR, the gradient, and 
%         the step size as inputs, and returns the constrained step as output. 
%         Both input VAR, input gradient and output must be formatted as VARs.
% 
%     VAR: Variable.
%         The intial value of the variable to optimize. VAR may be structured
%         in a variety of ways, and input parameters F, G, and C must be 
%         compatible with it.
% 
%     STEP_SIZES: Vector of positive numbers (optional).
%         The different step sizes for OPT to use while optimizing. 
%         Default value is STEP_SIZES = 2.^[0:-1:-40].
%         
% 
% Outputs
%     VAR: Variable.
%         The optimized (technically only improved) value of the variable.
% 
%     FVAL: Vector.
%         The function values at each step, including the initial point.
% 
%     SS_HIST: Vector.
%         The size of the step taken at each iteraton.
%
% Examples
%   See also test1



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
ss_ind(1) = 1; % Start with the largest step size.
k = 1; % Index variable.

% Iterate.
for k = 1 : max_iters

    % Evaluate the gradient at current location.
    grad = g(var);

    % Try step sizes in descending order. 
    % If we can decrease the objective function, then commit to taking the step.
    while true

        % Compute the step and the its function value. 
        var1 = c(var, grad, step_sizes(ss_ind(k)));
        fval1 = f(var1);

        % If function value decreases, take the step.
        if (fval1 < fval(k))
            % Update variables to prepare for next step.
            fval(k+1) = fval1; 
            var = var1; 
            if (ss_ind(k) >= 2) % Increase the next iteration's step size.
                ss_ind(k+1) = ss_ind(k) - 1;
            end
            break % Move on to next iteration.

        % If function value does not decrease, take next smaller step.
        else
            ss_ind(k) = ss_ind(k) + 1;

            % If we just tried smallest step, we are done! 
            if (ss_ind(k) > length(step_sizes))
                ss_hist(k) = 0;
                return
            end
        end
    end

%     % Done with a step, print status.
%     fprintf('%e | %e\n', fval(k), step_sizes(ss_ind(k)));

    % Record the size of the step taken.
    ss_hist(k) = step_sizes(ss_ind(k));
end


        
        


