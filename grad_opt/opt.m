function [var] = opt(f, g, c, var, varargin)
% VAR = OPT(F, G, C, VAR, STEP_SIZES)
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
% Examples



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

fval0 = f(var);
ss_ind = 1;

for k = 1 : 10
    grad = g(var);
    if (ss_ind >= 2)
        ss_ind = ss_ind - 1;
    end

    while true
        var1 = c(var, grad, step_sizes(ss_ind));
        fval1 = f(var1);
        if (fval1 < fval0)
            fval0 = fval1;
            var = var1;
            break
        else
            ss_ind = ss_ind + 1;
        end

        if (ss_ind > length(step_sizes))
            return
        end
    end
    fval0
end


        
        


