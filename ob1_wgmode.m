function [beta, Hz, Ey] = ob1_wgmode(omega, eps, order, varargin)
% Waveguide mode solver.
% 
% Casts the problem in terms of Hz components, solves for eigenvalues and 
% eigenvectors, and then selects the correct order of eigenvalue.
% Assumes that the mode is propagating in the x-direction.
% 
% If ORDER = 0, then asks user to select one of the plotted modes.
% 
% If user has selected an evanescent mode, then gives a warning.

N = length(eps); 


    % 
    % Setup the eigenvalue problem.
    %

% Get the (smoothed) x and y components of epsilon.
[eps_x, eps_y] = ob1_interp_eps(eps); 

% Build the derivative matrix.
Dy = spdiags([-ones(N,1), ones(N,1)], [-1 0], N, N);
Dy(1,N) = -1;

% Build the matrix for the eigenvalue problem.
my_diag = @(x) spdiags(x(:), 0, N, N);
A = my_diag(eps_y) * (-Dy' * my_diag(eps_x.^-1) * Dy + omega^2 * speye(N));  

% % Remove the wrap-around elements (we don't want periodic boundaries).
% A(N,1) = 0;
% A(1,N) = 0;


   %
   % Solve the eigenvalue problem.
   %

[V, D] = eig(full(A)); % Just use the dense eigenvalue solver.

% Sort modes based on descending beta^2.
[beta2, ind] = sort(real(diag(D)), 'descend');
beta = sqrt(beta2);

% Sort the mode profiles.
Hz = V(:,ind);

% Orient modes so that the largest magnitude element is always positive.
Hz = Hz * diag((max(Hz) > max(-Hz)) - (max(Hz) <= max(-Hz)));

% Calculate Ey.
Ey = ((1 / omega) * my_diag(eps_y(:).^-1) * Hz) * my_diag(beta);

% Normalize power (to unity).
power = diag(Ey' * Hz);
Hz = Hz * diag(power.^-0.5);
Ey = Ey * diag(power.^-0.5);


    %
    % Choose the desired mode.
    %

if (order == 0) % Cycle through modes one at a time and let user choose.
    fprintf('Enter "" (blank) to view next mode, ');
    fprintf('or the number of the mode to be selected.\n');

    for k = 1 : N
        my_mode_plot(Hz(:,k), Ey(:,k)); % Let user see the mode.
        prompt = sprintf('mode %d: ', k); % Tell the user which mode this is.
        choice = input(prompt, 's'); % Ask the user if they want it.

        if strcmp(choice, '') % User wants to see next mode.
            % Delete previous line so we can write on-top of it.
            fprintf(repmat('\b', 1, length(prompt) + length(choice) + 1));
        else % User has chosen a mode
            order = str2num(choice);
            beta = beta(order);
            Hz = Hz(:,order);
            Ey = Ey(:,order);
            break
        end
    end

else % User has already specified what mode they want, so just choose it.
    beta = beta(order);
    Hz = Hz(:,order);
    Ey = Ey(:,order);
end

if ~isreal(beta) % Not a propagating mode, warn user!
    warning('Selected mode is not propagating (non-real wave-vector)!');
end

    %
    % Print out results.
    %

err = norm(A * Hz - beta^2 * Hz);
power = Ey' * Hz;

% Plot result.
if isempty(varargin)
    fprintf('Mode number: %d\n', order);
    fprintf('Beta (wave-vector): %1.3f\n', beta);
    fprintf('Error in eigenvalue equation: %1.1e\n', err);
    fprintf('Power: %1.3f\n\n', power);

    my_mode_plot(Hz, Ey);
elseif ~strcmp(varargin{1}, 'noplot')
    error('Invalid option.');
end



function my_mode_plot(Hz, Ey)
% Visualize both components of the waveguide mode.
subplot 121; my_plot(Hz(:)); title('Hz');
subplot 122; my_plot(Ey(:)); title('Ey');

function my_plot(z)
% Plotting function.
area(1:length(z), z)
axis([1 length(z) -1.5*max(abs(z)) 1.5*max(abs(z))])
