function [] = paper_results()
% PAPER_RESULTS()
% 
% This is a script which generates the results found in the following paper
% (in submission):
%     Jesse Lu, and Jelena Vuckovic, "Objective-first design of nanophotonic
%     waveguide couplers," (2011)
% 
% A preprint can be found at:
% https://github.com/JesseLu/misc/blob/master/presentations/ob1_wg_paper/paper.pdf
% 
% Running the script in its original form does not actually calculate the 
% results, it simply pulls the pre-computed results from precomp_results.mat. 
% In order to recalculate the results, the first section of this script must
% be uncommented.

dims = [40 80]; % Dimensions of all the couplers.
eps_uniform = 9.0; % Uniform value of epsilon within design area.
eps_lims = [1 12.25]; % Limited range of epsilon.
num_iters = 400;

% Create the specifications for the various design problems.
specs = create_specs(dims, eps_lims, eps_uniform);

% for k = 1 : length(specs)
for k = 4 : 5
    eps{k} = solve(specs{k}, num_iters, 1e-6);
    simulate(specs{k}, eps{k}, [200 160]);
end

save('precomp_results.mat', 'eps', 'specs');
% Note: to make figures pretty, keep only central 160x100 pixels.


function [spec] = create_specs(dims, eps_lims, eps_uniform)

    %
    % Coupler to fundamental mode of wide, low-index waveguide.
    %

w = [10 0.6*dims(2)]; % Widths of the waveguides
eps = ones(dims);
eps(1:2,        (dims(2)-w(1))/2:(dims(2)+w(1))/2) = 12.25; % Input wg.
eps(end-1:end,  (dims(2)-w(2))/2:(dims(2)+w(2))/2) = 2.25; % Output wg.
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{1} = setup(0.15, eps, eps_lims, [1 1]);


    %
    % Coupler from fundamental to second-order mode of silicon waveguide.
    %

w = 16; % Width of the waveguide. Wider, so as to allow for second-order mode.
eps = ones(dims);
eps(:,(dims(2)-w)/2:(dims(2)+w)/2) = 12.25;
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{2} = setup(0.15, eps, eps_lims, [1 2]);
 
 
    %
    % Coupler to an "air-core" fiber mode.
    %

w = 10; % Width of input waveguide
eps = ones(dims);
eps(:,(dims(2)-w)/2:(dims(2)+w)/2) = 12.25; % Input waveguide.

% Create output waveguide.
l = 30; % An effective wavelength.
qw_air = round(l/4); % Quarter wavelength in air.
qw_si = round(qw_air/3.5); % Quarter wavelength in silicon.
eps_temp = [ones(1, qw_air), 12.25*ones(1, qw_si)]; % One period.

% Create quarter-wavelength stack.
epsilon1 = repmat(eps_temp, 1, ceil(dims(2)/2/(qw_air + qw_si))); 
epsilon1 = epsilon1(1:dims(2)/2); % Trim.
epsilon1 = [epsilon1(end:-1:1), epsilon1]; % Mirror, to create cavity.

eps(end-1:end,:)  = repmat(epsilon1, 2, 1); % Stretch out to fill space.
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{3} = setup(0.25, eps, eps_lims, [1 9]);


    %
    % Coupler to a metal-insulator-metal plasmonic waveguide.
    %

w = [10 2]; % Widths of input and output waveguides.
eps = ones(dims);
eps(1:2,        (dims(2)-w(1))/2:(dims(2)+w(1))/2) = 12.25; % Input wg.
eps(end-1:end, :) = -1.1; % Output plasmonic wg.
eps(end-1:end,  (dims(2)-w(2))/2:(dims(2)+w(2))/2) = 1; 
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{4} = setup(0.25, eps, eps_lims, [1 1]);


    %
    % Coupler to a plasmonic wire.
    %

w = [10 2]; % Widths of input and output waveguides.
eps = ones(dims);
eps(1:2,        (dims(2)-w(1))/2:(dims(2)+w(1))/2) = 12.25; % Input wg.
eps(end-1:end,  (dims(2)-w(2))/2:(dims(2)+w(2))/2) = -1.1; % Output wg.
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{5} = setup(0.25, eps, eps_lims, [1 1]);
