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

% 
%     %
%     % Generate the five designs found in the paper.
%     % This section must be uncommented if you want to re-generate the results
%     % yourself.
%     %
% 
% dims = [40 80]; % Dimensions of all the couplers.
% eps_uniform = 9.0; % Uniform value of epsilon within design area.
% eps_lims = [1 12.25]; % Limited range of epsilon.
% num_iters = 400; % Number of iterations to run the design optimization for.
% 
% % Create the specifications for the various design problems.
% specs = create_specs(dims, eps_lims, eps_uniform);
% 
% for k = 1 : length(specs)
%     eps{k} = solve(specs{k}, num_iters, 1e-6);
%     simulate(specs{k}, eps{k}, [200 160]);
% end
% 
% save('precomp_results.mat', 'eps', 'specs');
%

    %
    % Load data and generate figures.
    %

% Load the data.
results = load('precomp_results.mat', 'eps', 'specs');
specs = results.specs;
eps = results.eps;

% Generate the figures.
for k = 1 : length(specs)
    basename = ['fig/res', num2str(k), '/'];
    fprintf('\n\nGenerating plots used for result #%d...\n===\n', k);

    fprintf('Simulation results:\n'); % Simulate.
    [eff(k), eps_sim, Ex, Ey, Hz] = simulate(specs{k}, eps{k}, [160 100]);

    fprintf('\nInput/output mode profiles (figure 1)...\n');
    figure(1); subplot 111; % Generate image files.
    my_area_plot(specs{k}.in.Hz, [basename, 'a']);
    my_area_plot(specs{k}.out.Hz, [basename, 'b']);
    figure(1); % Plot images for user.
    my_plot_images({[basename, 'a.png'], 'Input mode (Hz)'}, ...
                   {[basename, 'b.png'], 'Output mode (Hz)'});
                   
    fprintf('\nDesign results (figure 2)...\n');
    figure(2); subplot 111; % Generate image files.
    if (min(eps{k}(:)) < 0) % Custom colormap if we have metallic devices.
        cmap = flipud([colormap('bone'); flipud(fliplr(colormap('bone')))]);
        r = (abs(min(eps{k}(:))) + 1) / (max(eps{k}(:)) - 1); % Ratio.
        n = size(cmap, 1);
        ind = round(n/2 * (1 - r))+2;
        cmap = interp1(cmap, ind:(n-ind)/64:n);
    else
        cmap = flipud(colormap('bone'));
    end
    my_imagesc(eps_sim, cmap, ...
        [min(eps{k}(:)), max(eps{k}(:))], [basename, 'c']);
    my_imagesc(abs(Hz), colormap('hot'), ...
        mean(max(abs(Hz(1:20,50)))) * [0 1], [basename, 'd']);
    my_imagesc(real(Hz), colormap('jet'), ...
        mean(max(real(Hz(1:20,50)))) * [-1 1], [basename, 'e']);
    figure(2); % Plot images for user.
    my_plot_images({[basename, 'c.png'], 'Relative permittivity'}, ...
                   {[basename, 'd.png'], '|Hz| (simulation)'}, ...
                   {[basename, 'e.png'], 'Re(Hz) (simulation)'});

    fprintf('\nPress enter to continue...\n'); pause; % Wait for user.
end

function my_plot_images(varargin)
% Plot the images for the user to see.
N = length(varargin);
for k = 1 : N
    subplot(1, N, k); 
    [im, map] = imread(varargin{k}{1}); 
    if ~isempty(map) % Try to convert to RGB values if needed (mapped data).
        try
            im = idx2rgb(im, map ); 
        catch 
            colormap('jet');
        end
    end
    image(im);
    title(varargin{k}{2});
    axis equal tight; 
end

function my_imagesc(z, map, lims, filename)
% Write out a mapped image.
z = (((z)-lims(1)) / diff(lims) * 63) + 1;
z = 1 * (z < 1) + 64 * (z > 64) + z .* ((z >= 1) & (z <= 64));
imwrite(z', map, [filename, '.png']);
imwrite([64:-1:1]', map, [filename, '_cbar.png']); % Colorbar.


function my_area_plot(z, filename)
% Normalized area plot.
h = area(1:length(z), z./max(abs(z)));
% set(h, 'FaceColor', [255 194 0]./256); % Tangerine.
axis([1 length(z) -3 3]);
set(gca, 'ytick', []); % No ticks wanted.
print(gcf, '-dpng', '-r150', [filename]); % Save image.
[im] = imread([filename, '.png']); % Reload image.
im = my_add_border(im(300:569,160:1059,:), 0); % Crop and add border.
imwrite(uint8(im), [filename, '.png'], 'png'); % Save.

function [A1] = my_add_border(A0, val)
% Add a one pixel border around image.
A1 = val * ones(size(A0));
A1(2:end-1,2:end-1,:) = A0(2:end-1,2:end-1,:);

function [spec] = create_specs(dims, eps_lims, eps_uniform)
% Create the specifications for all five results.

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
eps(end-1:end, :) = -2; % Output plasmonic wg.
eps(end-1:end,  (dims(2)-w(2))/2:(dims(2)+w(2))/2) = 1; 
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{4} = setup(0.25, eps, eps_lims, [1 1]);


    %
    % Coupler to a plasmonic wire.
    %

w = [10 2]; % Widths of input and output waveguides.
eps = ones(dims);
eps(1:2,        (dims(2)-w(1))/2:(dims(2)+w(1))/2) = 12.25; % Input wg.
eps(end-1:end,  (dims(2)-w(2))/2:(dims(2)+w(2))/2) = -2; % Output wg.
eps(3:end-2, 3:end-2) = eps_uniform; % Fill the design area with uniform eps.

spec{5} = setup(0.25, eps, eps_lims, [1 1]);
