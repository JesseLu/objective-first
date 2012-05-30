function cloak_results()
    dims = [60 100]; % Dimensions of all the couplers.
    eps_init = 9.0; % Initial value of epsilon within design area.
    eps_lims = [1 12.25]; % Limited range of epsilon.
    num_iters = 400; % Number of iterations to run the design optimization for.

%     % Create the specifications for the various design problems.
%     [specs, blocks] = create_specs(dims, eps_lims, eps_init);
% 
%     for k = 1 : length(specs)
%         eps{k} = solve(specs{k}, num_iters, 1e-6, blocks{k});
%     end
% 
%     save('precomp_cloak_results.mat', 'eps', 'specs');

    results = load('precomp_cloak_results.mat', 'eps', 'specs');
    specs = results.specs;
    eps = results.eps;

% Generate the figures.
for k = 1 : length(specs)
    basename = ['fig/cloak', num2str(k), '/'];
    try
        system(['mkdir ', basename]);
    end
    fprintf('\n\nGenerating plots used for result #%d...\n===\n', k);

    fprintf('Simulation results:\n'); % Simulate.
    [eff(k), eps_sim, Ex, Ey, Hz] = simulate(specs{k}, eps{k}, [160 100]);

    % fprintf('\nInput/output mode profiles (figure 1)...\n');
    figure(1); subplot 111; % Generate image files.
    my_area_plot(specs{k}.in.Hz, [basename, 'a']);
    my_area_plot(specs{k}.out.Hz, [basename, 'b']);
    figure(1); % Plot images for user.
    my_plot_images({[basename, 'a.png'], 'Input mode (Hz)'}, ...
                   {[basename, 'b.png'], 'Output mode (Hz)'});
                   
    % fprintf('\nDesign results (figure 2)...\n');
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
        mean(max(max(abs(Hz(1:20,:))))) * [0 1], [basename, 'd']);
    my_imagesc(real(Hz), colormap('jet'), ...
        mean(max(max(abs(real(Hz(1:20,:)))))) * [-1 1], [basename, 'e']);
    figure(2); % Plot images for user.
    my_plot_images({[basename, 'c.png'], 'Relative permittivity'}, ...
                   {[basename, 'd.png'], '|Hz| (simulation)'}, ...
                   {[basename, 'e.png'], 'Re(Hz) (simulation)'});

    fprintf('\nPress enter to continue...\n'); pause; % Wait for user.
end


function [specs, blocks] = create_specs(dims, eps_lims, eps_init)

	% Anti-reflection coating.
	eps0 = ones(dims);
	eps0(3:end,:) = 12.25;
	specs{1} = setup(0.10, eps0, eps_lims, [1 1], 'periodic');
    blocks{1} = false * ones(dims);

	% Wrap cloak. 
	block = false * ones(dims);
	block = my_circle(block, dims/2, 10, 1);
	eps0 = ones(dims);
	eps0(3:end-2,:) = eps_init;
	eps0 = my_circle(eps0, dims/2, 10, 1);
	eps0 = my_circle(eps0, dims/2, 8, -2);
	specs{2} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    blocks{2} = block;
	
	% Front/back cloak.
	block = false * ones(dims);
	block(dims(1)/2 + [-10:10],:) = true;
	eps0 = ones(dims);
	eps0 = my_circle(eps0, dims/2, 8, -2);
	eps0 = block .* eps0 + ~block .* eps_init;
	eps0([1:2, end-1:end],:) = 1;
	specs{3} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    blocks{3} = block;
	
	% Side cloak.
	block = false * ones(dims);
	block(:,dims(2)/2 + [-10:10]) = true;
	eps0 = ones(dims);
	eps0 = my_circle(eps0, dims/2, 8, -2);
	eps0 = block .* eps0 + ~block .* eps_init;
	eps0([1:2, end-1:end],:) = 1;
	specs{4} = setup(0.15, eps0, eps_lims, [1 1], 'periodic');
    blocks{4} = block;

	% Channeler
	block = false * ones(dims);
	block(dims(1)/2 + [-4:4],:) = true;
	eps0 = ones(dims);
	eps0(dims(1)/2 + [-2:2], :) = -20;
	eps0(dims(1)/2 + [-2:2], dims(2)/2 + [-1:1]) = 1;
	eps0 = block .* eps0 + ~block .* eps_init;
	eps0([1:2, end-1:end],:) = 1;
	specs{5} = setup(0.15, eps0, [1 12.25], [1 1], 'periodic');
    blocks{5} = block;




function [z] = my_circle(z, pos, radius, val)
    [x, y] = ndgrid(1:size(z,1), 1:size(z,2));
    r = sqrt((x-pos(1)).^2 + (y-pos(2)).^2);
    in = r <= radius;
    z = in * val + ~in .* z;

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


