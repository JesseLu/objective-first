
function [source] = backtrack(omega, h, d, kx_cutoff)
    n = length(h);
    kx = fftshift(2*pi*([0:n-1] - n/2) / n);
    ky = sqrt((omega)^2 - kx.^2); % Get ky vectors
    ky = ky .* (1 - 2*(~isreal(ky)));

    H = fft(h); % Get the plane wave components.
    H = H .* exp(i * ky * d); % Backtrack a distance d.
    H = H .* (abs(kx) <= kx_cutoff); % Cutoff.
    H(1)
    source = ifft(H);
