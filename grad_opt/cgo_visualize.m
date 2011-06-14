function cgo_visualize(fval, step_sizes)
% CGO_VISUALIZE(FVAL, STEP_SIZES)
% 
% Description
%     Visualize a c-go optimization.

subplot 211; semilogy(fval, '.-')
ylabel('function value');
xlabel('iterations');

subplot 212; semilogy(step_sizes, '.')
ylabel('step size');
xlabel('iterations');

