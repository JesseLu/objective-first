function my_plot(v, fval, ss_hist)
    
global DIMS_
dims = DIMS_;
N = prod(dims);
figure(1); 
plot_fields(dims, ...
    {'Re(Ey)', real(v.x(N+1:2*N))}, {'|Ey|', abs(v.x(N+1:2*N))}, {'p', v.p});

figure(2); cgo_visualize(fval, ss_hist);

drawnow


