function plot_spec_pressure(fname, new_figure)

% Plots stepped-pressure profile (without pscale) versus normalized toroidal flux used in SPEC
%   -fname   :      filename in HDF5 format 
%   -new_figure:    1 (0) to (not) open a new figure. =2 to erase existing
%                   figure
%   written by J.Loizu (2018)



data = read_spec_grid(fname);

pvol = data.pressure;

tfl  = data.tflux;

Nvol = data.Nvol;

p0   = zeros(1,10);

if new_figure==1
    figure
    hold on;
elseif new_figure==2
    hold off;
elseif new_figure==0
    hold on;
end


phi_plot = linspace(0, max(tfl), 1E6);
p_plot = zeros(0, 1E6);
temp = 1;

for i=1:length(tfl)
   [val, jj] = min(abs(phi_plot - tfl(i)));
   p_plot(temp:jj) = pvol(i);
   temp = jj;
end

plot(phi_plot, p_plot)

ylabel('p')
xlabel('\Psi / \Psi_{edge}')



