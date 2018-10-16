function plot_spec_hessian(data)


% Plots the hessian matrix elements
%   -data     : must be produced by calling read_spec_hessian(filename)
%   written by J.Loizu (2017)

H = data.Hmatrix;

figure

imagesc(H)

colorbar

%set(gca, 'CLim', [-0.1 0.1]);
