function plot_spec_hessian(data)

%
% PLOT_SPEC_HESSIAN( DATA )
% =========================
%
% Plots the hessian matrix elements
%
% INPUT
% -----
%   -data     : must be produced by calling read_spec_hessian(filename)
%
% written by J.Loizu (2017)

figure

imagesc(data)

colorbar

%set(gca, 'CLim', [-0.1 0.1]);
