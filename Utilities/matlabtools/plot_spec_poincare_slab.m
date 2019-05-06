function plot_spec_poincare_slab(data,nz0,newfig)

% Produces Poincare plots of the field lines on different sections (within one field period)
%
% INPUT
%   -data     : must be produced by calling read_spec_poincare(filename)
%   -nz0      : the toroidal plane to be shown
%   -newfig   : opens(=1) or not(=0) a new figure
%
%   written by J.Loizu (2015)


nptraj = size(data.R_lines,1);  % # of poincare trajectories (field lines)

nz     = size(data.R_lines,2);  % # of toroidal planes

nppts  = size(data.R_lines,3);  % # of toroidal transits per trajectory

try
 rpol   = data.rpol;            % poloidal extent of the slab is 2*pi*rpol
catch
 rpol   = 1;                    % in case rpol did not exist due to old version of SPEC
end

nth    = 1024;  %ploting options for the boundary
bcol   = 'r';
bthick = 3;

R    = squeeze(data.R_lines(:,nz0,:));

T    = mod(squeeze(data.th_lines(:,nz0,:)),2*pi);

if(newfig==1)
figure
end
hold on

for i=1:nptraj       %for each field line trajectory
 scatter(rpol*T(i,:),R(i,:),10,'.k')
 hold on
 set(gca,'FontSize',12)
 xlabel('y','FontSize',12)
 ylabel('R','FontSize',12)
 xlim([0 2*pi*rpol])
 ylim([-0.1 data.Rbc(1,end)+0.1])
end

Rb_u  = 0;
Tb_u  = 0;
Rb_d  = 0;
Tb_d  = 0;
dth   = 2*pi/nth;
theta = dth:dth:2*pi;
zeta  = (nz0-1)*(2*pi/nz);

for imn=1:data.mn     % get and plot the boundary
 alpha = double(data.im(imn))*theta-double(data.in(imn))*zeta;
 Rb_u  = Rb_u + data.Rbc(imn,end)*cos(alpha) + data.Rbs(imn,end)*sin(alpha);
 Rb_d  = Rb_d + data.Rbc(imn,1)*cos(alpha)   + data.Rbs(imn,1)*sin(alpha);
end

scatter(theta*rpol,Rb_u,bthick,'filled',bcol)
scatter(theta*rpol,Rb_d,bthick,'filled',bcol)
