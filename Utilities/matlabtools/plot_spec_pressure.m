function plot_spec_pressure(data, norm, newfig, varargin)

%
% PLOT_SPEC_PRESSURE( DATA, NORM, NEWFIG, VARARGIN )
% ==================================================
%
% Plots stepped-pressure profile versus normalized toroidal flux used in SPEC
%
% INPUT
% -----
%   -data   : data obtained from read_spec(fname)
%   -norm   : (0) plot p, (1) plot p / p0
%   -newfig : open a new figue (=1), plots on an existing one (=0) or overwrite last plot (=2)
%
% written by J.Loizu (2018)
% modified by A. Baillod (2019)

if ~any(norm==[0,1])
    error('InputError: invalid norm')
end

l = length(varargin);
if mod(l,2)~=0
    error('InputError: invalid number of argument')
end

opt.Color='b';
opt.LineWidth=2.0;
for ii=1:l/2
   value = varargin{2*ii};
   field = varargin{2*ii-1};
   
   opt.(field)=value;
end



pvol = data.input.physics.pressure * data.input.physics.pscale;

if norm
    pvol = pvol / pvol(1);
end


pmax = max(pvol);
pmin = min(pvol);

tfl  = data.input.physics.tflux;

Nvol = data.input.physics.Nvol;

p0   = zeros(1,10);

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

% First volume
p0(1:end) = pvol(1);
tmin      = 0;
tmax      = tfl(1);
tarr      = linspace(tmin,tmax,10);
plot(tarr,p0,'Color',opt.Color,'LineWidth',opt.LineWidth)
hold on
x         = [tfl(1),tfl(1)];
if Nvol>1
    y         = [pvol(2) pvol(1)];
    plot(x,y,'Color',opt.Color,'LineWidth',opt.LineWidth)
end

% Next volumes
for i=2:Nvol-1

 p0(1:end) = pvol(i);
 tmin      = tfl(i-1);
 tmax      = tfl(i);
 tarr      = linspace(tmin,tmax,10);
 plot(tarr,p0,'Color',opt.Color,'LineWidth',opt.LineWidth)
 x         = [tfl(i),tfl(i)];
 y         = [pvol(i+1) pvol(i)];
 plot(x,y,'Color',opt.Color,'LineWidth',opt.LineWidth)
 
end

if Nvol>1
    p0(1:end) = pvol(Nvol);
    tmin      = tfl(Nvol-1);
    tmax      = tfl(Nvol);
    tarr      = linspace(tmin,tmax,10);
    plot(tarr,p0,'Color',opt.Color,'LineWidth',opt.LineWidth)
    x         = [tfl(Nvol),tfl(Nvol)];
    y         = [0 pvol(Nvol)];
    plot(x,y,'Color',opt.Color,'LineWidth',opt.LineWidth)
end
    
if( norm )
    ylabel('$p / p_0$', 'Interpreter', 'latex')
else
    ylabel('$p$', 'Interpreter', 'latex')
end
xlabel('$\Psi / \Psi_{edge}$', 'Interpreter', 'latex')

if (pmin~=0 || pmax~=0)
    ylim([0, 1.1*pmax]);
end

set(gcf, 'Color', 'w')
set(gcf,'Position',[200 200 900 700])
set(gca,'FontSize',18)


