function plot_spec_pressure(data, newfig)

% Plots stepped-pressure profile (without pscale) versus normalized toroidal flux used in SPEC
%   -data   : data obtained from read_spec_data
%   -newfig : open a new figue (=1), plots on an existing one (=0) or
%   	       overwrite last plot (=2)
%   written by J.Loizu (2018)


pvol = data.input.physics.pressure * data.input.physics.pscale;
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

p0(1:end) = pvol(1);
tmin      = 0;
tmax      = tfl(1);
tarr      = linspace(tmin,tmax,10);
plot(tarr,p0,'b')
hold on
x         = [tfl(1),tfl(1)];
if Nvol>1
    y         = [pvol(2) pvol(1)];
    plot(x,y,'b')
end

for i=2:Nvol-1

 p0(1:end) = pvol(i);
 tmin      = tfl(i-1);
 tmax      = tfl(i);
 tarr      = linspace(tmin,tmax,10);
 plot(tarr,p0,'b')
 x         = [tfl(i),tfl(i)];
 y         = [pvol(i+1) pvol(i)];
 plot(x,y,'b')
 
end

if Nvol>1
    p0(1:end) = pvol(Nvol);
    tmin      = tfl(Nvol-1);
    tmax      = tfl(Nvol);
    tarr      = linspace(tmin,tmax,10);
    plot(tarr,p0,'b')
    x         = [tfl(Nvol),tfl(Nvol)];
    y         = [0 pvol(Nvol)];
    plot(x,y,'b')
end
    
ylabel('p')
xlabel('\Psi / \Psi_{edge}')

if (pmin~=0 || pmax~=0)
    ylim([0.9*pmin, 1.1*pmax]);
end



