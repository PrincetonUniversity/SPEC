function plot_spec_iterations( data, xl, yl )

% PLOT_SPEC_ITERATIONS( DATA, XL, YL )
% ====================================
%
% Plot the volumes interfaces at any iteration of SPEC. Change the
% iteration using a slider
%
% INPUT
% -----
%   -DATA: Obtained from read_spec( filename )
%   -xl  : xlimit
%   -yl  : ylimit
%
% OUTPUT
% ------
%   A magnificient plot!
%

    % Open figure
    figure( 'Position', [200 200 900 700], 'Color', 'w' )
    fig = gcf;
    ax = gca;
    ax.Position = [0.1300    0.2386    0.7750    0.6864];

    % Read size
    Nfp = double(data.input.physics.Nfp);
    [~, ~, Niter] = size(data.iterations.iRbc);
    
    if Niter==0
        error('No iterations available')
    end

    % Create some sliders for controling the iteration and the toroidal
    % angle
    sld_phi = uicontrol(fig, 'style', 'slider', 'Position', [200    25   500    20],'units','pixel', ...
                    'Value',0,'Max',2*pi / Nfp,'Min',0);
    sld_it = uicontrol(fig, 'style', 'slider', 'Position', [200 63 500 20],'units','pixel', ...
                    'Value',1,'Max', Niter,'Min',1, 'SliderStep', [1/(Niter-1), 1/(Niter-1)]);

    addlistener( sld_phi, 'ContinuousValueChange', @(sld_phi, event) updatePlot(sld_phi, sld_it, data, xl, yl) );
    addlistener( sld_it , 'ContinuousValueChange', @(sld_it , event) updatePlot(sld_phi, sld_it, data, xl, yl) );

    % Plot
    plot_iteration( data, 1, 0, xl, yl )
    mytitle = sprintf('%s=%2.3f, iteration %05i', '\phi', 0, 1);
    title(mytitle,'FontSize',18)
end




function updatePlot(sld_phi, sld_it, data, xl, yl)
%
% UPDATEPLOT( SLD_PHI, SLD_IT, DATA, XL, YL )
% ===========================================
%
% Update plot when a slider value has changed
%
% INPUTS
% ------
%   -sld_phi: handle to toroidal angle slider
%   -sld_it : handle to iteration slider
%   -data   : obtained from read_spec(filename)
%   -xl     : xlimit
%   -yl     : ylimit
%
    
    phi = sld_phi.Value;
    iter  = round(sld_it.Value);
    
    plot_iteration( data, iter, phi, xl, yl )
    
    mytitle = sprintf('%s=%2.3f, iteration %05i', '\phi', phi, iter);
    title(mytitle,'FontSize',18)
end


function plot_iteration( data, iter, phi, xl, yl )
%
% PLOT_ITERATION( DATA, ITER, PHI, XL, YL )
% =========================================
%
% Generate plot from SPEC data
%
% INPUTS
% ------
%   -data: obtained from read_spec( filename )
%   -iter: iteration number
%   -phi : toroidal angle
%   -xl  : xlimit
%   -yl  : ylimit
%

    Nfp = double(data.input.physics.Nfp);
    im = double(data.output.im);
    in = double(data.output.in) / Nfp;
    Mvol = double(data.output.Mvol);
    Ntor = double(data.input.physics.Ntor);
    
    % Erase previous plot
    hold off
    
    % Axis
    Ra = 0;
    Za = 0;
    for nn=0:Ntor
       Ra = Ra + data.iterations.iRbc(nn+1,1,iter) * cos( nn*Nfp*phi ); 
       Za = Za - data.iterations.iZbs(nn+1,1,iter) * sin( nn*Nfp*phi );
    end
    
    scatter( Ra, Za, 50, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r' )
    
    hold on
    % volume boundaries
    for ivol=2:Mvol+1
        Rmn = data.iterations.iRbc(:,ivol,iter);
        Zmn = data.iterations.iZbs(:,ivol,iter);
            
        tarr = linspace( 0, 2*pi, 1024 );
        R = zeros(1,1024);
        Z = zeros(1,1024);
        for ii=1:length(Rmn)
            arg = im(ii) * tarr - in(ii) * Nfp * phi;
            R = R + Rmn(ii) * cos( arg );
            Z = Z + Zmn(ii) * sin( arg );
        end
        
        scatter( R, Z )
        hold on
        axis equal
    end
    
    
    xlim(xl)
    ylim(yl)
end

