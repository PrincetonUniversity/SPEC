function plot_spec_iterations( data, xl, yl )

figure( 'Position', [200 200 900 700], 'Color', 'w' )
fig = gcf;
ax = gca;
ax.Position = [0.1300    0.2386    0.7750    0.6864];


Nfp = double(data.input.physics.Nfp);
Niter = size(data.iterations.iRbc);
Niter = Niter(3);

sld_phi = uicontrol(fig, 'style', 'slider', 'Position', [200    25   500    20],'units','pixel', ...
                'Value',0,'Max',2*pi / Nfp,'Min',0);
sld_it = uicontrol(fig, 'style', 'slider', 'Position', [200 63 500 20],'units','pixel', ...
                'Value',1,'Max', Niter,'Min',1, 'SliderStep', [1/(Niter-1), 1/(Niter-1)]);
            
el_phi = addlistener( sld_phi, 'ContinuousValueChange', @(sld_phi, event) updatePlot(sld_phi, sld_it, data, xl, yl) );
el_it  = addlistener( sld_it , 'ContinuousValueChange', @(sld_it , event) updatePlot(sld_phi, sld_it, data, xl, yl) );

plot_iteration( data, 1, 0, xl, yl )

mytitle = sprintf('%s=%2.3f, iteration %05i', '\phi', 0, 1);
title(mytitle,'FontSize',18)


end




function updatePlot(sld_phi, sld_it, data, xl, yl)
    
    phi = sld_phi.Value;
    iter  = round(sld_it.Value);
    
    plot_iteration( data, iter, phi, xl, yl )
    
    mytitle = sprintf('%s=%2.3f, iteration %05i', '\phi', phi, iter);
    title(mytitle,'FontSize',18)
end


function plot_iteration( data, iter, phi, xl, yl )

    Nfp = double(data.input.physics.Nfp);
    input.representation = 'hudson';
    input.im = double(data.output.im);
    input.in = double(data.output.in) / Nfp;
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
        input.Rmn = data.iterations.iRbc(:,ivol,iter);
        input.Zmn = data.iterations.iZbs(:,ivol,iter);
        
        plot_boundary_2d( input, 512, phi, Nfp, 0 )
    end
    
    
    xlim(xl)
    ylim(yl)
end

