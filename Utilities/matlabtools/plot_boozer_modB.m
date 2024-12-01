%% plot_boozer_modB( BDATA, NTHETA, NPHI, FILLED, NEWFIG  )
% =============================================
% 
% Plot of modB in Boozer coordinates
%
% INPUT
% -----
%   -bdata    : must be produced by calling read_boozer(filename,root)
%   -Ntheta   : number of meshpoints for theta array
%   -Nphi     : number of meshpoints for phi array
%   -filled   : (=1) pcolor, (=2) contourplot
%   -newfig   : opens (=1) or not (=0) a newfig or overwrites (=2)
%               previous figure
%
% ------------------------------------%
% Written by S.Guinchard (05/18/22)   %
% ------------------------------------%
function plot_boozer_modB(b, Ntheta, Nphi, filled, Newfig)

    modB = get_boozer_modB(b,Ntheta,Nphi);
    Theta = modB.Theta;
    Phi   = modB.Phi;
    modB  = modB.modB;

    switch Newfig

        case 0 
            hold on

        case 1
            figure

        case 2
            hold off

    end

    switch filled

        case 1
            pcolor(Phi, Theta, modB);
            shading interp
            colorbar
            hold on 
            contour(Phi, Theta, modB, 6, 'k', 'linewidth', 1)

        case 2
            
            contour(Phi, Theta, modB,linspace(min(min(modB)), max(max(modB)),20), 'linewidth', 1.5)
            colorbar
            colormap(jet)
    end
end