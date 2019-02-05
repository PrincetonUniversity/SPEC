function output_file = run_from_mhd(Nvol, a, lrad, nptr, r_mhd, p_mhd, mu, iota, psit, psip, L1, debug_figure)
    % Generate a SPEC input file with G=2, L=1 and run it for the given
    % input parameters
    %
    % INPUT
    % -----
    %
    %   Nvol:     number of volumes
    %   lrad:     Radial resolution
    %   nptr:     number of Poincarre points
    %   a:        Minor radiua
    %   r_mhd:    grid for MHD pressure plot
    %   p_mhd:    MHD pressure profile
    %   iota:     MHD rotational transform profile
    %   psit:     Toroidal flux as a function of r_mhd
    %   psip:     Poloidal flux as a function of r_mhd
    %   debug_figure: 1 to see additional figures, 0 otherwise
    %
    % OUTPUT
    % ------
    %   output_file: Name of the htf5 SPEC output file.
    %
    % Written by A.Baillod (2019)



    % =====================================================================
    %                    SPEC INPUT FILE GENERATION
    %                    --------------------------
    
    % Generate (constant) vector for SPEC input files lrad and nptr. 
    lrad = ones(1, Nvol) * lrad;
    nptr = ones(1, Nvol) * nptr;
    
    % Generate phit coordinate
    psit_edge = max(psit); 

    r_spec = linspace(0, a, Nvol+1);
    r_spec = r_spec(2:end);
    psit_spec = interp1(r_mhd(2:end), psit, r_spec, 'spline');
    
    % Interpolate poloidal flux
    psip_spec = interp1(r_mhd(2:end), psip, r_spec, 'spline');
    
    % Then interpolate the rotational transform profile
    iota_spec = interp1(r_mhd, iota, r_spec, 'spline');
    
    % And current profile
    mu = interp1(r_mhd, mu, r_spec, 'spline'); 
    
    if debug_figure
       figure
       plot(r_mhd, iota)
       hold on;
       plot(r_spec, iota_spec,'r*')
       yyaxis right
       plot(r_mhd(2:end), psit)
       plot(r_spec, psit_spec)
       yyaxis left
       for i=1:length(r_spec)
           plot([r_spec(i), r_spec(i)], [0, max(iota)], 'k--')
       end
    end
        
    
    % PRESSURE PROFILE
    % Need the pressure average (integral) over the volume to be equal
    % between the step pressure profile and the mhd profile
    
    % Allocate memory
    p_spec = zeros(1, Nvol);      
    p_spec_plot = zeros(1, length(r_mhd));

    
    % First get last volume
    if Nvol>1
        r_int = linspace(r_spec(Nvol-1), r_spec(Nvol), 100);
        p_int = interp1(r_mhd, p_mhd, r_int, 'spline');
        deltap = trapz(r_int, p_int) ./ (r_spec(Nvol)-r_spec(Nvol-1));
    else
        deltap = 0;
    end
    
    % First volume is special with point r=0
    r_int = linspace(0, r_spec(1), 100);
    p_int = interp1(r_mhd, p_mhd, r_int, 'spline');
    p_spec(1) = trapz(r_int, p_int ./ r_spec(1)) + deltap/(Nvol-1);
    % Some stuff for a debug figure
    if debug_figure
        [temp, jj] = min(abs(r_mhd - r_spec(1)));
        for kk=1:jj
            p_spec_plot(kk) = p_spec(1); 
        end
    end
    
    % Then loop over the other volumes
    for ii=2:Nvol-1
        
        % Integrate the pressure over the volume 
        r_int = linspace(r_spec(ii-1), r_spec(ii), 1E2);
        h = r_spec(ii) - r_spec(ii-1);
        p_int = interp1(r_mhd, p_mhd, r_int, 'spline');
        p_spec(ii) = trapz(r_int, p_int ./ h) + deltap/(Nvol-1);
        
        % For the figure
        if debug_figure
            [temp, jj] = min(abs(r_mhd - r_spec(ii-1)));
            [temp, jj2] = min(abs(r_mhd - r_spec(ii)));
            for kk=jj:jj2
                p_spec_plot(kk) = p_spec(ii); 
            end
        end
    end

    
    % Now generate the file
    inputname = 'G2L1_Nvol'+string(Nvol);
    
    if L1
        write_spec_input_L1('template_G2L1.sp', inputname+'.sp', Nvol,      ...
                            psit_spec, psip_spec, iota_spec, iota_spec, ...
                            psit_edge, p_spec, mu, lrad, nptr)
    else
        write_spec_input_L0('template_G2L0.sp', inputname+'.sp', Nvol, ...
                psit_edge, psit_spec, psip_spec, mu, p_spec, lrad, nptr);
    end
    
    % =====================================================================
    %                          DEBUG FIGURES
    %                          -------------
    
    % Pressure plot
    if debug_figure
        figure
        plot(r_mhd, p_mhd, '-')
        hold on;
        plot(r_spec, p_spec, '*')
        plot(r_spec, ones(1,Nvol), '*')
        plot(r_mhd, p_spec_plot, '-')
    end
    
    % Rotational profile plot
    if debug_figure
       figure
       plot(r_mhd, iota, '-')
       hold on;
       plot(r_spec, iota_spec, 'xk')
       plot(r_spec, ones(1,Nvol), '*')
    end
    
    % =====================================================================
    %                           RUN SPEC
    %                           --------
    command = char('./xspec '+inputname);
    system(command)
    
    % Generate the output filename for alter analysis
    output_file = inputname + '.sp.h5';
end








