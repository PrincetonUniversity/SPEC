classdef SPEC_Namelist
  properties (Access=public)
		lists
        physicslist
		numericlist
		locallist
		globallist
		diagnosticslist
		screenlist
        initial_guess
  end

  methods (Access=public)
      
    % =====================================================================
    % Class constructor
    function obj = SPEC_Namelist(inputfile)

        % Read namelist
        work = read_namelist( inputfile );
        obj.lists = fields(work);
             
        % Read Lboundary
        if isfield(work.physicslist, 'Lboundary')
            Lboundary = work.physicslist.Lboundary;
        else
            Lboundary = 0; % by default
        end
        
        if Lboundary==0
            % Build im, in arrays 
            
            %If boundary not given
            if ~isfield(work.physicslist, 'Rbc')
               if work.physicslist.Lfreebound==1
                   work.physicslist.Rbc = zeros(size(work.physicslist.Rwc));
                   work.physicslist.Rbs = zeros(size(work.physicslist.Rwc));
                   work.physicslist.Zbc = zeros(size(work.physicslist.Rwc));
                   work.physicslist.Zbs = zeros(size(work.physicslist.Rwc)); 
                   work.physicslist.shift.Rbc = work.physicslist.shift.Rwc;
                   work.physicslist.shift.Rbs = work.physicslist.shift.Rws;
                   work.physicslist.shift.Zbc = work.physicslist.shift.Zwc;
                   work.physicslist.shift.Zbs = work.physicslist.shift.Zws;
               else
                   error('Plasma boundary is not provided')
               end
            end
            
            if ~isfield(work.physicslist, 'Rwc') && work.physicslist.Lfreebound==1
                error('Computational boundary is not provided')
            end
            
            if ~isfield(work.physicslist, 'Vnc') && work.physicslist.Lfreebound==1
                work.physicslist.Vnc = zeros(size(
                        
            s = size(work.physicslist.Rbc);
            mpol =  s(2) - 1;
            ntor = (s(1) - 1) / 2.0;

            im = 0:mpol;
            in = -ntor:ntor;

            [work.physicslist.im, work.physicslist.in] = meshgrid(im,in);  
            
            
            % Need to build Rbc, Zbs, ...
            arr_str = {'Rbc', 'Rbs', 'Zbc', 'Zbs', 'Rwc', 'Rws', 'Zwc', 'Zws', 'Vnc', 'Vns', 'Bnc', 'Bns'};
            n = length(arr_str);
            ntor = zeros(1,n);
            mpol = zeros(1,n);
            for ii=1:n
                
                if ~isfield(work.physicslist, arr_str{ii})
                    ntor(ii) = 0;
                    mpol(ii) = 0;
                    
                    continue
                end
                
                if ~isfield(work.physicslist.shift, arr_str{ii})
                    ntor(ii) = 0;
                    mpol(ii) = 0;
                    
                    continue
                end 
                
                s_shift = size(work.physicslist.(arr_str{ii})) - work.physicslist.shift.(arr_str{ii});
                ntor(ii) = s_shift(1);
                mpol(ii) = s_shift(2);
                
            end
            
            %select poloidal and toroidal resolution for our arrays
            if isfield(work.physicslist, 'Mpol')
                Mpol = max(max(mpol), work.physicslist.Mpol);
                Ntor = max(max(ntor), work.physicslist.Ntor);
            else
                Mpol = max(mpol);
                Ntor = max(ntor);
            end
            
            
            bnd = struct;
            for ii=1:n
               fldname = arr_str{ii};
               bnd.(fldname) = zeros(Mpol+1, 2*Ntor+1);
                               
               if ~isfield(work.physicslist, arr_str{ii}) 
                   continue
               end
                
               if ~isfield(work.physicslist.shift, arr_str{ii})
                   continue
               end
               
               
               flds = size( work.physicslist.(fldname) );
               s = work.physicslist.shift.(fldname);
               for in=1:flds(1)
                       
                   nn = in - s(1);
                   
                   for im=1:flds(2)
                       
                       mm = im - s(2);
                       bnd.(fldname)(mm+1,nn+Ntor+1) = work.physicslist.(fldname)(in,im);
                       
                   end
               end               
            end
            
            % Finally build im, in grids
            im = 0:Mpol;
            in = -Ntor:Ntor;
            
            [in, im] = meshgrid(in, im);
            
            work.physicslist.im = im;
            work.physicslist.in = in;
            
            
            
        else
            
            arr_str = {'R0c', 'Z0s', 'bn','rhomn'};
            n = length(arr_str);
            ntor_arr = zeros(1,n);
            
            for ii=1:n
                if strcmp(arr_str{ii}, 'rhomn')
                    continue
                end
                
                ntor_arr(ii) = length(work.physicslist.(arr_str{ii}))-1;
            end
            
            ntor = max( [work.physicslist.Ntor, max(ntor_arr)] );
            
            for ii=1:n
                if strcmp(arr_str{ii}, 'rhomn')
                    continue
                end
                
                bnd.(arr_str{ii}) = zeros(1, ntor);
                bnd.(arr_str{ii})(1:ntor_arr(ii)+1) = work.physicslist.(arr_str{ii})(1:ntor_arr(ii)+1); 
            end
            
            % Now deal with rhomn
            s = size(work.physicslist.rhomn);
            ntor_rho = max( (s(1)-1)/2.0, work.physicslist.Ntor );
            mpol_rho = max( s(2)-1      , work.physicslist.Mpol );
            
            bnd.rhomn = zeros( mpol_rho+1, 2*ntor_rho+1 );
            
            for ii=1:s(1)
                
                nn = ii - work.physicslist.shift.rhomn(1);
                
                for jj=1:s(2)
                    
                    mm = jj - work.physicslist.shift.rhomn(2);
                    
                    bnd.rhomn( mm+1, nn+ntor_rho+1 ) = work.physicslist.rhomn( ii, jj );
                end
            end
            
            im = 0:mpol_rho;
            in = -ntor_rho:ntor_rho;
            
            [in, im] = meshgrid(in, im);
            
            work.physicslist.im = im;
            work.physicslist.in = in;           
            
        end
        
        % Store in object
        for ii=1:length(obj.lists)
           obj.(obj.lists{ii}) = work.(obj.lists{ii});
        end
        
        for ii=1:n
            fldname = arr_str{ii};
            obj.physicslist.(fldname) = bnd.(fldname);            
        end
        
        % Read initial guess
        try
          obj = obj.read_initial_guess( inputfile );
        catch
          warning('No initial guess available')
        end

    end
      
    function obj = read_initial_guess(obj, init_guess_file)
      %
      % Read geometrical initial guess from another input / .end file
      %
      % INPUT
      % -----
      %   - init_guess_file: filename from which initial guess is read
      %
      % OUTPUT
      % ------
      %   - obj: updated version of SPEC_Namelist object
      %

      fid2 = fopen(init_guess_file, 'r');

      save_line = false;
      category  = '';

      initial_guess = {};

      while( ~feof(fid2) ) % while it is not the end of the file

        tline = fgetl(fid2);
        tline = strtrim(tline);
        
        if( save_line ) % write line
          initial_guess{end+1} = tline;

        else % otherwise look for beginning of geometrical initial guess
          if( strcmp(tline, '/') )
            if( strcmp(category,'screenlist') )
              save_line = true;
            end

          else
            % read category
            if strcmp(tline(end-3:end),'list')
              category = tline(2:end);
            end

          end
        end
      end % end of while

      fclose(fid2); % Close file

      % ==================================
      % Now, format initial guess in arrays
      % Prepare format for reading
      nlines = length(initial_guess);
      
      if isfield(obj.physicslist, 'Lboundary');
          Lboundary = obj.physicslist.Lboundary;
      else
          Lboundary = 0; % by default
      end
      
      if nlines>0

          % Allocate memory
          Mpol = obj.physicslist.Mpol; 
          Ntor = obj.physicslist.Ntor;
          Nvol = obj.physicslist.Nvol;
          Mvol = Nvol + obj.physicslist.Lfreebound;

          if Lboundary==0
              Ric = zeros(Mpol+1,2*Ntor+1,Mvol); Zis = zeros(Mpol+1,2*Ntor+1,Mvol);
              Ris = zeros(Mpol+1,2*Ntor+1,Mvol); Zic = zeros(Mpol+1,2*Ntor+1,Mvol);

              % Read
              for iline=1:nlines

                  % Scan line
                  line_data = str2num( initial_guess{iline} );

                  % Find corresponding index
                  m = line_data(1); n = line_data(2);
                  im = m+1;
                  in = n+Ntor+1;

                  for ivol=1:Nvol
                      Ric(im,in,ivol) = line_data(ivol*4-1);
                      Zis(im,in,ivol) = line_data(ivol*4  );
                      Ris(im,in,ivol) = line_data(ivol*4+1);
                      Zic(im,in,ivol) = line_data(ivol*4+2);
                  end
              end

              init.Ric = Ric;
              init.Ris = Ris;
              init.Zis = Zis;
              init.Zic = Zic;
          else
              
             rhoi = zeros(Mpol  , 2*Ntor+1, Mvol);
             bin  = zeros(Ntor+1, Mvol);
             R0ic = zeros(Ntor+1, Mvol);
             Z0is = zeros(Ntor+1, Mvol);
              
             for iline=1:nlines
                 
                 line_data = str2num( initial_guess{iline} );
                 
                 m = line_data(1);
                 n = line_data(2);
                 
                 for ivol=1:Nvol
                    
                     if n>=0 && m==0
                        bin( n+1, ivol) = line_data(ivol*4-1);
                        R0ic(n+1, ivol) = line_data(ivol*4  );
                        Z0is(n+1, ivol) = line_data(ivol*4+1);
                     end
                     
                     if m>0
                        rhoi( m, n+Ntor+1, ivol ) = line_data( ivol*4+2 );
                     end
                 end
                 
             end
             
             
             init.rhoi = rhoi;
             init.bin  = bin;
             init.R0ic = R0ic;
             init.Z0is = Z0is;
              
          end
      else
          
          if Lboundary==0
              init.Ric = [];
              init.Ris = [];
              init.Zic = [];
              init.Zis = [];
              init.im  = [];
              init.in  = [];
              init.surfaces = cell(0);
          else
              init.rhoi = [];
              init.bin  = [];
              init.R0ic = [];
              init.Z0is = [];
          end
              
      end

      obj.initial_guess = init;
      
      
      if nlines>0 && Lboundary==0
          obj.initial_guess.surfaces = cell(1,Nvol);

          for ii=1:Nvol

              fsurf = fluxSurface( obj.physicslist.Nfp, obj.physicslist.Mpol, ...
                                   obj.physicslist.Ntor, obj.physicslist.Istellsym);

              work = obj.reshape_array_2to1d( squeeze(obj.initial_guess.Ric(:,:,ii)) );     
              fsurf = fsurf.set_array( 'rmnc', work.mode, work.im, work.in );

              work = obj.reshape_array_2to1d( squeeze(obj.initial_guess.Zis(:,:,ii)) );
              fsurf = fsurf.set_array( 'zmns', work.mode, work.im, work.in );

              if ~obj.physicslist.Istellsym
                  work = obj.reshape_array_2to1d( squeeze(obj.initial_guess.Ris(:,:,ii)) );
                  fsurf = fsurf.set_array( 'rmns', work.mode, work.im, work.in );

                  work = obj.reshape_array_2to1d( squeeze(obj.initial_guess.Zic(:,:,ii)) );
                  fsurf = fsurf.set_array( 'zmnc', work.mode, work.im, work.in );
              end

              obj.initial_guess.surfaces{ii} = fsurf;

          end


          % Build im, in arrays
          s = size(obj.initial_guess.Ric);
          mpol =  s(1) - 1;
          ntor = (s(2) - 1) / 2.0;

          im = 0:mpol;
          in = -ntor:ntor;
          [obj.initial_guess.in, obj.initial_guess.im] = meshgrid(in,im);
      end
      

    end
      
    function out = reshape_array_1to2d( obj, arr )
        
        out.Mpol = max(arr.im);
        out.Ntor = max(arr.in);
        out.arr  = zeros(Mpol+1,2*Ntor+1);
        
        nmn = length(in.im);
        
        ii=0;
        for mm=0:out.Mpol
            for nn=-out.Ntor:out.Ntor
                if mm==0 && nn<0
                    continue
                end
                
                ii = ii + 1;
                
                out.arr( mm+1, nn+out.Ntor+1 ) = arr.mode( ii );                
            end
        end
        
        
    end
    
    function out = reshape_array_2to1d( obj, arr )
       
        s = size(arr);
        Mpol =  s(1) - 1;
        Ntor = (s(2) - 1) /2.0 ;
        
        nmn = Ntor + 1 + Mpol*(2*Ntor+1);
        
        out.mode = zeros(1,nmn);
        out.im = zeros(1,nmn);
        out.in = zeros(1,nmn);
        
        ii=0;
        for mm=0:Mpol
            for nn=-Ntor:Ntor
                if mm==0 && nn<0
                    continue
                end
                
                ii = ii+1;
                out.im(ii) = mm;
                out.in(ii) = nn;
                
                out.mode(ii) = arr(mm+1,nn+Ntor+1);
            end
        end
    end
      
    function obj = read_axis_from_focus(obj, focus_data)

      Ntor = obj.physicslist.Ntor;

      axis = focus_data.get_axis_harmonics( Ntor );

      obj.physicslist.Rac = axis.Remn(Ntor+1:end);
      obj.physicslist.Ras = axis.Romn(Ntor+1:end);
      obj.physicslist.Zac = axis.Zemn(Ntor+1:end);
      obj.physicslist.Zas = axis.Zomn(Ntor+1:end);

      end
      
    function obj = interpolate_initial_guess( obj, interp_type, Linitialize )
      %
      % INTERPOLATE_INITIAL_GUESS( INTERP_TYPE, LINITIALIZE )
      % =====================================================
      %
      % Generate an initial guess by interpolation from computational boundary
      % (if free-boundary and Linitialize=2) or from plasma boundary
      % (if Linitialize=1).
      %
      % INPUTS
      % ------
      %   * interp_type: (1) Same interpolation as in SPEC
      %                  (2) Same as in SPEC, but plasma boundary is considered
      %                      as an internal plasma boundary
      %   * Linitialize  (0) No interpolation
      %                  (1) Inner plasma interfaces are interpolated
      %                  (2) Inner plasma and plasma-vacuum interfaces are
      %                      interpolated
      %
      % OUTPUTS
      % -------
      %   * obj: instance of SPEC_Namelist
      %
      % 
      
      error('deprecated')

      machine_precision = eps;
      vsmall = 100 * machine_precision;
      small  = 100 * vsmall;

      Nvol = obj.physicslist.Nvol;
      Mvol = Nvol + obj.physicslist.Lfreebound;

      % First find the axis
      if( obj.physicslist.Igeometry==3 && (obj.physicslist.Rac(1)<small) )

        switch Linitialize
        case -Mvol:-1
          vvol = Nvol + Linitialize;
        case 0
          vvol = 1;
        case 1
          vvol = Nvol;
        case 2
          vvol = Mvol;
        otherwise
          error( 'Invalid Linitialize')
        end

        error(' Interpolation of axis not yet implemented. Please provide one.')
      end

      % Get normalized toroidal flux
      tflux = obj.physicslist.tflux;
      switch interp_type
      case 1
        tflux = tflux / tflux(Nvol);
      case 2
        tflux = tflux / tflux(Mvol);
      otherwise
        error( 'Unknown interp_type input!')
      end

      % Allocate memory
      im   = obj.physicslist.im;
      mn   = length(im);
      psifactor = ones(mn,Mvol);
      inifactor = ones(mn,Mvol);

      % normalization
      if( obj.physicslist.Lfreebound )
        Rscale = obj.physicslist.Rwc(1);
      else
        Rscale = obj.physicslist.Rbc(1);
      end

      % Define psifactor and inifactor
      switch obj.physicslist.Igeometry
      case 1
        psifactor(1:mn,1:Nvol) = one;

      case 2
        for vvol = 1:Nvol
          for ii = 1:mn
            if( im(ii)==0 )
              psifactor(ii,vvol) = tflux(vvol)^(0.5);
            else
              psifactor(ii,vvol) = tflux(vvol)^(0.5*im(ii)-0.5);
            end % of if
          end % of for ii
        end % of for vvol

      case 3
        for vvol = 1:Nvol
          for ii = 1:mn
            if( im(ii)==0 )
              psifactor(ii,vvol) = Rscale;
              inifactor(ii,vvol) = Rscale * tflux(vvol)^0.5;
            else
              psifactor(ii,vvol) = Rscale * tflux(vvol)^(0.5*im(ii));
              inifactor(ii,vvol) = Rscale * tflux(vvol)^(0.5*im(ii));
            end
          end
        end

      otherwise
       error( 'invalid Igeometry for construction of psifactor' )

      end % of switch

      % Define from where we interpolated
      if Linitialize==1
        Rintc = obj.physicslist.Rbc';
        Rints = obj.physicslist.Rbs';
        Zintc = obj.physicslist.Zbc';
        Zints = obj.physicslist.Zbs';

      elseif Linitialize==2
        if obj.physicslist.Lfreebound
          Rintc = obj.physicslist.Rwc';
          Rints = obj.physicslist.Rws';
          Zintc = obj.physicslist.Zwc';
          Zints = obj.physicslist.Zws';
        else
          error( 'Linitialize=2 only valid in free-boundary ')
        end

      else
        disp( 'Linitialize is neither 1 or 2. No interpolation. ' )
        return

      end

      % Now perform interpolation
      Ric = zeros( mn, Mvol ); Ris = zeros( mn, Mvol );
      Zic = zeros( mn, Mvol ); Zis = zeros( mn, Mvol );

      if( Linitialize~=0 )
        switch obj.physicslist.Igeometry
        case 1
          for vvol=1:Nvol
            Ric(:,vvol) = Rintc * tflux(vvol) / tflux(Mvol);
            Ris(:,vvol) = Rints * tflux(vvol) / tflux(Mvol);
          end

        case 2
          if( obj.physicslist.Lfreebound )
            error( 'Cylindrical geometry is not compatible with free-boundary')
          end

          for vvol=1:Nvol-1
            Ric(:,vvol) = Rintc .* psifactor(:, vvol);
            Ris(:,vvol) = Rints .* psifactor(:, vvol);
          end

        case 3
          lvol = Nvol -1 + Linitialize;

          im = obj.physicslist.im;
          mn = length( im );

          Rc_axis = zeros(mn, 1); Rs_axis = zeros(mn, 1);
          Zc_axis = zeros(mn, 1); Zs_axis = zeros(mn, 1);
          counter = 0;
          for imn = 1:mn
            if im(imn)==0
              counter = counter + 1;
              Rc_axis(imn) = obj.physicslist.Rac(counter);
              Rs_axis(imn) = obj.physicslist.Ras(counter);
              Zc_axis(imn) = obj.physicslist.Zac(counter);
              Zs_axis(imn) = obj.physicslist.Zas(counter);
            end
            
            if counter==length(obj.physicslist.Rac)
                break
            end
          end

          for vvol = 1:lvol-1
            Ric(:, vvol) = Rc_axis + (Rintc - Rc_axis) .* (inifactor(:,vvol) / Rscale) ./ tflux(lvol).^(0.5*im');
            Ris(:, vvol) = Rs_axis + (Rints - Rs_axis) .* (inifactor(:,vvol) / Rscale) ./ tflux(lvol).^(0.5*im');
            Zic(:, vvol) = Zc_axis + (Zintc - Zc_axis) .* (inifactor(:,vvol) / Rscale) ./ tflux(lvol).^(0.5*im');
            Zis(:, vvol) = Zs_axis + (Zints - Zs_axis) .* (inifactor(:,vvol) / Rscale) ./ tflux(lvol).^(0.5*im');
          end

        otherwise
          error( 'Invalid geometry')

        end % of switch Igeometry

        if( Linitialize==2 )
          obj.physicslist.Rbc = Ric(:,Nvol);
          obj.physicslist.Rbs = Ris(:,Nvol);
          obj.physicslist.Zbc = Zic(:,Nvol);
          obj.physicslist.Zbs = Zis(:,Nvol);
        end

        obj.initial_guess.Ric = Ric;
        obj.initial_guess.Ris = Ris;
        obj.initial_guess.Zic = Zic;
        obj.initial_guess.Zis = Zis;

      end % of if Linitialize
    end % of function    
    
    
    % =====================================================================
    % Plotters
    function plot_BdotN( obj, Nt, Nphi )
       
        theta = linspace( 0, 2*pi, Nt  );
        phi   = linspace( 0, 2*pi, Nphi);
        
        [thgrid, phgrid] = meshgrid(theta, phi);
        
        im = obj.physicslist.im;
        in = obj.physicslist.in;
        
        nmn = length(im);
        
        BdotN = zeros(size(thgrid));
        Nfp = double(obj.physicslist.Nfp);
        
        for imn=1:nmn
    
            BdotN = BdotN + obj.physicslist.Vnc(imn) * cos(im(imn)*thgrid - in(imn)*Nfp*phgrid) ...
                          + obj.physicslist.Vns(imn) * sin(im(imn)*thgrid - in(imn)*Nfp*phgrid);
            
        end
        
        figure( 'Color', 'w', 'Position', [200 200 900 700] )
        pcolor( thgrid, phgrid, BdotN )
        shading(gca,'interp')
        colorbar
        
    end
    
    function plot_boundary( obj, phi, Nt, Mpol, Ntor, PorW, newfig )
      %
      % PLOT_BOUNDARY( PHI, NT, PorW, NEWFIG )
      % ====================================================
      %
      % Plots spec input initial guess.
      %
      % INPUTS
      % ------
      %   * PHI: Toroidal angle
      %   * NT: Number of points per surface
      %   * Mpol: Poloidal resolution. Set to 0 to use input Mpol
      %   * Ntor: Toroidal resolution. Set to 0 to use input Ntor
      %   * PorW: 'P' Plots the plasma boundary
      %           'W' Plots the computational boundary
      %   * NEWFIG: (0) Plots on gcf
      %             (1) Opens a new figure
      %             (2) Overplots on gcf
      %
      % Written by A. Baillod (2020)
      %
      
      if Mpol==0
          Mpol = obj.physicslist.Mpol;
      end
      
      if Ntor==0
          Ntor = obj.physicslist.Ntor;
      end
      
      im = obj.physicslist.im;
      in = obj.physicslist.in;
      
      if isfield(obj.physicslist,'Lboundary')
        Lboundary = obj.physicslist.Lboundary;
      else
        Lboundary = 0;
      end
      
      
      
      Nfp = obj.physicslist.Nfp;
      theta = linspace(0, 2*pi, Nt );
      R = zeros(1,Nt); 
      Z = zeros(1,Nt);

      
      if Lboundary==0
          switch PorW
              case 'P'
                  Rmnc = obj.physicslist.Rbc;
                  Rmns = obj.physicslist.Rbs;
                  Zmnc = obj.physicslist.Zbc;
                  Zmns = obj.physicslist.Zbs;
              case 'W'
                  Rmnc = obj.physicslist.Rwc;
                  Rmns = obj.physicslist.Rws;
                  Zmnc = obj.physicslist.Zwc;
                  Zmns = obj.physicslist.Zws;
              otherwise
                  error('Invalid input')
          end
          
          s = size(im);

          for ii=1:s(1)
              for jj=1:s(2)

                  if im(ii,jj) > Mpol
                      continue
                  end

                  if abs(in(ii,jj)) > Ntor
                      continue
                  end


                  arg = im(ii,jj) * theta - in(ii,jj) * Nfp * phi;
                  cosarg = cos(arg);
                  sinarg = sin(arg);

                  R = R + Rmnc(ii,jj) * cosarg ;
                  Z = Z + Zmns(ii,jj) * sinarg;

                  if ~obj.physicslist.Istellsym
                      R = R + Rmns(ii,jj) * sinarg ;
                      Z = Z + Zmnc(ii,jj) * cosarg ;
                  end

              end

          end
          
      else
          switch PorW
              case 'P'
                  rhomn = obj.physicslist.rhomn;
                  bn    = obj.physicslist.bn;
                  R0c   = obj.physicslist.R0c;
                  Z0s   = obj.physicslist.Z0s;
              case 'W'
                  error('Not implemented yet')
              otherwise
                  error('Invalid input')
          end
          
          alpha = double(obj.physicslist.twoalpha) / 2.0;
          
          s = size(rhomn);
          Nrho = (s(2)-1)/2.0;
          
          rhoij  = zeros(1, Nt);
          R0 = 0;
          Z0 = 0;
          b  = 0;
          
          lr = length(obj.physicslist.R0c);
          lz = length(obj.physicslist.Z0s);
          lb = length(obj.physicslist.bn);
          
          for nn=0:Ntor
             
              in = nn+1;
              
              if in<=lr
                  R0 = R0 + obj.physicslist.R0c(in) * cos(nn*Nfp*phi );
              end
              if in<=lz
                  Z0 = Z0 + obj.physicslist.Z0s(in) * sin(nn*Nfp*phi );
              end
              if in<=lb
                  b  = b  + obj.physicslist.bn(in)  * cos(nn*Nfp*phi );
              end
          end
          
          zetaij = b*sin(theta - alpha*Nfp*phi);
          
          for ii=1:s(1)
              for jj=1:s(2)
                  mm = ii-1;
                  nn = jj-1-Nrho;
                  
                  arg = mm*theta + nn*Nfp*phi - alpha*Nfp*phi;
                  
                  rhoij = rhoij + obj.physicslist.rhomn(ii,jj)*cos(arg);
              end
          end
          
          R = R0 + rhoij * cos(alpha*Nfp*phi) - zetaij * sin(alpha*Nfp*phi);
          Z = Z0 + rhoij * sin(alpha*Nfp*phi) + zetaij * cos(alpha*Nfp*phi);
          
          
      end
      
      switch newfig
          case 0
              hold on;
          case 1
              figure('Color', 'w', 'Position', [200 200 900 700])
          case 2
              hold off;
      end
      scatter(R,Z,50, 'filled')
      xlabel('$R$[m]', 'Interpreter', 'latex')
      ylabel('$Z$[m]', 'Interpreter', 'latex')
      set(gca, 'FontSize', 18)
        
    end
    
    function plot_initial_guess( obj, phi, Nt, newfig )
      %
      % PLOT_INITIAL_GUESS( PHI, NT, ONLY_BOUNDARY, NEWFIG )
      % ====================================================
      %
      % Plots spec input initial guess.
      %
      % INPUTS
      % ------
      %   * PHI: Toroidal angle
      %   * NT: Number of points per surface
      %   * NEWFIG: (0) Plots on gcf
      %             (1) Opens a new figure
      %             (2) Overplots on gcf
      %
      % Written by A. Baillod (2020)
      %

      % Open figure
      switch newfig
        case 0
          hold on
        case 1
          figure('Position',[200 200 900 700], 'Color','w')
          hold on;
        case 2
          hold off;
        otherwise
          error( 'Unknown value of newfig')
      end

      Nvol   = obj.physicslist.Nvol;

      % Prepare theta coordinate
      theta = linspace(0,2*pi,Nt)';
      
       % Map to real space
       if( ~isempty(obj.initial_guess) )
           newfig = 0;
           for ivol=1:Nvol-1
               obj.initial_guess.surfaces{ivol}.plot_poloidal_section( phi, newfig )
           end
       end
       
       % Plot plasma boundary
       R = zeros(1,Nt);
       Z = zeros(1,Nt);
       
       im = obj.physicslist.im; 
       in = obj.physicslist.in;
       mn = size(im);

       for ii=1:mn(1)
           for jj=1:mn(2)
               cosarg = cos(im(ii,jj)*theta - in(ii,jj)*Np*phi);
               sinarg = sin(im(ii,jj)*theta - in(ii,jj)*Np*phi);

               try
                 R(:) = R(:) + obj.physicslist.Rbc(ii,jj) * cosarg ...
                             + obj.physicslist.Rbs(ii,jj) * sinarg;
                 Z(:) = Z(:) + obj.physicslist.Zbc(ii,jj) * cosarg ...
                             + obj.physicslist.Zbs(ii,jj) * sinarg;
               catch
                 disp('ouch')
               end
           end
       end

       scatter(R, Z, 10, 'filled', 'MarkerFaceColor', [1 0 0])


       if obj.physicslist.Lfreebound
           R = zeros(1,Nt);
           Z = zeros(1,Nt);

           Rwc   = obj.physicslist.Rwc;
           Zwc   = obj.physicslist.Zwc;
           Rws   = obj.physicslist.Rws;
           Zws   = obj.physicslist.Zws;

           for ii=1:mn(1)
               for jj=1:mn(2)
                   cosarg = cos(im(ii,jj)*theta - in(ii,jj)*Np*phi);
                   sinarg = sin(im(ii,jj)*theta - in(ii,jj)*Np*phi);

                   R(:) = R(:) + Rwc(ii,jj) * cosarg + Rws(ii,jj) * sinarg;
                   Z(:) = Z(:) + Zwc(ii,jj) * cosarg + Zws(ii,jj) * sinarg;
               end
           end

           scatter(R, Z, 5, 'filled', 'MarkerFaceColor', [0 0 0])
       end

       %Plot the magnetic axis guess
       Rax = 0; Zax = 0;
       mn = length(obj.physicslist.Rac);
       Nfp = double(obj.physicslist.Nfp);
       Ntor = obj.physicslist.Ntor;
       for n=0:Ntor
        Rax = Rax  + obj.physicslist.Rac(n+1) * cos(n * Nfp * phi) ...
                   - obj.physicslist.Ras(n+1) * sin(n * Nfp * phi);
        Zax = Zax  + obj.physicslist.Zac(n+1) * cos(n * Nfp * phi) ...
                   - obj.physicslist.Zas(n+1) * sin(n * Nfp * phi);
       end

       scatter( Rax, Zax, 50, 'filled', 'MarkerFaceColor', 'g')


       axis equal;

    end
    
    function gif_initial_guess( obj, gifname, Nframes, frame_rate, varargin )

      l = length(varargin);
      if mod(l,2)~=0
          error('Incorrect number of input arguments')
      end
      
      focus = false; f= NaN;
      lines = NaN;
      for ii=1:l/2
         field = varargin{2*ii-1};
         value = varargin{2*ii  };
         
         switch field
             case 'focus'
                 focus = true;
                 f = value;
             case 'lines'
                 lines = value;
         end
      end
        
        
      figdir = pwd;
      GIFNAME = [figdir, '/', gifname ];

      Nfp = double(obj.physicslist.Nfp);
      
      if focus
        s = size(f.outputdata.ppr);
        Nframes = s(2)-1;
      else
        phi_array = linspace(0, 2.0*pi/Nfp, Nframes);
      end
      Nt = 1024;

      % determine xlim and ylim
      ind = find( double(obj.physicslist.im~=0) .* double(obj.physicslist.in~=0) );
      %xmin = obj.physicslist.Rac(1) - sum(abs(obj.physicslist.Rwc(ind))) - sum(abs(obj.physicslist.Rws(ind)));
      %xmax = obj.physicslist.Rac(1) + sum(abs(obj.physicslist.Rwc(ind))) + sum(abs(obj.physicslist.Rws(ind)));
      %ymin = obj.physicslist.Zac(1) - sum(abs(obj.physicslist.Zwc(ind))) - sum(abs(obj.physicslist.Zws(ind)));
      %ymax = obj.physicslist.Zac(1) + sum(abs(obj.physicslist.Zwc(ind))) + sum(abs(obj.physicslist.Zws(ind)));

      %ydiff = ymax-ymin;
      %xdiff = xmax-xmin;

      for iframe = 1:Nframes
          
          if focus
              if isnan(lines)              
                phi = f.plot_poincare( iframe, 1 );
              else
                phi = f.plot_poincare( iframe, 1, 'lines', lines );
              end
                  
              newfig = 0;
          else
              newfig = 1;
              phi = phi_array(iframe);
          end

        plot_initial_guess( obj, phi, Nt, newfig )

        axis equal;

        %xlim([xmin*0.9, xmax*1.1])
        %ylim([ymin*0.9, ymax*1.1])

        drawnow;

        % Capture the plot as an image
        fig = gcf;
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,32);

        % Write to the GIF File
        if iframe == 1
          imwrite(imind,cm,GIFNAME,'gif', 'Loopcount',inf);
        else
          imwrite(imind,cm,GIFNAME,'gif','WriteMode','append', 'DelayTime', frame_rate);
        end

        close( fig);
      end

      disp(' ')
      disp(['Gif saved @ : ',GIFNAME])

    end

    
    % =====================================================================
    % Getters
    function out = get_plasma_boundary(obj, m, n)%
		  % Return Rbc, Zbs, Rbs, Zbc for given m, n mode
	    %
	    % INPUT
	    % -----
	    %   - m: poloidal harmonic
	    %   - n: toroidal harmonic
	    %
	    % OUTPUT
	    % ------
	    %   - out: structure with field out.Rbc, out.Rbs, out.Zbc, out.Zbs
	    %

          Mpol = max(max(obj.physicslist.im));
          Ntor = max(max(obj.physicslist.in));
        
		  if m>Mpol || abs(n)>Ntor
			  out.Rbc = 0.0;
			  out.Rbs = 0.0;
			  out.Zbc = 0.0;
			  out.Zbs = 0.0;
		  else
			  out.Rbc = obj.physicslist.Rbc(m+1, n+Ntor+1);
			  out.Rbs = obj.physicslist.Rbs(m+1, n+Ntor+1);
			  out.Zbc = obj.physicslist.Zbc(m+1, n+Ntor+1);
			  out.Zbs = obj.physicslist.Zbs(m+1, n+Ntor+1);
		  end

    end

    function out = get_wall(obj, m, n)
      %
		  % Return Rwc, Zws, Rws, Zwc for given m, n mode
	    %
	    % INPUT
	    % -----
	    %   - m: poloidal harmonic
	    %   - n: toroidal harmonic
	    %
	    % OUTPUT
	    % ------
	    %   - out: structure with field out.Rwc, out.Rws, out.Zwc, out.Zws
	    %
          Mpol = max(max(obj.physicslist.im));
          Ntor = max(max(obj.physicslist.in));

		  if m>Mpol || abs(n)>Ntor
			  out.Rwc = 0.0;
			  out.Rws = 0.0;
			  out.Zwc = 0.0;
			  out.Zws = 0.0;
		  else
			  out.Rwc = obj.physicslist.Rwc(m+1, n+Ntor+1);
			  out.Rws = obj.physicslist.Rws(m+1, n+Ntor+1);
			  out.Zwc = obj.physicslist.Zwc(m+1, n+Ntor+1);
			  out.Zws = obj.physicslist.Zws(m+1, n+Ntor+1);
		  end
    end

    function out = get_vnc_bns(obj, m, n)
      % Return Vnc, Vns, Bnc, Bns for given m, n mode
	    %
	    % INPUT
	    % -----
	    %   - m: poloidal harmonic
	    %   - n: toroidal harmonic
	    %
	    % OUTPUT
	    % ------
	    %   - out: structure with field out.Vnc, out.Vns, out.Bnc, out.Bns
	    %
          Mpol = max(max(obj.physicslist.im));
          Ntor = max(max(obj.physicslist.in));

		  if m>Mpol || abs(n)>Ntor
			  out.Vnc = 0.0;
			  out.Vns = 0.0;
			  out.Bnc = 0.0;
			  out.Bns = 0.0;
		  else
			  out.Vnc = obj.physicslist.Vnc(m+1, n+Ntor+1);
			  out.Vns = obj.physicslist.Vns(m+1, n+Ntor+1);
			  out.Bnc = obj.physicslist.Bnc(m+1, n+Ntor+1);
			  out.Bns = obj.physicslist.Bns(m+1, n+Ntor+1);
		  end
    end
    
    function out = get_fourier_harmonics(obj, m, n, field)
        
      Mpol = max(max(obj.physicslist.im));
      Ntor = max(max(obj.physicslist.in));
        

      if( strcmp(field, 'Ric') || strcmp(field, 'Ris') || strcmp(field, 'Zic') || strcmp(field, 'Zis'))
        if( isempty(obj.initial_guess) )
          error('No initial guess')
        else
          out = obj.initial_guess.(field)(ii);
        end
      else
        out = obj.physicslist.(field)(m+1, n+Ntor+1);
      end
    end
    

	% =====================================================================
    % Setters
    function obj = truncate_fourier_series( obj, min_value )
        
       nmn = length(obj.physicslist.im);
       
       ind = cell(1,12);
       
       ind{1} = find(obj.physicslist.Rbc<min_value);
       ind{2} = find(obj.physicslist.Rbs<min_value);
       ind{3} = find(obj.physicslist.Zbc<min_value);
       ind{4} = find(obj.physicslist.Zbs<min_value);
       
       ind{5} = find(obj.physicslist.Rwc<min_value);
       ind{6} = find(obj.physicslist.Rws<min_value);
       ind{7} = find(obj.physicslist.Zwc<min_value);
       ind{8} = find(obj.physicslist.Zws<min_value);
       
       ind{9} = find(obj.physicslist.Vnc<min_value);
       ind{10} = find(obj.physicslist.Vns<min_value);
       ind{11} = find(obj.physicslist.Bnc<min_value);
       ind{12} = find(obj.physicslist.Bns<min_value);
       
       common_index = 1:nmn;
       for ii=1:12
          common_index = intersect( common_index, ind{ii} );
       end
       
       new_im = obj.physicslist.im;
       new_im(common_index) = [];
       obj.physicslist.Mpol = max(new_im);
       
       new_in = obj.physicslist.in;
       new_in(common_index) = [];
       obj.physicslist.Ntor = max(abs(new_in));
       
       ind1 = find( obj.physicslist.im > obj.physicslist.Mpol );
       ind2 = find( abs(obj.physicslist.in) > obj.physicslist.Ntor );
       
       ind = union( ind1, ind2 );
        
       obj.physicslist.im(ind) = [];
       obj.physicslist.in(ind) = [];
       obj.physicslist.Rbc(ind) = [];
       obj.physicslist.Rbs(ind) = [];
       obj.physicslist.Zbc(ind) = [];
       obj.physicslist.Zbs(ind) = [];
       
       obj.physicslist.Rwc(ind) = [];
       obj.physicslist.Rws(ind) = [];
       obj.physicslist.Zwc(ind) = [];
       obj.physicslist.Zws(ind) = [];
       
       obj.physicslist.Vnc(ind) = [];
       obj.physicslist.Vns(ind) = [];
       obj.physicslist.Bnc(ind) = [];
       obj.physicslist.Bns(ind) = []; 
       
       if ~isempty( obj.initial_guess )
           obj.initial_guess.Ric(ind,:) = [];
           obj.initial_guess.Ris(ind,:) = [];
           obj.initial_guess.Zic(ind,:) = [];
           obj.initial_guess.Zis(ind,:) = [];
       end
    end
    
    function obj = set_fourier_harmonics(obj, im, in, values, field, varargin)

      l = length(varargin);
      if( mod(l,2)~=0 )
        error('Invalid number of argument')
      end

      set_to_zero=true;
      for iv=1:l/2
        fieldarg=varargin{2*iv-1};
        value=varargin{2*iv};
        switch fieldarg
	  case 'set_to_zero'
	    set_to_zero=value;
	  otherwise
	    error('Unknown input')
        end
      end


      mn = length(im);
      mn_old = length(obj.physicslist.Rbc);
      
      
      if set_to_zero
        obj.physicslist.(field) = zeros(size(obj.physicslist.Rbc));
      end

      for ii = 1:mn
        m = im(ii); n=in(ii);
        imn = obj.find_ii_given_mn(m, n);

        if( isempty(imn) ) % Add mode
          obj.physicslist.im(end+1) = m;
          obj.physicslist.in(end+1) = n;

          obj.physicslist.Rbc(end+1) = 0.0;
          obj.physicslist.Rbs(end+1) = 0.0;
          obj.physicslist.Zbc(end+1) = 0.0;
          obj.physicslist.Zbs(end+1) = 0.0;

          obj.physicslist.Rwc(end+1) = 0.0;
          obj.physicslist.Rws(end+1) = 0.0;
          obj.physicslist.Zwc(end+1) = 0.0;
          obj.physicslist.Zws(end+1) = 0.0;

          obj.physicslist.Vnc(end+1) = 0.0;
          obj.physicslist.Vns(end+1) = 0.0;
          obj.physicslist.Bnc(end+1) = 0.0;
          obj.physicslist.Bns(end+1) = 0.0;

          if( ~isempty(obj.initial_guess) )
            obj.initial_guess.Ric(end+1, :) = 0.0;
            obj.initial_guess.Ris(end+1, :) = 0.0;
            obj.initial_guess.Zic(end+1, :) = 0.0;
            obj.initial_guess.Zis(end+1, :) = 0.0;
          end

          obj.physicslist.(field)(end) = values(ii);

        else
          obj.physicslist.(field)(imn) = values(ii);

        end
        

      end
      
      if ~(strcmp(field, 'Vnc') || strcmp(field, 'Vns') || strcmp(field, 'Bnc') || strcmp(field, 'Bns') )
        obj = obj.update_flux_surfaces();
      end

    end
    
    function obj = convert_to_standard_representation( obj )
        
       Lb = obj.physicslist.Lboundary;
       if Lb==0
           error('Can only convert from Lboundary=1 to Lboundary=0')
       end
       
       obj.physicslist.Lboundary=0;
       obj.physicslist.twoalpha=0;
       
       % Map computational boundary
       if( obj.physicslist.Lfreebound ) 
          error('Not implemented for freeboundary cases yet')
       end
       
       % Map computational boundary
       rhomn = obj.physicslist.rhomn;
       bn = obj.physicslist.bn;
       R0c = obj.physicslist.R0c;
       Z0s = obj.physicslist.Z0s;
       alpha = obj.physicslist.twoalpha / 2.0;
       
       [Rmn, Zmn, Mpol, Ntor] = map_rho_to_RZ( rhomn, bn, R0c, Z0s, alpha, 0);
       
       obj.physicslist.Mpol = Mpol;
       obj.physicslist.Ntor = Ntor;
       
       obj.physicslist.Rbc = Rmn;
       obj.physicslist.Zbs = Zmn;
       obj.physicslist.Rbs = zeros(size(Rmn));
       obj.physicslist.Zbc = zeros(size(Zmn));
       
       obj.physicslist = rmfield(obj.physicslist, 'rhomn');
       obj.physicslist = rmfield(obj.physicslist, 'bn');
       obj.physicslist = rmfield(obj.physicslist, 'R0c');
       obj.physicslist = rmfield(obj.physicslist, 'Z0s');
       
       % Map initial guess
       if ~isempty(obj.initial_guess )
           Nvol = obj.physicslist.Nvol;
           
           obj.initial_guess.Ric = zeros( Mpol+1, 2*Ntor+1, Nvol );
           obj.initial_guess.Ris = zeros( Mpol+1, 2*Ntor+1, Nvol );
           obj.initial_guess.Zic = zeros( Mpol+1, 2*Ntor+1, Nvol );
           obj.initial_guess.Zis = zeros( Mpol+1, 2*Ntor+1, Nvol );
           
           % Plasma boundary
           obj.initial_guess.Ric(:,:,Nvol) = Rmn;
           obj.initial_guess.Zis(:,:,Nvol) = Zmn;
           
           % Inner interfaces
           for ivol=1:Nvol-1
               [Rmn, Zmn, Mpol, Ntor] = map_rho_to_RZ( obj.initial_guess.rhoi(:,:,ivol), ...
                                                       obj.initial_guess.bin(:,ivol), ...
                                                       obj.initial_guess.R0ic(:,ivol), ...
                                                       obj.initial_guess.Z0is(:,ivol), ...
                                                       alpha, 0);
               obj.initial_guess.Ric(1:Mpol+1,1:2*Ntor+1,ivol) = Rmn;
               obj.initial_guess.Zis(1:Mpol+1,1:2*Ntor+1,ivol) = Zmn;
           end
           
           obj.initial_guess = rmfield(obj.initial_guess, 'rhoi');
           obj.initial_guess = rmfield(obj.initial_guess, 'bin');
           obj.initial_guess = rmfield(obj.initial_guess, 'R0ic');
           obj.initial_guess = rmfield(obj.initial_guess, 'Z0is');
       end
        
        
    end
    
    
    
    function obj = set_initial_guess_harmonics( obj, volume, im, in, values, field)
        
        
      mn = length(im);
      mn_old = length(obj.physicslist.Rbc);
      Mvol = obj.physicslist.Nvol + obj.physicslist.Lfreebound;
      
      if isempty(obj.initial_guess)
         obj.initial_guess.Ric = zeros( mn_old, Mvol );
         obj.initial_guess.Ris = zeros( mn_old, Mvol );
         obj.initial_guess.Zic = zeros( mn_old, Mvol );
         obj.initial_guess.Zis = zeros( mn_old, Mvol );
      end

      %obj.initial_guess.(field) = zeros(mn_old, Mvol);
      tmp = obj.initial_guess.(field);
      tmp( 1:mn_old, volume ) = zeros( mn_old, 1 );

      for ii = 1:mn
        m = im(ii); n=in(ii);
        imn = obj.find_ii_given_mn(m, n);

        if( isempty(imn) ) % Add mode
          obj.physicslist.im(end+1) = m;
          obj.physicslist.in(end+1) = n;

          obj.physicslist.Rbc(end+1) = 0.0;
          obj.physicslist.Rbs(end+1) = 0.0;
          obj.physicslist.Zbc(end+1) = 0.0;
          obj.physicslist.Zbs(end+1) = 0.0;

          obj.physicslist.Rwc(end+1) = 0.0;
          obj.physicslist.Rws(end+1) = 0.0;
          obj.physicslist.Zwc(end+1) = 0.0;
          obj.physicslist.Zws(end+1) = 0.0;

          obj.physicslist.Vnc(end+1) = 0.0;
          obj.physicslist.Vns(end+1) = 0.0;
          obj.physicslist.Bnc(end+1) = 0.0;
          obj.physicslist.Bns(end+1) = 0.0;

          obj.initial_guess.Ric(end+1, :) = 0.0;
          obj.initial_guess.Ris(end+1, :) = 0.0;
          obj.initial_guess.Zic(end+1, :) = 0.0;
          obj.initial_guess.Zis(end+1, :) = 0.0;

          tmp( end+1, : ) = 0.0;
          tmp( end  , volume ) = values(ii);

        else
            tmp( imn, volume ) = values(ii);

        end
        
      obj.initial_guess.(field) = tmp;
        

      end
        
        
    end  
    
    
    % =====================================================================
    % Utility
    function compare_fourier_harmonics(obj, spec_nm, field)

      nm_in = length(obj.physicslist.im);
      nm_ou = length(spec_nm.physicslist.im);

      nm = max([nm_in, nm_ou]);

      for ii = 1:nm

        if( nm_in>nm_ou)
          m = obj.physicslist.im(ii);
          n = obj.physicslist.in(ii);
        else
          m = spec_nm.physicslist.im(ii);
          n = spec_nm.physicslist.in(ii);
        end
        
        ind_in = obj.find_ii_given_mn( m, n);
        ind_ou = spec_nm.find_ii_given_mn( m, n);

        if( isempty(ind_in) )
          if( ~isempty(ind_ou) )
            value = spec_nm.physicslist.(field)(ind_ou);
            if( value ~= 0 )
              disp( [field, '(', num2str(m), ',', num2str(n) ') missmatch. Obj has value 0 while spec_nm has value ', num2str(value)]);
            end
          end
        else
          val_in = obj.physicslist.(field)(ind_in);
          if( isempty(ind_ou) )
            disp( [field, '(', num2str(m), ',', num2str(n) ') missmatch. Obj has value ', num2str(val_in), ' while spec_nm has value 0']);
          else
            val_ou = spec_nm.physicslist.(field)(ind_ou);

            delta = abs(val_in - val_ou);
            if( delta>1E-16 )
              disp( [field, '(', num2str(m), ',', num2str(n) ') missmatch. Obj has value ', num2str(val_in), ' while spec_nm has value ', num2str(val_ou)])
            end
          end
        end
        


      end




    end
    
    function obj = add_volume( obj, tflux )
        
        Nvol = obj.physicslist.Nvol;
        if Nvol==1
            error('Not implemented for Nvol=1')
        end
       
        obj.physicslist.Nvol = obj.physicslist.Nvol + 1;
        
        % Find between which interface to put it
        [tmp, ind] = min(abs(tflux - obj.physicslist.tflux));
        tmp = obj.physicslist.tflux(ind);
        
        if tmp-tflux>0
            ind = ind-1;
        end
        
        disp(newline)
        disp(newline)
        fprintf('Adding a volume between interface %u and %u ...\n', ind, ind+1)
        fprintf('======================================================\n')
        
        % Build new profiles
        old_tflux = obj.physicslist.tflux;
        obj.physicslist.tflux = [obj.physicslist.tflux(1:ind), tflux, obj.physicslist.tflux(ind+1:end)];
        
        pflux = obj.physicslist.pflux;
        if ind==0
            newpflux = pflux(ind+1)/2;
        else
            newpflux = (pflux(ind)+pflux(ind+1))/2;
        end
        obj.physicslist.pflux = [pflux(1:ind), newpflux, pflux(ind+1:end)];
        
        Lrad = obj.physicslist.Lrad;
        obj.physicslist.Lrad = [Lrad(1:ind), Lrad(ind+1), Lrad];
        
        hel = obj.physicslist.helicity;
        obj.physicslist.helicity = [hel(1:ind), 0, hel(ind+1:end)];
        
        pressure = obj.physicslist.pressure;
        if ind==0
            newpressure = pressure(ind+1)*1.3;
        else
            newpressure = (pressure(ind)+pressure(ind+1))/2;
        end
        obj.physicslist.pressure = [pressure(1:ind), newpressure, pressure(ind+1:end)];
        
        adiabatic = obj.physicslist.adiabatic;
        if ind==0
            newadiab = adiabatic(ind+1);
        else
            newadiab = (adiabatic(ind)+adiabatic(ind+1))/2;
        end
        obj.physicslist.adiabatic = [adiabatic(1:ind), newadiab, adiabatic(ind+1:end)];
        
        mu = obj.physicslist.mu;
        obj.physicslist.mu = [mu(1:ind), mu(ind+1), mu(ind+1:end)];
        
        Ivolume = obj.physicslist.Ivolume;
        obj.physicslist.Ivolume = [Ivolume(1:ind), Ivolume(ind+1), Ivolume(ind+1:end)];
        
        Isurf = obj.physicslist.Isurf;
        if ind==0
            newIsurf = Isurf(1);
        else
            newIsurf = Isurf(ind);
        end
        obj.physicslist.Isurf = [Isurf(1:ind), newIsurf, Isurf(ind+1:end)];
        
        iota = obj.physicslist.iota;
        if ind==0
            newiota = iota(ind+1);
        else
            newiota = (iota(ind)+iota(ind+1))/2;
        end
        obj.physicslist.iota = [iota(1:ind), newiota, iota(ind+1:end)];
        
        oita = obj.physicslist.oita;
        if ind==0
            newoita = oita(ind+1);
        else
            newoita = (oita(ind)+oita(ind+1))/2;
        end
        obj.physicslist.oita = [oita(1:ind), newoita, oita(ind+1:end)];
        
        disp(newline)
        fprintf('New profiles have been generated with averaged values. Please modify if necessary.\n')
        fprintf(' tflux     = %.8E \n', tflux )
        fprintf(' pflux     = %.8E \n', newpflux )
        fprintf(' Lrad      = %.8E \n', Lrad(ind+1) )
        fprintf(' helicity  = %.8E \n', 0 )
        fprintf(' pressure  = %.8E \n', newpressure )
        fprintf(' adiabatic = %.8E \n', newadiab )
        fprintf(' mu        = %.8E \n', mu(ind+1) )
        fprintf(' Ivolume   = %.8E \n', Ivolume(ind+1) )
        fprintf(' Isurf     = %.8E \n', newIsurf )
        fprintf(' iota      = %.8E \n', newiota )
        fprintf(' oita      = %.8E \n', newoita )
        disp(newline)
        
        
        % Interpolate initial guess
        if( isempty(obj.initial_guess) )
            disp('No initial guess to build, returning ...')
            disp(newline)
            return
        end
        
        
        
        % First select closest boundaries harmonics 
%         
%         nmn = length( obj.physicslist.im );
%         Mpol = max(obj.physicslist.im);
%         Ntor = max(obj.physicslist.in);
%         stellsym = false;
%         
%         newGuess = fluxSurface( obj.physicslist.Nfp, Mpol, Ntor, stellsym );
%         
%         for imn = 1:nmn
%            
%             m = obj.physicslist.im(imn);
%             n = obj.physicslist.in(imn);
%             
%             x = [0, old_tflux( 1:end-1 )];
%             
%             y = obj.initial_guess.Ric( imn, 1:end-1 );
%             y = [0, y];
%             if m==0
%                y(1) = obj.physicslist.Rac( n+1 );
%             end
%             p = polyfit( x, y./x, 1 );
%             value = polyval( p, tflux ) * tflux;
%             newGuess = newGuess.set_value( 'rmnc', value, m, n );
%             
%             y = obj.initial_guess.Ris( imn, 1:end-1 );
%             y = [0, y];
%             if m==0
%                y(1) = obj.physicslist.Ras( n+1 );
%             end
%             p = polyfit( x, y./x, 1 );
%             value = polyval( p, tflux ) * tflux;
%             newGuess = newGuess.set_value( 'rmns', value, m, n );
%             
%             y = obj.initial_guess.Zic( imn, 1:end-1 );
%             y = [0, y];
%             if m==0
%                y(1) = obj.physicslist.Zac( n+1 );
%             end
%             p = polyfit( x, y./x, 1 );
%             value = polyval( p, tflux ) * tflux;
%             newGuess = newGuess.set_value( 'zmnc', value, m, n );
%             
%             y = obj.initial_guess.Zis( imn, 1:end-1 );
%             y = [0, y];
%             if m==0
%                y(1) = obj.physicslist.Zas( n+1 );
%             end
%             p = polyfit( x, y./x, 1 );
%             value = polyval( p, tflux ) * tflux;
%             newGuess = newGuess.set_value( 'zmns', value, m, n );
%         end
%         
%         obj.initial_guess.Ric = [obj.initial_guess.Ric(:,1:ind), newGuess.rmnc', obj.initial_guess.Ric(:,ind+1:end)];
%         obj.initial_guess.Ris = [obj.initial_guess.Ris(:,1:ind), newGuess.rmns', obj.initial_guess.Ris(:,ind+1:end)];
%         obj.initial_guess.Zic = [obj.initial_guess.Zic(:,1:ind), newGuess.zmnc', obj.initial_guess.Zic(:,ind+1:end)];
%         obj.initial_guess.Zis = [obj.initial_guess.Zis(:,1:ind), newGuess.zmns', obj.initial_guess.Zis(:,ind+1:end)];
        
        
        
        nmn = length(obj.physicslist.im);
        
        if ind==0 %Then first volume, need regularisation
            l = length(obj.physicslist.Rac);
            if l>obj.physicslist.Ntor+1
                l = obj.physicslist.Ntor+1;
            end
            
            Ric_low = zeros(nmn, 1);
            Ric_low(1:l) = obj.physicslist.Rac(1:l);
            Ris_low = zeros(nmn, 1);
            Ris_low(1:l) = obj.physicslist.Ras(1:l);
            Zic_low = zeros(nmn, 1);
            Zic_low(1:l) = obj.physicslist.Zac(1:l);
            Zis_low = zeros(nmn, 1);
            Zis_low(1:l) = obj.physicslist.Zas(1:l);
            
            delta_tflux = obj.physicslist.tflux(ind+2);
            tflux_norm = tflux / delta_tflux;
            
        else
            delta_tflux = obj.physicslist.tflux(ind+2) - obj.physicslist.tflux(ind);
            tflux_norm = (tflux - obj.physicslist.tflux(ind)) / delta_tflux;
            Ric_low = obj.initial_guess.Ric(:,ind);
            Ris_low = obj.initial_guess.Ris(:,ind);
            Zic_low = obj.initial_guess.Zic(:,ind);
            Zis_low = obj.initial_guess.Zis(:,ind);
            
        end
        
        
        if ind==Nvol-1
            Ric_high = obj.physicslist.Rbc';
            Ris_high = obj.physicslist.Rbs';
            Zic_high = obj.physicslist.Zbc';
            Zis_high = obj.physicslist.Zbs';
            
        else
            Ric_high = obj.initial_guess.Ric(:,ind+1);
            Ris_high = obj.initial_guess.Ris(:,ind+1);
            Zic_high = obj.initial_guess.Zic(:,ind+1);
            Zis_high = obj.initial_guess.Zis(:,ind+1);
            
        end
        
        % Now build regularisation factor
        fj = zeros(nmn, 1);
        
        if ind==0 % first volume, need to regularize non-zero m-harmonics
            
            for ii=1:nmn
               im = obj.physicslist.im(ii);
               
               if im==0
                   fj(ii) = sqrt(tflux_norm);
               else
                   fj(ii) = tflux_norm.^(im / 2.0);
               end
            end
            
        else % Not in first volume, them linear interpolation in radius
            
            fj = ones(nmn, 1) * sqrt(tflux_norm);
            
        end
        
        
        % Interpolate
        new_Ric = Ric_low + (Ric_high - Ric_low) .* tflux_norm;
        new_Ris = Ris_low + (Ris_high - Ris_low) .* tflux_norm;
        new_Zic = Zic_low + (Zic_high - Zic_low) .* tflux_norm;
        new_Zis = Zis_low + (Zis_high - Zis_low) .* tflux_norm;
        
        
        obj.initial_guess.Ric = [obj.initial_guess.Ric(:,1:ind), new_Ric, obj.initial_guess.Ric(:,ind+1:end)];
        obj.initial_guess.Ris = [obj.initial_guess.Ris(:,1:ind), new_Ris, obj.initial_guess.Ris(:,ind+1:end)];
        obj.initial_guess.Zic = [obj.initial_guess.Zic(:,1:ind), new_Zic, obj.initial_guess.Zic(:,ind+1:end)];
        obj.initial_guess.Zis = [obj.initial_guess.Zis(:,1:ind), new_Zis, obj.initial_guess.Zis(:,ind+1:end)];
        
       
    end
    
    function obj = remove_volume( obj, ivol )
        
       obj.physicslist.Nvol      = obj.physicslist.Nvol-1;
       obj.physicslist.Lrad      = [obj.physicslist.Lrad(1:ivol-1),      obj.physicslist.Lrad(ivol+1:end)     ];
       obj.physicslist.tflux     = [obj.physicslist.tflux(1:ivol-1),     obj.physicslist.tflux(ivol+1:end)    ];
       obj.physicslist.pflux     = [obj.physicslist.pflux(1:ivol-1),     obj.physicslist.pflux(ivol+1:end)    ];
       obj.physicslist.helicity  = [obj.physicslist.helicity(1:ivol-1),  obj.physicslist.helicity(ivol+1:end) ];
       obj.physicslist.pressure  = [obj.physicslist.pressure(1:ivol-1),  obj.physicslist.pressure(ivol+1:end) ];
       obj.physicslist.adiabatic = [obj.physicslist.adiabatic(1:ivol-1), obj.physicslist.adiabatic(ivol+1:end)];
       obj.physicslist.mu        = [obj.physicslist.mu(1:ivol-1),        obj.physicslist.mu(ivol+1:end)       ];
       obj.physicslist.Ivolume   = [obj.physicslist.Ivolume(1:ivol-1),   obj.physicslist.Ivolume(ivol+1:end)  ];
       obj.physicslist.Isurf     = [obj.physicslist.Isurf(1:ivol-1),     obj.physicslist.Isurf(ivol+1:end)    ];
       obj.physicslist.iota      = [obj.physicslist.iota(1:ivol),        obj.physicslist.iota(ivol+2:end)     ];
       obj.physicslist.oita      = [obj.physicslist.oita(1:ivol),        obj.physicslist.oita(ivol+2:end)     ];
       obj.diagnosticslist.nPtrj = [obj.diagnosticslist.nPtrj(1:ivol-1), obj.diagnosticslist.nPtrj(ivol+1:end)];
       
       
       obj.initial_guess.Ric = [obj.initial_guess.Ric(:,1:ivol-1), obj.initial_guess.Ric(:,ivol+1:end)];
       obj.initial_guess.Ris = [obj.initial_guess.Ris(:,1:ivol-1), obj.initial_guess.Ris(:,ivol+1:end)];
       obj.initial_guess.Zic = [obj.initial_guess.Zic(:,1:ivol-1), obj.initial_guess.Zic(:,ivol+1:end)];
       obj.initial_guess.Zis = [obj.initial_guess.Zis(:,1:ivol-1), obj.initial_guess.Zis(:,ivol+1:end)];
       
       
    end
    
    
    % =====================================================================
    % Write method
    function write_input_file(obj, filename )
      %
	  % write_input_file(filename, (force) )
	  %
	  % Write namelist in an input file. Be careful: if the file filename
	  % provided as input exist, it will be overwritten!
	  %
	  % INPUT
	  % -----
	  %   - filename: path where to save the input file.
	  %   - options: structure containing optional input
      %       * force: set to true to force overwritting the input file
      %


        nlists = length(obj.lists);
        S = struct;
        for ii=1:nlists
            S.(obj.lists{ii}) = obj.(obj.lists{ii});
        end
        
        S.physicslist = rmfield( S.physicslist, 'im' );
        S.physicslist = rmfield( S.physicslist, 'in' );
        try
            S.physicslist = rmfield( S.physicslist, 'PlasmaBoundary' );
        catch
            disp('No field PlasmaBoundary...')
        end
        
        try
            S.physicslist = rmfield( S.physicslist, 'ComputationalBoundary' );
        catch
            disp('No field ComputationalBoundary...')
        end
        
        if isfield(S.physicslist, 'Lboundary')
            Lboundary = S.physicslist.Lboundary;
        else
            Lboundary = 0; %by default
        end
        
        Ntor = max(max(obj.physicslist.in));
        S.shift.Rbc = [1, Ntor+1];
        S.shift.Rbs = [1, Ntor+1];
        S.shift.Zbc = [1, Ntor+1];
        S.shift.Zbs = [1, Ntor+1];
        S.shift.Rwc = [1, Ntor+1];
        S.shift.Rws = [1, Ntor+1];
        S.shift.Zws = [1, Ntor+1];
        S.shift.Zwc = [1, Ntor+1];
        S.shift.Vnc = [1, Ntor+1];
        S.shift.Vns = [1, Ntor+1];
        S.shift.Bnc = [1, Ntor+1];
        S.shift.Bns = [1, Ntor+1];
        S.shift.rhomn = [1, Ntor+1];
        
        if Lboundary==0        

            initialguess = cell(1,1);
            if ~isempty(obj.initial_guess)
                s = size(obj.initial_guess.Ric);
                initialguess = cell(1, s(1)*s(2));
                n=0;
                for ii=1:s(1)
                    for jj=1:s(2)
                        n = n+1;
                        initialguess{n} = sprintf( '%i   %i   ', obj.initial_guess.im(ii,jj), obj.initial_guess.in(ii,jj) );
                        for ivol=1:s(3)
                           initialguess{n} = sprintf( '%s   %0.12E   %0.12E   %0.12E   %0.12E', initialguess{n}, ...
                                                        obj.initial_guess.Ric(ii,jj,ivol), ...
                                                        obj.initial_guess.Zis(ii,jj,ivol), ...
                                                        obj.initial_guess.Ris(ii,jj,ivol), ...
                                                        obj.initial_guess.Zic(ii,jj,ivol)     );
                        end
                    end
                end
            end
            
        else
            
            nrho = max(max(obj.physicslist.in)); 
            S.shift.rhomn = [ 1, nrho+1 ];
            
            % Now prepare initial guess...
            if isfield(obj.initial_guess, 'bin')
                sb   = size(obj.initial_guess.bin);
                srho = size(obj.initial_guess.rhoi);
                ntor = (srho(2)-1)/2.0;
                initialguess = cell(1, sb(1) + srho(1)*srho(2));

                % Modes m=0, n
                for ii=1:sb(1)
                   nn = ii-1;
                   initialguess{ii} = sprintf('0   %i   ', nn);
                   for ivol=1:sb(2)
                      initialguess{ii} = sprintf('%s   %0.12E   %0.12E   %0.12E   %0.12E', initialguess{ii}, ...
                                                 obj.initial_guess.bin( ii, ivol ), ...
                                                 obj.initial_guess.R0ic(ii, ivol ), ...
                                                 obj.initial_guess.Z0is(ii, ivol ), 0.0 ); 
                   end
                end

                it = sb(1);
                for ii=1:srho(2)

                    nn = ii-ntor-1;

                    for jj=1:srho(1)

                        it = it+1;

                        mm = jj;

                        initialguess{it} = sprintf('%i   %i   ', mm, nn);

                        for ivol=1:srho(3)

                            initialguess{it} = sprintf('%s   %0.12E   %0.12E   %0.12E   %0.12E', initialguess{ii},...
                                                       0.0, 0.0, 0.0, obj.initial_guess.rhoi(ii, jj, ivol));

                        end
                    end
                end
            else
                initialguess = cell(0);
                
            end
            
            
            
        end
        
        write_namelist( S, filename, initialguess );
        

%       % Read options
%       if( isempty(varargin) )
%         force = false;
%       else
%         if( isfield(varargin{1}, 'force') )
%           force = varargin{1}.force;
%         else
%           force = false;
%         end
%       end
% 
% 			% Check if user want to erase existing file
% 			if (exist(filename, 'file') && (~force))
%                 prompt = ['The file ', filename, ' already exist and will be overwritten. Continue? Y/N [Y]:  '];
% 				c = input(prompt, 's');
%                 if( isempty(c) )
%                   c = 'Y';
%                 end
%                 
% 				while ~( (c=='Y') || (c=='N') )
% 					c = input('Please answer Y or N:  ', 's');
% 				end
% 
% 				if strcmp(c,'N') %exit routine
% 					disp('ABORTED')
% 					return;
%         else
%           disp(['Overwriting file ', filename, '.'])
% 				end
% 			end
% 
% 			% Open file
% 			fid = fopen(filename, 'w');
% 
% 
% 			special_list = {'Rbc', 'Rbs', 'Zbc', 'Zbs', 'im', 'in', ...
%         				    'Rwc', 'Rws', 'Zwc', 'Zws', ...
%         					'Vnc', 'Vns', 'Bnc', 'Bns', ...
%                             'PlasmaBoundary', 'ComputationalBoundary'};
% 			for ii=1:numel(obj.lists)
% 				list = obj.lists{ii};
% 
% 				if strcmp(list, 'physicslist')
% 					obj.write_list(fid, list, special_list)
% 
% 					% write Rbc, ...
% 					obj.write_4_values(fid, 'Rbc', 'Zbs', 'Rbs', 'Zbc');
% 					obj.write_4_values(fid, 'Rwc', 'Zws', 'Rws', 'Zwc');
% 					obj.write_4_values(fid, 'Vnc', 'Vns', 'Bnc', 'Bns');
% 
% 					fprintf(fid, '/ \n');
% 
% 				else
% 					obj.write_list(fid, list, special_list)
% 
% 					fprintf(fid, '/ \n');
%                 end
% 
%             end
% 
% 
%             im = obj.physicslist.im; in = obj.physicslist.in;
%             mn = length(im);
% 
%             if( ~isempty(obj.initial_guess) )
%               Ric = obj.initial_guess.Ric;
%               s = size(Ric);
% 
%               for ii=1:mn
%                   fprintf(fid, '   %i   %i', im(ii), in(ii));
%                   for jj=1:s(2)
%                       fprintf(fid, '   %19.15E   %19.15E   %19.15E   %19.15E', ...
%                               obj.initial_guess.Ric(ii,jj), ...
%                               obj.initial_guess.Zis(ii,jj), ...
%                               obj.initial_guess.Ris(ii,jj), ...
%                               obj.initial_guess.Zic(ii,jj)); % Print in file
%                   end
%                   fprintf(fid, '\n');
%               end
%             end
% 
% 			% close file
% 			fclose(fid);

        end %of function write_input_file
        
    
    function ii = find_ii_given_mn(obj, m, n)

      im = obj.physicslist.im;
      in = obj.physicslist.in;
      ii = find(double(im==m) .* double(in==n));
        end

  end % of public methods
	% ----------------------------------------------------------------------
	%                           PRIVATE METHODS
	% ======================================================================
	%
	methods (Access = private)
      
        function obj = update_flux_surfaces( obj )

            im = obj.physicslist.im;
            in = obj.physicslist.in;
            Mpol = max(max(im)); Ntor = max(max(in));

            PlasmaBoundary = fluxSurface(obj.physicslist.Nfp, Mpol, Ntor, 0);
            PlasmaBoundary = PlasmaBoundary.set_array('rmnc', obj.physicslist.Rbc, im, in);
            PlasmaBoundary = PlasmaBoundary.set_array('rmns', obj.physicslist.Rbs, im, in);
            PlasmaBoundary = PlasmaBoundary.set_array('zmnc', obj.physicslist.Zbc, im, in);
            PlasmaBoundary = PlasmaBoundary.set_array('zmns', obj.physicslist.Zbs, im, in);
            obj.physicslist.PlasmaBoundary = PlasmaBoundary;

            CompBoundary = fluxSurface(obj.physicslist.Nfp, Mpol, Ntor, 0);
            CompBoundary = CompBoundary.set_array('rmnc', obj.physicslist.Rwc, im, in);
            CompBoundary = CompBoundary.set_array('rmns', obj.physicslist.Rws, im, in);
            CompBoundary = CompBoundary.set_array('zmnc', obj.physicslist.Zwc, im, in);
            CompBoundary = CompBoundary.set_array('zmns', obj.physicslist.Zws, im, in);
            obj.physicslist.ComputationalBoundary = CompBoundary;
        end

		function [value, field] = extract_4_values(obj, tline, n, m, str1, str2, str3, str4)
			% Extract rbc, rbs, zbc, zbs components from input file. Only
			% used in the class constructor
			work = tline;


			% Remove (n,m) from line
            p_open  = find(work=='(');
            p_close = find(work==')');

            nfield = length(p_open);
            if(nfield>4)
               error(['Bad format on line: ', work,', need only four fields']);
            end

            work = tline(1:p_open(1)-1);
            for ii=2:nfield
                work = [work, tline(p_close(ii-1)+1:p_open(ii)-1)];
            end
            work = [work, tline(p_close(end)+1:end)];


			nm_str = ['(',num2str(n),',',num2str(m),')'];
			work = erase(work, nm_str);

			% Remove = signs as well
			work = erase(work, '=');

			% Now read line
            format = '';
            for ii=1:nfield
               format = [format, '%s%f64']; 
            end
			work = textscan(work, format);

            for ii=1:nfield
				field{ii} = work{2*ii-1}{1};
				value{ii} = work{2*ii};
            end
            
            
            
            
		end % of function extract_4_values


    function [n, m] = get_nm_from_field(obj, field);
      % Read (n,m) from field Rbc(n,m), Rwc(n,m) or Vnc(n,m)

      s = extractBefore(field, ')');
      s = extractAfter( s,     '(');

      n_str = extractBefore(s, ',');
      m_str = extractAfter(s, ',');

      n = str2num(char(n_str));
      m = str2num(char(m_str));
    end


		function write_4_values(obj, fid, field1, field2, field3, field4)

      im = obj.physicslist.im;
      in = obj.physicslist.in;

      for ii=1:length(im)
          m = im(ii);
          n = in(ii);

          form = ['  ', field1, '(%i,%i) =  %19.15E  ', ...
                        field2, '(%i,%i) =  %19.15E  ', ...
                        field3, '(%i,%i) =  %19.15E  ', ...
                        field4, '(%i,%i) =  %19.15E \n'];


          value1 = obj.physicslist.(field1)(ii);
          value2 = obj.physicslist.(field2)(ii);
          value3 = obj.physicslist.(field3)(ii);
          value4 = obj.physicslist.(field4)(ii);

          fprintf(fid, form, n, m, value1, ...
	                 n, m, value2, ...
                             n, m, value3, ...
                             n, m, value4     );
      end % of for loop
		end % of function write_4_values


		function write_list(obj, fid, list, special_list)
      if( isempty(obj.(list)) )
          return
      end
			fields = fieldnames(obj.(list));

			% Beginning of list
			% End of list - add a \
			fprintf(fid, '&%s \n', list);
			for ii = 1:numel(fields)
				f = fields(ii);
				f = f{1};
                % deal with special_list outside this method
				if any(strcmp(f,special_list))
					continue;
				end

				% write line
				value = obj.(list).(f);
				l     = length(value);

				fprintf(fid, '  %s  =    ', f);
				for jj=1:l
					c = class(value(jj));

					% Few logical test to decide the printing format
          if strcmp(c,'logical') % bool
            if value(jj)
              fprintf(fid, 'T    ');
            else
              fprintf(fid, 'F    ');
            end
          else % not a bool
            if floor(value(jj))==value(jj) %Integer value
              fprintf(fid, '%i    ',value(jj));
            else % double
              fprintf(fid, '%19.15E    ',value(jj));
            end
          end
				end
				fprintf(fid,'\n');
			end % of for loop

		end % of function write_list
  end % of private methods
end % of SPEC_Namelist class
