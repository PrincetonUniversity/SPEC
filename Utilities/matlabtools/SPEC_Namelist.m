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
    function obj = SPEC_Namelist(inputfile, varargin)

        if( ~isempty(varargin) )
          input_category = varargin{1};
        end

        fid = fopen(inputfile,'r'); % open template file
        tline = fgetl(fid);

        obj.lists = {'physicslist', 'numericlist', 'locallist', ...
                     'globallist', 'diagnosticslist', 'screenlist'};

        for ii=numel(obj.lists)
            list = obj.lists{ii};
            obj.(list) = struct;
        end
        category = '';

        boundary.Rbc = []; boundary.Rbs = []; boundary.Zbc = []; boundary.Zbs = []; ib=0; bm = []; bn = [];
        wall.Rwc     = []; wall.Rws     = []; wall.Zwc     = []; wall.Zws     = []; iw=0; wm = []; wn = [];
        Nfield.Vns   = []; Nfield.Vnc   = []; Nfield.Bnc   = []; Nfield.Bns   = []; iv=0; vm = []; vn = [];

        iline = 0;
        while ~isempty(tline) || feof(fid)

          % Remove trailing and leading spaces
          tline = strtrim(tline);

          if feof(fid)
             break;
          end

          iline = iline + 1		;
          % jump / lines between categories.
          if strcmp(tline, '/')
             if strcmp(category,'screenlist')
                break; % end of file, ignore geometry initial guess
             else
                 tline = fgetl(fid);
                 continue;
             end
          end

          % read category
          if strcmp(tline(end-3:end),'list')
             category = tline(2:end);
             tline = fgetl(fid);
             continue;
          end


          % Field
          try
              field = extractBefore(tline,'=');  %read field
              field = char(field);
              field = field(~isspace(field));    %remove spaces
          catch
              disp(['Error in field reading at line ', num2str(iline)])
              disp(['Field is ', field])
          end


          % Content
          if length(field)<3

            try
              content_str = extractAfter(tline,'=');
              content = char(content_str);
              content = content(~isspace(content));
              if strcmp(content,'F')
                  content = false;
              elseif strcmp(content,'T')
                  content = true;
              else
                  content = str2num(char(content_str));
              end

              % Store in structure
              obj.(category).(field) = content;

            catch
                error(['Error reading line ', tline])

            end


          else
            if strcmp(field(1:3),'Rbc') %need to split Rbc, Rbs, Zbc, Zbs
              ib = ib + 1;
              [n, m] = obj.get_nm_from_field(field);

              [value, field]  = obj.extract_4_values(tline, n, m, 'Rbc', 'Zbs', 'Rbs', 'Zbc');

              boundary.Rbc(ib) = 0.0;
              boundary.Zbc(ib) = 0.0;
              boundary.Zbs(ib) = 0.0;
              boundary.Rbs(ib) = 0.0;
              
              nfield = length(field);
              for ii=1:nfield
                  boundary.(field{ii})(ib) = value{ii};
              end

              bm(ib) = m;
              bn(ib) = n;

            elseif strcmp(field(1:3),'Rwc') % same as above
              iw = iw + 1;
              [n, m] = obj.get_nm_from_field(field);

              [value, field]  = obj.extract_4_values(tline, n, m, 'Rwc', 'Zws', 'Rws', 'Zwc');

              wall.Rwc(iw) = 0.0;
              wall.Zwc(iw) = 0.0;
              wall.Zws(iw) = 0.0;
              wall.Rws(iw) = 0.0;
              
              nfield = length(field);
              for ii=1:nfield
                  wall.(field{ii})(iw) = value{ii};
              end

              wm(iw) = m;
              wn(iw) = n;

            elseif strcmp(field(1:3),'Vns') || strcmp(field(1:3),'Vnc') || strcmp(field(1:3),'Bns') || strcmp(field(1:3),'Bnc')% same as above
              iv = iv + 1;
              [n, m] = obj.get_nm_from_field(field);

              [value, field]  = obj.extract_4_values(tline, n, m, 'Vns', 'Bns', 'Vnc', 'Bnc');

              Nfield.Vnc(iv) = 0.0;
              Nfield.Bnc(iv) = 0.0;
              Nfield.Vns(iv) = 0.0;
              Nfield.Bns(iv) = 0.0;
              
              nfield = length(field);
              for ii=1:nfield
                  Nfield.(field{ii})(iv) = value{ii};
              end

              vm(iv) = m;
              vn(iv) = n;

            else % easier, just read content
              content_str = extractAfter(tline,'=');

              content = char(content_str);
              content = content(~isspace(content));
              if strcmp(content,'F')
                  content = false;
              elseif strcmp(content,'T')
                  content = true;
              else
                  content = str2num(char(content_str));
              end

              % Store in structure
              if isempty(category)
                obj.(input_category).(field) = content;
              else
                obj.(category).(field) = content;
              end
            end
          end

          % get next line
          tline = fgetl(fid);
        end

        % Store rbc, rbs ... in structure

        boundaryfields = fieldnames(boundary);
        wallfields     = fieldnames(wall);
        Nfieldfields   = fieldnames(Nfield);

        for ii=1:4
            obj.physicslist.(boundaryfields{ii}) = boundary.(boundaryfields{ii});
            obj.physicslist.(wallfields{ii})     = wall.(wallfields{ii});
            obj.physicslist.(Nfieldfields{ii})   = Nfield.(Nfieldfields{ii});
        end

        fclose(fid);                % close template file

        % Reformat the Rbc, Zbs, ... arrays
        if(isfield(obj.physicslist, 'Mpol'))
          Mpol = obj.physicslist.Mpol;
          Ntor = obj.physicslist.Ntor;
        else
          Mpol = max([wm,bm,vm]); Ntor = max([wn,bn,vn]);
        end

        mn   = 1 + Ntor +  Mpol * ( 2 *  Ntor + 1 );

        Rbc = zeros(1,mn); Rbs = zeros(1,mn); Zbc = zeros(1,mn); Zbs = zeros(1,mn);
        Rwc = zeros(1,mn); Rws = zeros(1,mn); Zwc = zeros(1,mn); Zws = zeros(1,mn);
        Vnc = zeros(1,mn); Vns = zeros(1,mn); Bnc = zeros(1,mn); Bns = zeros(1,mn);
        im = zeros(1,mn);
        in = zeros(1,mn);

        im(1:1+Ntor) = 0;
        in(1:1+Ntor) = 0:Ntor;
        count = Ntor+1;
        for ii=1:Mpol
          for jj=-Ntor:Ntor
            count = count + 1;
            im( count ) = ii;
            in( count ) = jj;
          end
        end

        for imn = 1:mn
          m = im(imn);
          n = in(imn);

          b_ind = find(double(bm==m) .* double(bn==n));
          if( ~isempty(b_ind) )
            Rbc( imn ) = obj.physicslist.Rbc( b_ind );
            Rbs( imn ) = obj.physicslist.Rbs( b_ind );
            Zbc( imn ) = obj.physicslist.Zbc( b_ind );
            Zbs( imn ) = obj.physicslist.Zbs( b_ind );
          end

          w_ind = find(double(wm==m) .* double(wn==n));
          if( ~isempty(w_ind) )
            Rwc( imn ) = obj.physicslist.Rwc( w_ind );
            Rws( imn ) = obj.physicslist.Rws( w_ind );
            Zwc( imn ) = obj.physicslist.Zwc( w_ind );
            Zws( imn ) = obj.physicslist.Zws( w_ind );
          end

          v_ind = find(double(vm==m) .* double(vn==n));
          if( ~isempty(v_ind) )
            Vnc( imn ) = obj.physicslist.Vnc( v_ind );
            Vns( imn ) = obj.physicslist.Vns( v_ind );
            Bnc( imn ) = obj.physicslist.Bnc( v_ind );
            Bns( imn ) = obj.physicslist.Bns( v_ind );
          end
        end

        obj.physicslist.im = im;
        obj.physicslist.in = in;

        obj.physicslist.Rbc = Rbc;
        obj.physicslist.Rbs = Rbs;
        obj.physicslist.Zbc = Zbc;
        obj.physicslist.Zbs = Zbs;
        
        obj.physicslist.Rwc = Rwc;
        obj.physicslist.Rws = Rws;
        obj.physicslist.Zwc = Zwc;
        obj.physicslist.Zws = Zws;

        obj.physicslist.Vnc = Vnc;
        obj.physicslist.Vns = Vns;
        obj.physicslist.Bnc = Bnc;
        obj.physicslist.Bns = Bns;
        
        obj = obj.update_flux_surfaces();

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

        field = extractBefore(tline,'=');  %read field
        field = char(field);
        field = field(~isspace(field));    %remove spaces
        content_str = extractAfter(tline,'=');
        content = str2num(content_str);

        % Read initial guess for coordinate axis
        switch field
        case 'Rac'
          obj.physicslist.Rac = content;
        case 'Ras'
          obj.physicslist.Ras = content;
        case 'Zac'
          obj.physicslist.Zac = content;
        case 'Zas'
          obj.physicslist.Zas = content;
        end


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

      % Allocate memory
      im = obj.physicslist.im; in = obj.physicslist.in;
      mn = length( im );
      Nvol = obj.physicslist.Nvol;
      Mvol = Nvol + obj.physicslist.Lfreebound;

      Ric = zeros(mn,Mvol); Zis = zeros(mn,Mvol);
      Ris = zeros(mn,Mvol); Zic = zeros(mn,Mvol);

      format = '%f%f';
      for ii=1:Nvol
         format = [format, '%f%f%f%f'];
      end


      % Read
      for iline=1:nlines

          % Scan line
          line_data = textscan( initial_guess{iline}, format );

          % Find corresponding index
          m = line_data{1}; n = line_data{2};
          imn = obj.find_ii_given_mn(m, n);

          for ivol=1:Nvol
              Ric(imn,ivol) = line_data{ivol*4-1};
              Zis(imn,ivol) = line_data{ivol*4  };
              Ris(imn,ivol) = line_data{ivol*4+1};
              Zic(imn,ivol) = line_data{ivol*4+2};
          end
      end

      init.Ric = Ric;
      init.Ris = Ris;
      init.Zis = Zis;
      init.Zic = Zic;

      obj.initial_guess = init;

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
    
    
    % =====================================================================
    % Plotters
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
              'Invalid input'
      end
      
      Nfp = obj.physicslist.Nfp;
      theta = linspace(0, 2*pi, Nt );
      R = zeros(1,Nt); 
      Z = zeros(1,Nt);
      for ii=1:length(im)
         
          if im(ii) > Mpol
              continue
          end
          
          if in(ii) > Ntor
              continue
          end
          
         arg = im(ii) * theta - in(ii) * Nfp * phi;
         cosarg = cos(arg);
         sinarg = sin(arg);
         
         R = R + Rmnc(ii) * cosarg + Rmns(ii) * sinarg;
         Z = Z + Zmnc(ii) * cosarg + Zmns(ii) * sinarg;
                               
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


%        if isempty(obj.initial_guess)
%           error('No initial guess. First read with read_initial_guess');
%        end


       % Open figure
       switch newfig
       case 0
         hold on
       case 1
         figure
         hold on;
       case 2
         hold off;
       otherwise
         error( 'Unknown value of newfig')
       end

       Nvol   = obj.physicslist.Nvol;
       im = obj.physicslist.im; in = obj.physicslist.in;
       mn = length(im);
       Np = obj.physicslist.Nfp;

      % Prepare theta coordinate
      theta = linspace(0,2*pi,Nt)';

       % Map to real space
       if( ~isempty(obj.initial_guess) )
           for ivol=1:Nvol
               R = zeros(1,Nt);
               Z = zeros(1,Nt);

               for imn=1:mn
                   cosarg = cos(im(imn)*theta - in(imn)*Np*phi);
                   sinarg = sin(im(imn)*theta - in(imn)*Np*phi);

                 
                   R(:) = R(:) + obj.initial_guess.Ric(imn,ivol) * cosarg ...
                             + obj.initial_guess.Ris(imn,ivol) * sinarg;
                   Z(:) = Z(:) + obj.initial_guess.Zic(imn,ivol) * cosarg ...
                             + obj.initial_guess.Zis(imn,ivol) * sinarg;
               end

               scatter(R, Z, 3, 'filled', 'MarkerFaceColor', 'b')
           end
       end
       
       % Plot plasma boundary
       R = zeros(1,Nt);
       Z = zeros(1,Nt);

       for imn=1:mn
           cosarg = cos(im(imn)*theta - in(imn)*Np*phi);
           sinarg = sin(im(imn)*theta - in(imn)*Np*phi);

           try
             R(:) = R(:) + obj.physicslist.Rbc(imn) * cosarg ...
                         + obj.physicslist.Rbs(imn) * sinarg;
             Z(:) = Z(:) + obj.physicslist.Zbc(imn) * cosarg ...
                         + obj.physicslist.Zbs(imn) * sinarg;
           catch
             disp('ouch')
           end
       end

       scatter(R, Z, 3, 'filled', 'MarkerFaceColor', [1 0 0])


       if obj.physicslist.Lfreebound
           R = zeros(1,Nt);
           Z = zeros(1,Nt);

           Rwc   = obj.physicslist.Rwc;
           Zwc   = obj.physicslist.Zwc;
           Rws   = obj.physicslist.Rws;
           Zws   = obj.physicslist.Zws;

           for imn=1:mn
               cosarg = cos(im(imn)*theta - in(imn)*Np*phi);
               sinarg = sin(im(imn)*theta - in(imn)*Np*phi);

               R(:) = R(:) + Rwc(imn) * cosarg + Rws(imn) * sinarg;
               Z(:) = Z(:) + Zwc(imn) * cosarg + Zws(imn) * sinarg;
           end

           scatter(R, Z, 5, 'filled', 'MarkerFaceColor', [0 0 0])
       end

       %Plot the magnetic axis guess
       Rax = 0; Zax = 0;
       mn = length(obj.physicslist.Rac);
       Nfp = double(obj.physicslist.Nfp);
       for n=0:mn-1
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

		  ii = obj.find_ii_given_mn(m, n);

		  if isnan(ii)
			  out.Rbc = 0.0;
			  out.Rbs = 0.0;
			  out.Zbc = 0.0;
			  out.Zbs = 0.0;
		  else
			  out.Rbc = obj.physicslist.Rbc(ii);
			  out.Rbs = obj.physicslist.Rbs(ii);
			  out.Zbc = obj.physicslist.Zbc(ii);
			  out.Zbs = obj.physicslist.Zbs(ii);
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
		  ii = obj.find_ii_given_mn(m, n);

		  if isnan(ii)
			  out.Rwc = 0.0;
			  out.Rws = 0.0;
			  out.Zwc = 0.0;
			  out.Zws = 0.0;
		  else
			  out.Rwc = obj.physicslist.Rwc(ii);
			  out.Rws = obj.physicslist.Rws(ii);
			  out.Zwc = obj.physicslist.Zwc(ii);
			  out.Zws = obj.physicslist.Zws(ii);
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

		  ii = obj.find_ii_given_mn(m, n);

		  if isnan(ii)
			  out.Vnc = 0.0;
			  out.Vns = 0.0;
			  out.Bnc = 0.0;
			  out.Bns = 0.0;
		  else
			  out.Vnc = obj.physicslist.Vnc(ii);
			  out.Vns = obj.physicslist.Vns(ii);
			  out.Bnc = obj.physicslist.Bnc(ii);
			  out.Bns = obj.physicslist.Bns(ii);
		  end
    end
    

    function out = get_fourier_harmonics(obj, m, n, field)

      ii = obj.find_ii_given_mn( m, n);

      if( strcmp(field, 'Ric') || strcmp(field, 'Ris') || strcmp(field, 'Zic') || strcmp(field, 'Zis'))
        if( isempty(obj.initial_guess) )
          error('No initial guess')
        else
          out = obj.initial_guess.(field)(ii);
        end
      else
        out = obj.physicslist.(field)(ii);
      end
    end
    

	% =====================================================================
    % Setters
    function obj = set_fourier_harmonics(obj, im, in, values, field)

      mn = length(im);
      mn_old = length(obj.physicslist.Rbc);

      obj.physicslist.(field) = zeros(1, mn_old);

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
        
        obj = obj.update_flux_surfaces();

      end


    end
    
    
    function obj = set_initial_guess_harmonics( obj, volume, im, in, values, field)
        
      mn = length(im);
      mn_old = length(obj.physicslist.Rbc);

      Mvol = obj.physicslist.Nvol + obj.physicslist.Lfreebound;
      obj.initial_guess.(field) = zeros(mn_old, Mvol);

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

          obj.initial_guess.(field)(end, volume) = values(ii);

        else
          obj.initial_guess.(field)(imn, volume) = values(ii);

        end
        

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
    function write_input_file(obj, filename, varargin )
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

      % Read options
      if( isempty(varargin) )
        force = false;
      else
        if( isfield(varargin{1}, 'force') )
          force = varargin{1}.force;
        else
          force = false;
        end
      end

			% Check if user want to erase existing file
			if (exist(filename, 'file') && (~force))
        prompt = ['The file ', filename, ' already exist and will be overwritten. Continue? Y/N [Y]:  '];
				c = input(prompt, 's');
        if( isempty(c) )
          c = 'Y';
        end
				while ~( strcmp(c,"Y") || strcmp(c,"N") )
					c = input(['Please answer Y or N:  '], 's');
				end

				if strcmp(c,'N') %exit routine
					disp('ABORTED')
					return;
        else
          disp(['Overwriting file ', filename, '.'])
				end
			end

			% Open file
			fid = fopen(filename, 'w');


			special_list = {'Rbc', 'Rbs', 'Zbc', 'Zbs', 'im', 'in', ...
        				    'Rwc', 'Rws', 'Zwc', 'Zws', ...
        					'Vnc', 'Vns', 'Bnc', 'Bns', ...
                            'PlasmaBoundary', 'ComputationalBoundary'};
			for ii=1:numel(obj.lists)
				list = obj.lists{ii};

				if strcmp(list, 'physicslist')
					obj.write_list(fid, list, special_list)

					% write Rbc, ...
					obj.write_4_values(fid, 'Rbc', 'Zbs', 'Rbs', 'Zbc');
					obj.write_4_values(fid, 'Rwc', 'Zws', 'Rws', 'Zwc');
					obj.write_4_values(fid, 'Vnc', 'Vns', 'Bnc', 'Bns');

					fprintf(fid, '/ \n');

				else
					obj.write_list(fid, list, special_list)

					fprintf(fid, '/ \n');
                end

            end


            im = obj.physicslist.im; in = obj.physicslist.in;
            mn = length(im);

            if( ~isempty(obj.initial_guess) )
              Ric = obj.initial_guess.Ric;
              s = size(Ric);

              for ii=1:mn
                  fprintf(fid, '   %i   %i', im(ii), in(ii));
                  for jj=1:s(2)
                      fprintf(fid, '   %19.15E   %19.15E   %19.15E   %19.15E', ...
                              obj.initial_guess.Ric(ii,jj), ...
                              obj.initial_guess.Zis(ii,jj), ...
                              obj.initial_guess.Ris(ii,jj), ...
                              obj.initial_guess.Zic(ii,jj)); % Print in file
                  end
                  fprintf(fid, '\n');
              end
            end

			% close file
			fclose(fid);

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
            Mpol = max(im); Ntor = max(in);

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
