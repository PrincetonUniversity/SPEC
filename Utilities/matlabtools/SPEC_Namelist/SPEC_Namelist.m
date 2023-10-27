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
    
    properties (Access=private)
        mpol = 0;
        ntor = 0;
        array_size = [0, 0];
        Mvol = 0;
        nvol = 0;
        verbose = true;
        lboundary = 0;
        
    end
    
    methods (Access=public)
        % Class constructor
        function obj = SPEC_Namelist( filename, varargin )
            %
            % SPEC_NAMELIST( FILENAME )
            % =========================
            %
            % Class constructor of the class SPEC_Namelist. This class
            % reads the SPEC input file, checks that the data is correctly
            % formatted, allow some plottings and easy changes in the input
            % file, and provides a routine to write SPEC input file.
            %
            % INPUTS
            % ------
            %   -FILENAME: SPEC input file (.sp)
            %   -VArArGIN: Any couple of input:
            %       - 'Liniguess': set to true to read initial guess, false
            %                      to skip it. default: true
            %                       
            %
            % OUTPUTS
            % -------
            %   -OBJ: An instance of the class SPEC_Namelist
            %
            %
        
            % read optional input
            l = length(varargin);
            if mod(l,2)~=0
                error('InputError: invalid number of inputs')
            end
            
            opt.Liniguess = true; % Decide whether or not we read initial guess
            opt.verbose = true;   % Print additional warnings
            for ii=1:l/2
               opt.(varargin{2*ii-1}) = varargin{2*ii}; 
            end
            
            obj.verbose = opt.verbose;
            
            % read input file
            work = read_namelist( filename );
            obj.lists = fields(work);
            for ii=1:length(obj.lists)
                obj.(obj.lists{ii}) = work.(obj.lists{ii});
            end
            
            % Check that the size of arrays makes sense, fills with zeros
            % otherwise
            obj = obj.initialize_structure();
            
            % Find the largest Fourier resolution used in the input file;
            % reformat all spectral quantities to have the same resolution
            obj = obj.set_fourier_resolution();
            
            % read initial guess
            if opt.Liniguess
                obj = obj.read_initial_guess( filename );
            else
                obj.initial_guess = struct([]);
            end
            
        end
        
        function obj = read_initial_guess( obj, filename )
            %
            % READ_INITIAL_GUESS( FILENAME )
            % ==============================
            %
            % read the initial guess from the input file filename. 
            %
            % If the Fourier resolution of the initial guess is larger than
            % the Fourier resolution of the other physical qunatities, these
            % are extended with zeros. 
            %
            % If the initial guess requires less Fourier harmonics, then it 
            % is extended with zeros to match the Fourier resolution of the
            % other physical quantities
            %
            % If no initial guess are provided in filename, return an empty
            % structure.
            
            
            % First, open file and read relevant portion. We look for the
            % end of the category "screenlist" and then read each line.
            % Each line is then saved as a string in a structure.
            fid = fopen( filename, 'r' );
            save_line = false; % This is switched to true once 
                               % we are at the end of the screelist
            category  = '';    % Save the name of the category

            initial_guess_str = {};

            while( ~feof(fid) ) % while it is not the end of the file

                tline = fgetl(fid);
                tline = strtrim(tline);
                
                if isempty( tline )
                    continue
                end

                if( save_line ) % write line
                    initial_guess_str{end+1} = tline;

                else % otherwise look for beginning of initial guess
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

            fclose(fid); % Close file
            
            % Check if structure is empty - i.e. no initial guess is
            % provided
            if isempty( initial_guess_str )
                obj.initial_guess = struct([]); % generate empty structure
                return                
            end
                        
            % Check size. The number of elements should be 2 + 4 * Nvol.
            % The first two elements are the mode numbers m and n, then we
            % have Rbc, Zbs, Rbs, Zbc for each volume interface
            l = length(str2num( initial_guess_str{1} ));
            l = l-2; % remove m, n
            if mod(l,4)~=0
                error('Invalid number of modes in initial guess')
            end
            
            if l/4<obj.nvol
                error('Not enough volumes in initial guess')
            end
            
            if l/4>obj.nvol
                error('Too many volumes in initial guess')
            end
            
            % read format of initial guess
            tmp = SPEC_Namelist( filename, 'Liniguess', false );
            
            if tmp.physicslist.lboundary~=obj.physicslist.lboundary
                error(['The initial guess from file %s does not use the',...
                       'same boundary representation.'], filename)
            end
            
            % Check if there is an initial guess
            nlines = length(initial_guess_str);
            if nlines<1
                obj.initial_guess = struct([]); % generate empty structure
                return
            end

            % read mpol, ntor
            mpol_in = 0; ntor_in = 0;
            for iline=1:nlines
                % Scan line
                line_data = str2num( initial_guess_str{iline} );

                % Find corresponding index
                mm = line_data(1); nn = line_data(2);
                mpol_in = max([mm, mpol_in]);
                ntor_in = max([abs(nn), ntor_in]);
            end
            
            % Check if resolution is smaller or larger than inner
            % resolution. This changes obj.mpol and obj.ntor if
            % necessary
            if (mpol_in>obj.mpol) || (ntor_in>obj.ntor)
               obj = obj.change_fourier_resolution( mpol_in, ntor_in );
            end
            
            % Now format initial guess in a structure
            switch obj.physicslist.lboundary
                case 0 % rmn, zmn representation
                    
                    % Allocate memory
                    ric = zeros(2*obj.ntor+1, obj.mpol+1, obj.Mvol);
                    ris = zeros(2*obj.ntor+1, obj.mpol+1, obj.Mvol);
                    zic = zeros(2*obj.ntor+1, obj.mpol+1, obj.Mvol);
                    zis = zeros(2*obj.ntor+1, obj.mpol+1, obj.Mvol);
                    
                    % Fill initial guess arrays
                    for iline=1:nlines

                        % Scan line
                        line_data = str2num( initial_guess_str{iline} );

                        % Find corresponding index
                        m = line_data(1); n = line_data(2);
                        im = m+1;
                        in = n+obj.ntor+1;

                        for ivol=1:obj.nvol
                            ric(in,im,ivol) = line_data(ivol*4-1);
                            zis(in,im,ivol) = line_data(ivol*4  );
                            ris(in,im,ivol) = line_data(ivol*4+1);
                            zic(in,im,ivol) = line_data(ivol*4+2);
                        end
                    end
                    
                    % Fill structure
                    obj.initial_guess.ric = ric;
                    obj.initial_guess.zis = zis;
                    obj.initial_guess.ris = ris;
                    obj.initial_guess.zic = zic;
                    
                case 1 % Henneberg representation
                                  
                    rhoi = zeros(2*obj.ntor+1, obj.mpol+1, obj.nvol);
                    bin  = zeros(obj.ntor+1, obj.nvol);
                    r0ic = zeros(obj.ntor+1, obj.nvol);
                    z0is = zeros(obj.ntor+1, obj.nvol);

                    for iline=1:nlines

                        line_data = str2num( initial_guess_str{iline} );

                        m = line_data(1); n = line_data(2);
                        im = m+1;
                        in = n+obj.ntor+1;

                        for ivol=1:obj.nvol

                            if n>=0 && m==0
                                bin( n+1, ivol) = line_data(ivol*4-1);
                                r0ic(n+1, ivol) = line_data(ivol*4  );
                                z0is(n+1, ivol) = line_data(ivol*4+1);
                            end

                            if m>0
                                rhoi(in, im, ivol ) = line_data( ivol*4+2 );
                            end
                        end

                    end


                    obj.initial_guess.rhoi = rhoi;
                    obj.initial_guess.bin  = bin;
                    obj.initial_guess.r0ic = r0ic;
                    obj.initial_guess.z0is = z0is;
                    
                otherwise
                    error('Invalid lboundary!')
            end
            
        end
        
        
        function obj = set_boundary_from_namelist( obj, filename, boundary )
           %
           % SET_BOUNDARY_FROM_NAMELIST( FILENAME )
           % ======================================
           %
           % Read the plasma or computational boundary from another SPEC
           % namelist and use it for the current instance
           %
           % INPUT
           % -----
           %   -filename: Filename of the other SPEC input namelist
           %   -boundary: PB for Plasma Boundary or CB for computational
           %              boundary
           %
           % OUTPUT
           % ------
           %   -obj: Updated output of SPEC Namelist
           
           
           % Read input SPEC Namelist
           nm = SPEC_Namelist( filename );
           
           % Get relevant boundary
           switch boundary
               case 'PB'
                   Remn = nm.physicslist.rbc;
                   Romn = nm.physicslist.rbs;
                   Zemn = nm.physicslist.zbc;
                   Zomn = nm.physicslist.zbs;
               case 'CB'
                   Remn = nm.physicslist.rwc;
                   Romn = nm.physicslist.rws;
                   Zemn = nm.physicslist.zwc;
                   Zomn = nm.physicslist.zws;
               otherwise
                   error('InputError: Invalid boundary')
           end
           
           % Check size
           if any(size(Remn)~=size(Romn))
               error('Size mismatch')
           end
           if any(size(Remn)~=size(Zemn))
               error('Size mismatch')
           end
           if any(size(Remn)~=size(Zomn))
               error('Size mismatch')
           end
           
           % Get Fourier resolution
           s = size(Remn);
           ntor_in = (s(1)-1) / 2.0;
           mpol_in =  s(2)-1;
           
           % Change FOurier resolution to the largest one between obj and
           % nm
           if mpol_in>obj.mpol
               obj = obj.change_fourier_resolution( mpol_in, obj.ntor );
           elseif mpol_in<obj.mpol
               nm = nm.change_fourier_resolution( obj.mpol, ntor_in );
               
               % nm poloidal resolution is now mpol_in
               mpol_in = obj.mpol;
           end
           
           if ntor_in>obj.ntor
               obj = obj.change_fourier_resolution( obj.mpol, ntor_in );
           elseif ntor_in<obj.ntor
               nm = nm.change_fourier_resolution( mpol_in, obj.ntor );
           end
           
           % Check that the size is correct
           if any(size(obj.physicslist.rbc)~=size(nm.physicslist.rbc))
               error('Size mismatch !')
           end
           
           % Overwrite relevant boundary in obj
           switch boundary
               case 'PB'
                   obj.physicslist.rbc = nm.physicslist.rbc;
                   obj.physicslist.rbs = nm.physicslist.rbs;
                   obj.physicslist.zbc = nm.physicslist.zbc;
                   obj.physicslist.zbs = nm.physicslist.zbs;
               case 'CB'
                   obj.physicslist.rwc = nm.physicslist.rwc;
                   obj.physicslist.rws = nm.physicslist.rws;
                   obj.physicslist.zwc = nm.physicslist.zwc;
                   obj.physicslist.zws = nm.physicslist.zws;
           end
        end
        
        % =================================================================
        % Getters
        function out = get_fourier_harmonics( obj, field, m, n, varargin )
            %
            % GET_FOUriEr_HArMONICS( FIELD, M, N, (IVOL) )
            % ============================================
            %
            % returns the fourier harmonics associated to the field,
            % poloidal mode number m and toroidal mode number n
            %
            % INPUTS
            % ------
            %   -field: 'rbc', 'zws', 'ric', ...
            %   -m    : Poloidal mode number
            %   -n    : Toroidal mode number
            %   -(ivol): Optional argument, required for fourier harmonics
            %            of the initial guess
            
            
            
            if m>obj.mpol
                warning('M is larger than the Fourier resolution')
                out = 0;
                return
            end 
            
            if abs(n)>obj.ntor
                warning('M is larger than the Fourier resolution')
                out = 0;
                return
            end
            
            % Need to make the difference between initial guess harmonics
            % or harmonics from the physics list. rhoi, bin, R0ic and Z0is
            % are related to the Henneberg representation
            if     strcmp(field, 'ric') || strcmp(field, 'ris') ...
                || strcmp(field, 'zic') || strcmp(field, 'zis') ...
                || strcmp(field, 'rhoi') || strcmp(field, 'bin') ...
                || strcmp(field, 'r0ic') || strcmp(field, 'z0is')
                
                if isempty(obj.initial_guess)
                    error('No initial guess available')
                end
                
                ivol = varargin{1};
            
                if ivol<1 || ivol>obj.nvol
                    error('Invalid ivol')
                end                   
                
                cat = 'initial_guess';
                
            else
                cat = 'physicslist';
                ivol = 1;
                
            end
                
            % Read relevant harmonic
            out = obj.(cat).(field)(n+obj.ntor+1, m+1, ivol);
            
        end
        
        % =================================================================
        % Setter
        function obj = set_harmonics_to_zero( obj, field )
            %
            % SET_HARMONICS_TO_ZERO( FIELD )
            % ==============================
            %
            % Set all Fourier harmonics of the field to zero
            %
            % INPUT
            % -----
            %   -FIELD: Field to set harmonics to zero
            %
            % OUTPUT
            % ------
            %   -OBJ: Updated instance of SPEC_Namelist\
            %
            
            
            if     strcmp(field, 'ric') || strcmp(field, 'ris') ...
                || strcmp(field, 'zic') || strcmp(field, 'zis') ...
                || strcmp(field, 'rhoi') || strcmp(field, 'bin') ...
                || strcmp(field, 'r0ic') || strcmp(field, 'z0is')
            
                cat = 'initial_guess';
            else
                cat = 'physicslist';
            end
            
            obj.(cat).(field) = zeros( obj.array_size );
            
        end
        
        
        
        function obj = set_fourier_harmonics( obj, field, im, in, value, varargin )
            %
            % SET_FOURiER_HARMONICS( FIELD, M, N, (IVOL) )
            % ============================================
            %
            % Sets the fourier harmonics associated to the field,
            % poloidal mode number m and toroidal mode number n
            %
            % INPUTS
            % ------
            %   -field: 'rbc', 'zws', 'ric', ...
            %   -im   : Poloidal mode number, size 1xmn
            %   -in   : Toroidal mode number, size 1xmn
            %   -value: Value for the mode, size 1xmn
            %   -(ivol): Optional argument, required for fourier harmonics
            %            of the initial guess
            
            
            % First check the input size
            mn = length(im);
            if length(in)~=mn
                error('The array in has not the same length as im')
            end
            if length(value)~=mn
                error('The array value has not the same length as im')
            end
            
            % If some modes are greater than the actual resolution,
            % increase the resolution accordingly
            mpol_in = max(im);
            if mpol_in>obj.mpol
                warning('Poloidal resolution has to be increased ...')
                obj = obj.change_fourier_resolution( mpol_in, obj.ntor );
            end 

            ntor_in = max(abs(in));
            if ntor_in>obj.ntor
                warning('Toroidal resolution has to be increased ...')
                obj = obj.change_fourier_resolution( obj.mpol, ntor_in );
            end
            
            % Loop over the input harmonics
            for imn=1:mn
                m = im(imn);
                n = in(imn);

                
                % Check if the quantity is located in the initial guess
                % structure or in the physicslist structure.
                if     strcmp(field, 'ric') || strcmp(field, 'ris') ...
                    || strcmp(field, 'zic') || strcmp(field, 'zis') ...
                    || strcmp(field, 'rhoi') || strcmp(field, 'bin') ...
                    || strcmp(field, 'r0ic') || strcmp(field, 'z0is')

                    if isempty(obj.initial_guess)
                        warning('No initial guess available... filling with zeros')
                    
                        % If the initial guess is empty (i.e. no initial guess
                        % are available), then we create one filled with zeros
                        if obj.lboundary == 0
                            obj.initial_guess.ric = zeros( obj.array_size(1), obj.array_size(2), obj.nvol );
                            obj.initial_guess.ris = zeros( obj.array_size(1), obj.array_size(2), obj.nvol );
                            obj.initial_guess.zic = zeros( obj.array_size(1), obj.array_size(2), obj.nvol );
                            obj.initial_guess.zis = zeros( obj.array_size(1), obj.array_size(2), obj.nvol );                    
                        else
                            obj.initial_guess.rhoi = zeros( obj.array_size(1), obj.array_size(2), obj.nvol );
                            obj.initial_guess.bin  = zeros( obj.ntor+1, obj.nvol );
                            obj.initial_guess.r0ic = zeros( obj.ntor+1, obj.nvol );
                            obj.initial_guess.z0is = zeros( obj.ntor+1, obj.nvol );
                        end
                    end

                    ivol = varargin{1};

                    if ivol<1 || ivol>obj.nvol
                        error('Invalid ivol')
                    end                   

                    obj.initial_guess.(field)(n+obj.ntor+1, m+1, ivol) = value(imn);

                else
                    
                    obj.physicslist.(field)(n+obj.ntor+1, m+1) = value(imn);
                    

                end
            end
            
        end
        
        function obj = truncate_fourier_series( obj, mpol, ntor )
            %
            % TRUNCATE_FOURiER_SERiES( mpol, ntor )
            % =====================================
            %
            % Truncates all spectral quantities to the requested poloidal
            % and toroidal resolution This can also be used to increase
            % the Fourier resolution.
            %
            % INPUTS
            % ------
            %   -mpol: Poloidal resolution
            %   -ntor: Toroidal resolution
            %
            % OUTPUT
            % ------
            %   -OBJ: Truncated instance of SPEC_Namelist
            %
           
            obj = obj.change_fourier_resolution( mpol, ntor );
            
        end
        
        function obj = change_boundary_representation( obj, new_lboundary )
            % 
            % CHANGE_BOUNDArY_rEPrESENTATION( NEW_lbOUNDArY )
            % ===============================================
            %
            % Changes from the hudson representation to the henneberg's one
            % and vice versa.
            %
            % INPUT
            % -----
            %   -new_lboundary: New value for lboundary
            
            
            if obj.lboundary==0 && new_lboundary==1
                
                obj.physicslist.rhomn = zeros(2*obj.ntor+1, obj.mpol);
                obj.physicslist.shift.rhomn = [obj.ntor+1, 1];
                obj.physicslist.bn    = zeros(obj.ntor+1, 1);
                obj.physicslist.r0c   = zeros(obj.ntor+1, 1);
                obj.physicslist.z0s   = zeros(obj.ntor+1, 1);
                
                obj.physicslist = rmfield( obj.physicslist, 'rbc' );
                obj.physicslist = rmfield( obj.physicslist, 'rbs' );
                obj.physicslist = rmfield( obj.physicslist, 'zbc' );
                obj.physicslist = rmfield( obj.physicslist, 'zbs' );
                
            elseif obj.lboundary==1 && new_lboundary==0
                obj.physicslist.rbc = zeros(2*obj.ntor+1, obj.mpol);
                obj.physicslist.rbs = zeros(2*obj.ntor+1, obj.mpol);
                obj.physicslist.zbc = zeros(2*obj.ntor+1, obj.mpol);
                obj.physicslist.zbs = zeros(2*obj.ntor+1, obj.mpol);
                
                obj.physicslist = rmfield( obj.physicslist, 'rhomn' );
                obj.physicslist = rmfield( obj.physicslist, 'bn' );
                obj.physicslist = rmfield( obj.physicslist, 'r0c' );
                obj.physicslist = rmfield( obj.physicslist, 'z0s' );
                
                obj.physicslist.shift.rbc = [obj.ntor+1, 1];
                obj.physicslist.shift.rbs = [obj.ntor+1, 1];
                obj.physicslist.shift.zbc = [obj.ntor+1, 1];
                obj.physicslist.shift.zbs = [obj.ntor+1, 1];
                
            end
            
            obj.physicslist.lboundary = new_lboundary;
            obj.lboundary = new_lboundary;
            
            
                       
            
        end
        
        % =================================================================
        % Plotters
        function plot_plasma_boundary( obj, nt, phi, newfig, varargin )
           %
           % PLOT_PLASMA_BOUNDArY
           % ====================
           %
           % Plot the boundary given by rbc, zbs, rbs, zbc
           %
           % INPUTS
           % ------
           %   -NT:     Number of poloidal points
           %   -PHI:    Toroidal angle
           %   -NEWFIG: =0: plot on gca
           %            =1: plot on a new figure
           %            =2: erase and plot on gca
           %   -varargin: Any input you could give to plot()
           %
           %
           
           if obj.lboundary==0
               rbc = obj.physicslist.rbc;
               rbs = obj.physicslist.rbs;
               zbc = obj.physicslist.zbc;
               zbs = obj.physicslist.zbs;

               obj.plot_surface_lb0( rbc, zbs, rbs, zbc, nt, phi, newfig, varargin{:} )
           else
               rhomn = obj.physicslist.rhomn;
               bn    = obj.physicslist.bn;
               r0c   = obj.physicslist.r0c;
               z0s   = obj.physicslist.z0s;
               
               obj.plot_surface_lb1( rhomn, bn, r0c, z0s, nt, phi, newfig, varargin{:} )
           end
        end
        
        function plot_computational_boundary( obj, nt, phi, VorB, newfig, varargin )
           %
           % PLOT_PLASMA_BOUNDArY
           % ====================
           %
           % Plot the boundary given by rbc, zbs, rbs, zbc
           %
           % INPUTS
           % ------
           %   -NT:     Number of poloidal points
           %   -PHI:    Toroidal angle
           %   -VorB:   ='V'
           %            ='B'
           %            ='F'
           %            ='N'
           %   -NEWFIG: =0: plot on gca
           %            =1: plot on a new figure
           %            =2: erase and plot on gca
           %   -varargin: Any input you could give to plot()
           %
           %
           
           if obj.lboundary == 0
               rwc = obj.physicslist.rwc;
               rws = obj.physicslist.rws;
               zwc = obj.physicslist.zwc;
               zws = obj.physicslist.zws;

               obj.plot_surface_lb0( rwc, zws, rws, zwc, nt, phi, newfig, varargin{:} )
           else
               error('Not implemented')
           end

           if VorB~='N'
               newfig = 0;
               obj.plot_normal_field( floor(nt/10), phi, VorB, newfig )
           end
        end
           
        function plot_initial_guess( obj, nt, phi, newfig, varargin )
           %
           % PLOT_PLASMA_BOUNDArY( NT, PHI, NEWFIG, VArArGIN )
           % =================================================
           %
           % Plot the boundary given by rbc, zbs, rbs, zbc
           %
           % INPUTS
           % ------
           %   -NT: Number of poloidal points
           %   -PHI: Toroidal angle
           %   -NEWFIG: =0: plot on gca
           %            =1: plot on a new figure
           %            =2: erase and plot on gca
           %   -varargin: Any input you could give to plot()
           %
           %
           
           if isempty(obj.initial_guess)
               error('No initial guess is provided')
           end
           
           if obj.lboundary == 0
               for ivol=1:obj.nvol
                   ric = obj.initial_guess.ric(:,:,ivol);
                   ris = obj.initial_guess.ris(:,:,ivol);
                   zic = obj.initial_guess.zic(:,:,ivol);
                   zis = obj.initial_guess.zis(:,:,ivol);

                   if ivol>=2
                       newfig=0;
                   end

                   obj.plot_surface_lb0( ric, zis, ris, zic, nt, phi, newfig, varargin{:} )
               end
           else
               for ivol=1:obj.nvol
                  rhoi = obj.initial_guess.rhoi(:,:,ivol);
                  bin  = obj.initial_guess.bin(:, ivol);
                  r0ic = obj.initial_guess.r0ic(:, ivol);
                  z0is = obj.initial_guess.z0is(:, ivol);
                  
                  if ivol>=2
                      newfig=0;
                  end
                  
                  obj.plot_surface_lb1( rhoi, bin, r0ic, z0is, nt, phi, newfig, varargin{:} )
               end
           end
        end
        
        % =================================================================
        % Write method
        function write_input_file(obj, filename )
            %
            % WriTE_INPUT_FILE( filename )
            % ============================
            %
            % Write namelist in an input file. Be careful: if the file filename
            % provided as input exist, it will be overwritten!
            %
            % INPUT
            % -----
            %   - filename: path where to save the input file.
            %

            % Set minimal toroidal resolution to one, otherwise writting
            % routine does not detect arrays. This is not ideal and should
            % be fixed
            if obj.ntor==0
                obj = obj.change_fourier_resolution( obj.mpol, 1 );
            end
            
            % Create a structure with the different lists
            nlists = length(obj.lists);
            S = struct;
            for ii=1:nlists
                S.(obj.lists{ii}) = obj.(obj.lists{ii});
            end
            
            % Create a shift quantity - this tells the writing routine how
            % much each index has to be shifted. In MATLAB, indices start
            % at 1, while in the FORTRAN Namelist, we want to write
            S.shift.rbc = [obj.ntor+1, 1];
            S.shift.rbs = [obj.ntor+1, 1];
            S.shift.zbc = [obj.ntor+1, 1];
            S.shift.zbs = [obj.ntor+1, 1];
            S.shift.rwc = [obj.ntor+1, 1];
            S.shift.rws = [obj.ntor+1, 1];
            S.shift.zws = [obj.ntor+1, 1];
            S.shift.zwc = [obj.ntor+1, 1];
            S.shift.vnc = [obj.ntor+1, 1];
            S.shift.vns = [obj.ntor+1, 1];
            S.shift.bnc = [obj.ntor+1, 1];
            S.shift.bns = [obj.ntor+1, 1];
            S.shift.rhomn = [obj.ntor+1, 1];

            % remove unnecessary fields
            if obj.lboundary == 0
                if isfield( S.physicslist, 'rhomn' )
                    S.physicslist = rmfield(S.physicslist, 'rhomn');
                end
                if isfield( S.physicslist, 'bn' )
                    S.physicslist = rmfield(S.physicslist, 'bn');
                end
                if isfield( S.physicslist, 'r0c' )
                    S.physicslist = rmfield(S.physicslist, 'r0c');
                end
                if isfield( S.physicslist, 'z0s' )
                    S.physicslist = rmfield(S.physicslist, 'z0s');
                end
                
            else % lboundary==1
                if isfield( S.physicslist, 'rbc' )
                    S.physicslist = rmfield(S.physicslist, 'rbc');
                end
                if isfield( S.physicslist, 'rbs' )
                    S.physicslist = rmfield(S.physicslist, 'rbs');
                end
                if isfield( S.physicslist, 'zbc' )
                    S.physicslist = rmfield(S.physicslist, 'zbc');
                end
                if isfield( S.physicslist, 'zbs' )
                    S.physicslist = rmfield(S.physicslist, 'zbs');
                end            
                
            end
            
            if obj.physicslist.lfreebound==0 % Then no need to freeboundary info
                if isfield( S.physicslist, 'rwc' )
                    S.physicslist = rmfield(S.physicslist, 'rwc');
                end
                if isfield( S.physicslist, 'rws' )
                    S.physicslist = rmfield(S.physicslist, 'rws');
                end
                if isfield( S.physicslist, 'zwc' )
                    S.physicslist = rmfield(S.physicslist, 'zwc');
                end
                if isfield( S.physicslist, 'zws' )
                    S.physicslist = rmfield(S.physicslist, 'zws');
                end 
                if isfield( S.physicslist, 'vnc' )
                    S.physicslist = rmfield(S.physicslist, 'vnc');
                end
                if isfield( S.physicslist, 'vns' )
                    S.physicslist = rmfield(S.physicslist, 'vns');
                end
                if isfield( S.physicslist, 'bnc' )
                    S.physicslist = rmfield(S.physicslist, 'bnc');
                end
                if isfield( S.physicslist, 'bns' )
                    S.physicslist = rmfield(S.physicslist, 'bns');
                end 
                
            end
            
            
            % Build initial guess strings
            if obj.physicslist.lboundary==0        
                if ~isempty(obj.initial_guess)
                    s = size(obj.initial_guess.ric);
                    initialguess = cell(1, s(1)*s(2));
                    iline=0;
                    for ii=1:s(1)
                        for jj=1:s(2)
                            iline = iline+1;
                            mm = jj-1;
                            nn = ii-obj.ntor-1;
                            
                            initialguess{iline} = sprintf( '%i   %i   ', mm, nn );
                            for ivol=1:s(3)
                               initialguess{iline} = sprintf( '%s   %0.12E   %0.12E   %0.12E   %0.12E', initialguess{iline}, ...
                                                            obj.initial_guess.ric(ii,jj,ivol), ...
                                                            obj.initial_guess.zis(ii,jj,ivol), ...
                                                            obj.initial_guess.ris(ii,jj,ivol), ...
                                                            obj.initial_guess.zic(ii,jj,ivol)     );
                            end
                        end
                    end
                else
                    initialguess = cell(0);
                end

            else

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
                                                     obj.initial_guess.r0ic(ii, ivol ), ...
                                                     obj.initial_guess.z0is(ii, ivol ), 0.0 ); 
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

            % Call namelist writer
            write_namelist( S, filename, initialguess );
       

        end
    end
    
    methods (Access=private)
        
        function obj = initialize_structure( obj )
            %
            % INITIALIzE_STrUCTUrE( OBJ )
            % ===========================
            %
            % Use to check that all required inputs are correctly set;
            % check that the sizes of arrays are correct. raise errors in
            % case crucial informations is missing.
            %
            % 
                       
            % Fill some important inputs
            if ~isfield(obj.physicslist, 'nvol')
                error('Missing nvol')
            end
            if ~isfield(obj.physicslist, 'lfreebound')
                error('Missing lfreebound')
            end
            
            obj.Mvol = obj.physicslist.nvol + obj.physicslist.lfreebound;
            obj.nvol = obj.physicslist.nvol;
            
            if isfield(obj.physicslist, 'lboundary')
                obj.lboundary = obj.physicslist.lboundary;
            else
                warning('lboundary not provided. Setting with 0...')
                obj.lboundary = 0;
            end
            
            % PHYSICSLIST
            % -----------
            % lrad
            if ~isfield(obj.physicslist, 'lrad') && obj.verbose
               warning('Missing lrad. Filling with 4...')
               obj.physicslist.lrad = ones(1,obj.Mvol) * 4;
            else
                % Fill potential missing elements with 4s
               obj.physicslist.lrad(end+1:obj.Mvol) = 4;                
            end
            
            % tflux
            if ~isfield(obj.physicslist, 'tflux') && obj.verbose
               warning('Missing tflux. Filling equal radial distances')
               obj.physicslist.tflux = (1:obj.Mvol).^2;
            else
                if length(obj.physicslist.tflux)~=obj.Mvol
                    error('Invalid number of tflux elements')
                end                
            end
            
            % pflux
            if ~isfield(obj.physicslist, 'pflux')
               warning('Missing pflux. Filling with 0...')
               obj.physicslist.pflux = zeros(1,obj.Mvol);
            else
                % Fill potential missing elements with zeros
               obj.physicslist.pflux(end+1:obj.Mvol) = 0;                
            end
            
            % helicity
            if ~isfield(obj.physicslist, 'helicity')
               warning('Missing lrad. Filling with 0...')
               obj.physicslist.helicity = zeros(1,obj.Mvol);
            else
                % Fill potential missing elements with zeros
               obj.physicslist.helicity(end+1:obj.Mvol) = 0;                
            end
            
            % pscale
            if ~isfield(obj.physicslist, 'pscale')
                warning('Missing pscale. Setting to zero...')
                obj.physicslist.pscale = 0.0;
            end
            
            % ladiabatic
            if ~isfield(obj.physicslist, 'ladiabatic')
                warning('Missing ladiabatic. Setting to zero...')
                obj.physicslist.ladiabatic = 0.0;
            end
            
            % Pressure
            if ~isfield(obj.physicslist, 'pressure')
               warning('Missing pressure. Filling with 0...')
               obj.physicslist.pressure = zeros(1,obj.Mvol);
            else
                % Fill potential missing elements with zeros
               obj.physicslist.pressure(end+1:obj.Mvol) = 0;                
            end
            
            % Adiabatic
            if obj.physicslist.ladiabatic==1
                if ~isfield(obj.physicslist, 'adiabatic')
                   warning('Missing adiabatic. Filling with 0...')
                   obj.physicslist.adiabatic = zeros(1,obj.Mvol);
                else
                    % Fill potential missing elements with zeros
                   obj.physicslist.adiabatic(end+1:obj.Mvol) = 0;                
                end
            end
            
            % mu
            if ~isfield(obj.physicslist, 'mu')
               warning('Missing mu. Filling with 0...')
               obj.physicslist.mu = zeros(1,obj.Mvol);
            else
                % Fill potential missing elements with zeros
               obj.physicslist.mu(end+1:obj.Mvol) = 0;                
            end
            
            % lconstraint
            if ~isfield(obj.physicslist, 'lconstraint')
               warning('Missing lconstraint. Setting to 0...')
               obj.physicslist.lconstraint = 0;              
            end
            
            if ~any(obj.physicslist.lconstraint==[0,1,2,3])
                error('Invalid lconstraint')
            end
            
            
            % Ivolume, Isurf
            if obj.physicslist.lconstraint==3
                if ~isfield(obj.physicslist, 'Ivolume')
                   warning('Missing Ivolume. Filling with 0...')
                   obj.physicslist.Ivolume = zeros(1,obj.Mvol);
                else
                    % Fill potential missing elements with zeros
                   obj.physicslist.Ivolume(end+1:obj.Mvol) = 0;                
                end
                if ~isfield(obj.physicslist, 'Isurf')
                   warning('Missing Isurf. Filling with 0...')
                   obj.physicslist.Isurf = zeros(1,obj.Mvol);
                else
                    % Fill potential missing elements with zeros
                   obj.physicslist.Isurf(end+1:obj.Mvol) = 0;                
                end                
            end
            
            if obj.physicslist.lconstraint==1
                if ~isfield(obj.physicslist, 'iota')
                   warning('Missing iota. Filling with sqrt(2)...')
                   obj.physicslist.iota = sqrt(2)*ones(1,obj.Mvol);
                else
                    % Fill potential missing elements with zeros
                   obj.physicslist.iota(end+1:obj.Mvol) = sqrt(2);                
                end
                if ~isfield(obj.physicslist, 'oita')
                   warning('Missing oita. Filling with sqrt(2)...')
                   obj.physicslist.oita = sqrt(2)*ones(1,obj.Mvol);
                else
                    % Fill potential missing elements with zeros
                   obj.physicslist.oita(end+1:obj.Mvol) = sqrt(2);                
                end                    
            end
            
            % mupftol
            if ~isfield(obj.physicslist, 'mupftol')
               warning('Missing mupftol. Setting to 1E-12...')
               obj.physicslist.mupftol = 1E-12;              
            end
            
            % mupfits
            if ~isfield(obj.physicslist, 'mupfits')
               warning('Missing mupfits. Setting to 128...')
               obj.physicslist.mupfits = 128;              
            end
            
            % Check geometry
            if ~isfield(obj.physicslist, 'mpol')
                error('Missing mpol information')
            end
            if ~isfield(obj.physicslist, 'ntor')
                error('Missing ntor information')
            end
            if ~isfield(obj.physicslist, 'lboundary')
               warning('Missing lboundary. Setting to zero')
               obj.physicslist.lboundary = 0;
            end
            
            mpol_in = obj.physicslist.mpol;
            ntor_in = obj.physicslist.ntor;
            
            if obj.lboundary==0
                if ~isfield(obj.physicslist, 'rbc')
                    obj.physicslist.shift.rbc = [ntor_in+1, 1];
                    obj.physicslist.rbc = zeros(2*ntor_in+1, mpol_in);
                end
                if ~isfield(obj.physicslist, 'rbs')
                    obj.physicslist.shift.rbs = [ntor_in+1, 1];
                    obj.physicslist.rbs = zeros(2*ntor_in+1, mpol_in);
                end
                if ~isfield(obj.physicslist, 'zbc')
                    obj.physicslist.shift.zbc = [ntor_in+1, 1];
                    obj.physicslist.zbc = zeros(2*ntor_in+1, mpol_in);
                end
                if ~isfield(obj.physicslist, 'zbs')
                    obj.physicslist.shift.zbs = [ntor_in+1, 1];
                    obj.physicslist.zbs = zeros(2*ntor_in+1, mpol_in);
                end
                
                % Check that sizes are consistent with each others
                if any(size(obj.physicslist.rbc)~=size(obj.physicslist.rbs))
                    error('Size mismatch between rbc and rbs')
                end
                if any(size(obj.physicslist.rbc)~=size(obj.physicslist.zbc))
                    error('Size mismatch between rbc and zbc')
                end
                if any(size(obj.physicslist.rbc)~=size(obj.physicslist.zbs))
                    error('Size mismatch between rbc and zbs')
                end
                
                
            else %lboundary==1
                if ~isfield( obj.physicslist, 'bn' )
                    obj.physicslist.bn = zeros( ntor_in+1, 1 );
                end
                
                % Now fill missing elements with zeros
                obj.physicslist.bn(end+1:ntor_in+1) = 0.0;
                
                if ~isfield( obj.physicslist, 'r0c' )
                    obj.physicslist.r0c = zeros( ntor_in+1, 1 );
                end
                
                % Now fill missing elements with zeros
                obj.physicslist.r0c(end+1:ntor_in+1) = 0.0;
                
                if ~isfield( obj.physicslist, 'z0s' )
                    obj.physicslist.z0s = zeros( ntor_in+1, 1 );
                end
                
                % Now fill missing elements with zeros
                obj.physicslist.z0s(end+1:ntor_in+1) = 0.0;
                
                if ~isfield( obj.physicslist, 'rhomn' )
                    obj.physicslist.shift.rhomn = [ntor_in+1, 1];
                    obj.physicslist.rhomn = zeros( 2*ntor_in+1, mpol_in+1 );
                end
                                
                % Check sizes                
                if any(length(obj.physicslist.bn)~=ntor_in+1)
                    obj.physicslist.bn(end+1:ntor_in+1) = 0;
                    obj.physicslist.bn = obj.physicslist.bn(ntor_in+1);
                end
                if any(length(obj.physicslist.r0c)~=ntor_in+1)
                    obj.physicslist.r0c(end+1:ntor_in+1) = 0;
                    obj.physicslist.r0c = obj.physicslist.r0c(ntor_in+1);
                end
                if any(length(obj.physicslist.z0s)~=ntor_in+1)
                    obj.physicslist.z0s(end+1:ntor_in+1) = 0;
                    obj.physicslist.z0s = obj.physicslist.z0s(ntor_in+1);
                end
                
                
            end
            
             

            if ~isfield(obj.physicslist, 'rwc')
                obj.physicslist.shift.rwc = [ntor_in+1, 1];
                obj.physicslist.rwc = zeros(2*ntor_in+1, mpol_in);
            end
            if ~isfield(obj.physicslist, 'rws')
                obj.physicslist.shift.rws = [ntor_in+1, 1];
                obj.physicslist.rws = zeros(2*ntor_in+1, mpol_in);
            end
            if ~isfield(obj.physicslist, 'zwc')
                obj.physicslist.shift.zwc = [ntor_in+1, 1];
                obj.physicslist.zwc = zeros(2*ntor_in+1, mpol_in);
            end
            if ~isfield(obj.physicslist, 'zws')
                obj.physicslist.shift.zws = [ntor_in+1, 1];
                obj.physicslist.zws = zeros(2*ntor_in+1, mpol_in);
            end

            % Check that sizes are consistent with each others
            if any(size(obj.physicslist.rwc)~=size(obj.physicslist.rws))
                error('Size mismatch between rwc and rws')
            end
            if any(size(obj.physicslist.rwc)~=size(obj.physicslist.zwc))
                error('Size mismatch between rwc and zwc')
            end
            if any(size(obj.physicslist.rwc)~=size(obj.physicslist.zws))
                error('Size mismatch between rwc and zws')
            end            
                
                
            
            if ~isfield(obj.physicslist, 'vnc')
                obj.physicslist.shift.vnc = [ntor_in+1, 1];
                obj.physicslist.vnc = zeros(2*ntor_in+1, mpol_in);
            end
            if ~isfield(obj.physicslist, 'vns')
                obj.physicslist.shift.vns = [ntor_in+1, 1];
                obj.physicslist.vns = zeros(2*ntor_in+1, mpol_in);
            end
            if ~isfield(obj.physicslist, 'bnc')
                obj.physicslist.shift.bnc = [ntor_in+1, 1];
                obj.physicslist.bnc = zeros(2*ntor_in+1, mpol_in);
            end
            if ~isfield(obj.physicslist, 'bns')
                obj.physicslist.shift.bns = [ntor_in+1, 1];
                obj.physicslist.bns = zeros(2*ntor_in+1, mpol_in);
            end
            
            % Check that sizes are consistent with each others
            if any(size(obj.physicslist.vnc)~=size(obj.physicslist.vns))
                error('Size mismatch between vnc and vns')
            end
            if any(size(obj.physicslist.bnc)~=size(obj.physicslist.bnc))
                error('Size mismatch between bnc and bns')
            end
            
            % DIAGNOSTICSLIST
            % ---------------
            if ~isfield(obj.diagnosticslist, 'nppts')
                warning('Missing nppts. Setting to zero')
                obj.diagnosticslist.nppts = 0;
            end
            
            if ~isfield(obj.diagnosticslist, 'nptrj')
                warning('Missing nptrj. Setting to zero')
                obj.diagnosticslist.nptrj = 0;
            end            
        end
        
        
        function obj = change_fourier_resolution( obj, mpol_new, ntor_new )
            %
            % CHANGE_FOUriEr_rESOLUTION( mpol_NEW, ntor_NEW )
            % ===============================================
            %
            % Change inner Fourier resolution of an instance of
            % SPEC_Namelist.
            %
            % INPUTS
            % ------
            %   -mpol_new: New poloidal resolution
            %   -ntor_new: New toroidal resolution
            %
            % OUTPUT
            % ------
            %   -obj: Updated instance of SPEC_Namelist
            %
       
            if mpol_new<1
                error('InputError, mpol_new should be larger than 0')
            end
            if ntor_new<0
                error('InputError, ntor_new should be larger or equal to zero')
            end
            
            obj.mpol = mpol_new;
            obj.ntor = ntor_new;
            obj.array_size = [2*obj.ntor+1, obj.mpol+1];
            
            % Check that all arrays have the same size; otherwise, fill
            % with zeros the missing elements
            if obj.lboundary == 0
                if any(size(obj.physicslist.rbc)~=obj.array_size)
                    obj = obj.reshape_array( 'rbc' );
                    obj = obj.reshape_array( 'rbs' );
                    obj = obj.reshape_array( 'zbc' );
                    obj = obj.reshape_array( 'zbs' );
                end

                if any(size(obj.physicslist.rwc)~=obj.array_size)
                    obj = obj.reshape_array( 'rwc' );
                    obj = obj.reshape_array( 'rws' );
                    obj = obj.reshape_array( 'zwc' );
                    obj = obj.reshape_array( 'zws' );
                end
            else
                if any(size(obj.physicslist.rhomn)~=obj.array_size)
                    obj = obj.reshape_array( 'rhomn' );
                end
                
                if any(size(obj.physicslist.bn)~=[obj.ntor+1, 1])
                    obj.physicslist.bn(end+1:obj.ntor+1) = 0;
                    obj.physicslist.bn = obj.physicslist.bn(1:obj.ntor+1);
                end
                if any(size(obj.physicslist.r0c)~=[obj.ntor+1, 1])
                    obj.physicslist.r0c(end+1:obj.ntor+1) = 0;
                    obj.physicslist.r0c = obj.physicslist.r0c(1:obj.ntor+1);
                end
                if any(size(obj.physicslist.z0s)~=[obj.ntor+1, 1])
                    obj.physicslist.z0s(end+1:obj.ntor+1) = 0;
                    obj.physicslist.z0s = obj.physicslist.z0s(1:obj.ntor+1);
                end
                
            end
                
            
            if any(size(obj.physicslist.vnc)~=obj.array_size)
                obj = obj.reshape_array( 'vnc' );
                obj = obj.reshape_array( 'vns' );
            end
            
            if any(size(obj.physicslist.bnc)~=obj.array_size)
                obj = obj.reshape_array( 'bnc' );
                obj = obj.reshape_array( 'bns' );           
            end
        end
        
        
        function obj = set_fourier_resolution( obj )
            %
            % SET_FOUriEr_rESOLUTION( OBJ )
            % =============================
            %
            % Set the obj internal Fourier resolution to the maximal value
            % required to store the given data; then, reshapees all arrays
            % to have the size obj.array_size.
            %
            % Check beforehand that the resolution of rbc is the same as
            % rbs, zbc and zbs. Do something similar for vnc, vns, bnc, bns
            % and rwc, rws, zws, zwc.
            %
           
            
            % Find largest Fourier resolution in the input file
            mpol_in = obj.physicslist.mpol;
            ntor_in = obj.physicslist.ntor;
            
            if obj.lboundary == 0
                s_bc = size(obj.physicslist.rbc);
                shift = obj.physicslist.shift.rbc(1);
                mpol_bc = s_bc(2)-1       ;
                ntor_bc = max([abs(1-shift), s_bc(1)-shift]);
            else
                s_bc = size(obj.physicslist.rhomn);
                shift = obj.physicslist.shift.rhomn(1);
                mpol_bc = s_bc(2)-1;
                
                lbn = length( obj.physicslist.bn  );
                lrc = length( obj.physicslist.r0c );
                lzs = length( obj.physicslist.z0s );
                ntor_bc = max([abs(1-shift), s_bc(1)-shift, lbn, lrc, lzs]);
            end
            
            s_wc = size(obj.physicslist.rwc);
            shift = obj.physicslist.shift.rwc(1);
            mpol_wc = s_wc(2)-1       ;
            ntor_wc = max([abs(1-shift), s_wc(1)-shift]);
            
            s_vb = size(obj.physicslist.vnc);
            shift = obj.physicslist.shift.vnc(1);
            mpol_vb = s_vb(2)-1       ;
            ntor_vb = max([abs(1-shift), s_vb(1)-shift]);
            
            mpol_new = max([mpol_in, mpol_bc, mpol_wc, mpol_vb]);
            ntor_new = max([ntor_in, ntor_bc, ntor_wc, ntor_vb]);
            obj = obj.change_fourier_resolution( mpol_new, ntor_new );
        end
        
        
        function obj = reshape_array( obj, field )
            %
            % rESHAPE_ArrAY( OBJ, FIELD )
            % =======================================
            %
            % reshape the input array in a an array of size obj.array_size
            % and fills missing elements with zeros
            %
            % INPUTS
            % ------
            %   -field: Field, in obj.physicslist, to be modified
            %
            % OUTPUT
            % ------
            %   -obj: Updated instance of SPEC_Namelist
            %
            
            array = obj.physicslist.(field);
            shift = obj.physicslist.shift.(field);
            
            s = size(array);   
            
            new_array = zeros(obj.array_size);
            for ii=1:s(1)
                for jj=1:s(2)
                    nn = ii-shift(1);
                    mm = jj-shift(2);
                    
                    if abs(nn)>obj.ntor || mm>obj.mpol
                        continue
                    end
                    
                    new_array(nn+obj.ntor+1, mm+1) = array(ii, jj);
                end
            end
            
            obj.physicslist.(field) = new_array;
        end
        
        
        function plot_surface_lb0( obj, rmnc, zmns, rmns, zmnc, nt, phi, ...
                               newfig, varargin )
            %
            % PLOT_SUrFACE( rMNC, zMNS, rMNS, zMNC, NT, NEWFIG, VArArGIN )
            % ===========================================================
            %
            % Plot a surface parametrized by the standard representation
            %
            % INPUTS
            % ------
            %   -rmnc: Even Fourier modes of r, format (2*ntor+1, mpol+1)
            %   -zmns: Odd  Fourier modes of z, format (2*ntor+1, mpol+1)
            %   -rmns: Odd  Fourier modes of r, format (2*ntor+1, mpol+1)
            %   -zmnc: Even Fourier modes of z, format (2*ntor+1, mpol+1)
            %   -NT  : Number of poloidal points
            %   -PHI : Toroidal angle
            %   -newfig: =0: plot on gca
            %            =1: plot on a new figure
            %            =2: erase and plot on gca
            %   - varargin: Any input you could give to plot()
            
            switch newfig
                case 0
                    hold on
                case 1
                    figure('Color','w','Position',[200 200 900 700])
                    hold on
                case 2
                    hold off
                otherwise
                    error('InputError: invalid newfig')
            end
            
            s = size(rmnc);
            
            if any(s~=size(rmns))
                error('InputError: rmnc has not the same size as rmns')
            end
            if any(s~=size(zmns))
                error('InputError: rmnc has not the same size as zmns')
            end
            if any(s~=size(zmnc))
                error('InputError: rmnc has not the same size as zmnc')
            end    
            if nt<1
                error('InputError: nt should be larger than 1')
            end
            
            N = (s(1)-1) / 2.0;
            
            tarr = linspace( 0, 2*pi, nt );
            r    = zeros( 1, nt );
            z    = zeros( 1, nt );
            nfp  = double(obj.physicslist.nfp);
            
            
            for in=1:s(1)
                nn = in-1-N;
                for im=1:s(2)
                    mm = im-1;
                    
                    arg = mm*tarr - nn*nfp*phi;
                    
                    r = r + rmnc(in,im) * cos(arg) + rmns(in,im) * sin(arg);
                    z = z + zmnc(in,im) * cos(arg) + zmns(in,im) * sin(arg);
                end
            end
            
            plot( r, z, varargin{:} )
            
            axis equal
        end
        
        
        function plot_surface_lb1( obj, rhomn, bn, r0c, z0s, nt, phi, newfig, varargin )
           %
           % PLOT_SUrFACE_lb1( rHOMN, bn, r0C, z0S, NT, PHI, NEWFIG, VArArGIN )
           % ==================================================================
           %
           % Plots a surface using the Henneberg representation
           %
           % INPUTS
           % ------
           %   -rhomn:  rho_mn harmonics, format (2*ntor+1, mpol+1)
           %   -bn:     b_n harmonics, format (ntor+1, 1)
           %   -r0c:    r_0c harmonics, format (ntor+1, 1)
           %   -z0s:    z_0s harmonics, format (ntor+1, 1)
           %   -nt:     Number of poloidal points
           %   -phi:    Toroidal angle
           %   -newfig: =0: plot on gca
           %            =1: plot on a new figure
           %            =2: erase and plot on gca
           %   -varargin: Optionnal input arguments, used as inputs to
           %              plot()
           %
           %
           
           

            switch newfig
                case 0
                    hold on
                case 1
                    figure('Color','w','Position',[200 200 900 700])
                    hold on
                case 2
                    hold off
                otherwise
                    error('InputError: invalid newfig')
            end

            s = size(rhomn);
            N = (s(1)-1) / 2.0;

            if any(length(bn)~=N+1)
                error('InputError: bn has not the size ntor+1')
            end
            if any(length(r0c)~=N+1)
                error('InputError: r0c has not the size ntor+1')
            end
            if any(length(z0s)~=N+1)
                error('InputError: z0s has not the size ntor+1')
            end
            if nt<1
                error('InputError: nt should be larger than 1')
            end
            
            tarr = linspace(0, 2*pi, nt);
            rho = zeros(1, nt);
            alpha = obj.physicslist.twoalpha / 2.0;
            nfp = double(obj.physicslist.nfp);

            for im=1:s(2)
               mm = im-1;
               for in=1:s(1)
                  nn = in-N-1;
                  rho = rho + rhomn(in,im) * cos(mm*tarr + nn*nfp*phi - alpha*nfp*phi);
               end
            end
            
            r0 = 0;
            z0 = 0;
            b  = 0;
            for nn=0:N
                r0 = r0 + r0c(nn+1) * cos(nn*nfp*phi);
                b  = b  + bn(nn+1 ) * cos(nn*nfp*phi);
                z0 = z0 + z0s(nn+1) * sin(nn*nfp*phi);
            end
            
            zeta = b.*sin( tarr - alpha*nfp*phi );
            
            r = r0 + rho * cos(alpha*nfp*phi) - zeta * sin(alpha*nfp*phi);
            z = z0 + rho * sin(alpha*nfp*phi) + zeta * cos(alpha*nfp*phi);
            
            plot(r, z, varargin{:})
            axis equal
        end

        function plot_normal_field( obj, nt, phi, VorB, newfig )
            %
            % PLOT_NORMAL_FIELD( NT, PHI, VORB, NEWFIG )
            % ====================================
            %
            % Plots the normal field on the computational boundary as a
            % vector field
            %
            % INPUTS
            % ------
            %   -NT: Number of poloidal points
            %   -PHI: Toroidal angle
            %   -VORB: ='V': Only vacuum field
            %          ='B': Only plasma field
            %          ='F': Add vacuum to plasma field
            %   -NEWFIG: =0: plot on gca
            %            =1: plot on a new figure
            %            =2: erase and plot on gca

            % First, build coordinate
            tarr = linspace( 0, 2*pi, nt );
            r    = zeros( 1, nt );
            z    = zeros( 1, nt );
            
            rmnc = obj.physicslist.rwc;
            rmns = obj.physicslist.rws;
            zmnc = obj.physicslist.zwc;
            zmns = obj.physicslist.zws;

            Nfp = double(obj.physicslist.nfp);

            s = size(rmnc);
            
            if any(s~=size(rmns))
                error('InputError: rmnc has not the same size as rmns')
            end
            if any(s~=size(zmns))
                error('InputError: rmnc has not the same size as zmns')
            end
            if any(s~=size(zmnc))
                error('InputError: rmnc has not the same size as zmnc')
            end    
            if nt<1
                error('InputError: nt should be larger than 1')
            end
            
            N = (s(1)-1) / 2.0;
            nfp  = double(obj.physicslist.nfp);
            
            
            for in=1:s(1)
                nn = in-1-N;
                for im=1:s(2)
                    mm = im-1;
                    
                    arg = mm*tarr - nn*nfp*phi;
                    
                    r = r + rmnc(in,im) * cos(arg) + rmns(in,im) * sin(arg);
                    z = z + zmnc(in,im) * cos(arg) + zmns(in,im) * sin(arg);
                end
            end

            % Then evaluate norm of normal field
            switch VorB
                case 'V'
                    fmnc = obj.physicslist.vnc;
                    fmns = obj.physicslist.vns;
                case 'B'
                    fmnc = obj.physicslist.bnc;
                    fmns = obj.physicslist.bns;
                case 'F'
                    fmnc = obj.physicslist.bnc + obj.physicslist.vnc;
                    fmns = obj.physicslist.bns + obj.physicslist.vns;
                otherwise
                    error('InputError: Invalid VorB')
            end

            bnorm = zeros(1, nt);
            for mm=0:obj.mpol
                for nn=-obj.ntor:obj.ntor
                    if mm==0 && nn<0
                        continue
                    end
                    im = mm+1;
                    in = nn+obj.ntor+1;
                    arg = mm*tarr - nn*Nfp*phi;
                    
                    bnorm = bnorm + fmnc(in,im) * cos( arg ) ...
                                  + fmns(in,im) * sin( arg );
                end
            end
            
            bnorm = bnorm / max(bnorm);

            % Evaluate normal direction
            bnormal = zeros(2, nt);
            dt  = 1e-5;
            for it = 1:nt
                
                rj = 0; zj = 0; tj = tarr(it)-dt;
                rk = 0; zk = 0; tk = tarr(it)+dt;
                for in=1:s(1)
                    nn = in-1-N;
                    for im=1:s(2)
                        mm = im-1;

                        argj = mm*tj - nn*nfp*phi;
                        rj = rj + rmnc(in,im) * cos(argj) + rmns(in,im) * sin(argj);
                        zj = zj + zmnc(in,im) * cos(argj) + zmns(in,im) * sin(argj);

                        argk = mm*tk - nn*nfp*phi;
                        rk = rk + rmnc(in,im) * cos(argk) + rmns(in,im) * sin(argk);
                        zk = zk + zmnc(in,im) * cos(argk) + zmns(in,im) * sin(argk);
                    end
                end

                bnormal(1,it) =   zk-zj ;
                bnormal(2,it) = -(rk-rj);

                bnormal(:,it) = bnorm(it) * bnormal(:,it) / sqrt(bnormal(1,it)^2+bnormal(2,it)^2);
            end

            
            
            % Plot
            quiver( r, z, bnormal(1,:), bnormal(2,:) )
            
            

        end
    end
    
end


















