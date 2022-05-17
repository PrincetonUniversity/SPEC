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
        Mpol = 0;
        Ntor = 0;
        array_size = [0, 0];
        Mvol = 0;
        Nvol = 0;
        verbose = true;
        Lboundary = 0;
        
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
            %   -VARARGIN: Any couple of input:
            %       - 'Liniguess': set to true to read initial guess, false
            %                      to skip it. default: true
            %                       
            %
            % OUTPUTS
            % -------
            %   -OBJ: An instance of the class SPEC_Namelist
            %
            %
        
            % Read optional input
            l = length(varargin);
            if mod(l,2)~=0
                error('InputError: invalid number of inputs')
            end
            
            opt.Liniguess = true; % Decide whether or not we read initial guess
            opt.verbose = true;   % Print additional warnings
            for ii=1:l/2
               opt.(varargin{2*ii-1}) = varargin{2*ii}; 
            end
            
            % Fill class arguments
            obj.verbose = opt.verbose;
            
            % Read input file
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
            
            % Read initial guess
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
            % Read the initial guess from the input file filename. 
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
            save_line = false;            % This is switched to true once 
                                          % we are at the end of the
                                          % screelist
            category  = '';  % Save the name of the category

            initial_guess_str = {};

            while( ~feof(fid) ) % while it is not the end of the file

                tline = fgetl(fid);
                tline = strtrim(tline);
                
                if isempty(tline)
                    continue
                end

                if( save_line ) % write line
                    initial_guess_str{end+1} = tline;

                else % otherwise look for beginning of initial guess
                    % Check if beginning of new category
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
            l = l-2; % Remove m, n
            if mod(l,4)~=0
                error('Invalid number of modes in initial guess')
            end
            
            if l/4<obj.Nvol
                error('Not enough volumes in initial guess')
            end
            
            if l/4>obj.Nvol
                error('Too many volumes in initial guess')
            end
            
            % Read format of initial guess
            tmp = SPEC_Namelist( filename, 'Liniguess', false );
            
            if tmp.physicslist.Lboundary~=obj.physicslist.Lboundary
                error(['The initial guess from file %s does not use the',...
                       'same boundary representation.'], filename)
            end
            
            % Check if there is an initial guess...
            nlines = length(initial_guess_str);
            if nlines<1
                obj.initial_guess = struct([]); % generate empty structure
                return
            end

            % Read Mpol, Ntor
            Mpol_in = 0; Ntor_in = 0;
            for iline=1:nlines
                % Scan line
                line_data = str2num( initial_guess_str{iline} );

                % Find corresponding index
                mm = line_data(1); nn = line_data(2);
                Mpol_in = max([mm, Mpol_in]);
                Ntor_in = max([abs(nn), Ntor_in]);
            end
            
            % Check if resolution is smaller or larger than inner
            % resolution. This changes obj.Mpol and obj.Ntor if
            % necessary
            if (Mpol_in>obj.Mpol) || (Ntor_in>obj.Ntor)
               obj = obj.change_fourier_resolution( Mpol_in, Ntor_in );
            end
            
            % Now format initial guess in a structure
            switch obj.Lboundary
                case 0 % Rmn, Zmn representation
                    
                    % Allocate memory
                    Ric = zeros(2*obj.Ntor+1, obj.Mpol+1, obj.Mvol);
                    Ris = zeros(2*obj.Ntor+1, obj.Mpol+1, obj.Mvol);
                    Zic = zeros(2*obj.Ntor+1, obj.Mpol+1, obj.Mvol);
                    Zis = zeros(2*obj.Ntor+1, obj.Mpol+1, obj.Mvol);
                    
                    % Fill initial guess arrays
                    for iline=1:nlines

                        % Scan line
                        line_data = str2num( initial_guess_str{iline} );

                        % Find corresponding index
                        m = line_data(1); n = line_data(2);
                        im = m+1;
                        in = n+obj.Ntor+1;

                        for ivol=1:obj.Nvol
                            Ric(in,im,ivol) = line_data(ivol*4-1);
                            Zis(in,im,ivol) = line_data(ivol*4  );
                            Ris(in,im,ivol) = line_data(ivol*4+1);
                            Zic(in,im,ivol) = line_data(ivol*4+2);
                        end
                    end
                    
                    % Fill structure
                    obj.initial_guess.Ric = Ric;
                    obj.initial_guess.Zis = Zis;
                    obj.initial_guess.Ris = Ris;
                    obj.initial_guess.Zic = Zic;
                    
                case 1 % Henneberg representation
                                  
                    error('Henneberg representation not yet implemented')
                    
                otherwise
                    error('Invalid Lboundary!')
            end
            
        end
        
        % =================================================================
        % Getters
        function out = get_fourier_harmonics( obj, field, m, n, varargin )
            %
            % GET_FOURIER_HARMONICS( FIELD, M, N, (IVOL) )
            % ============================================
            %
            % Returns the fourier harmonics associated to the field,
            % poloidal mode number m and toroidal mode number n
            %
            % INPUTS
            % ------
            %   -field: 'Rbc', 'Zws', 'Ric', ...
            %   -m    : Poloidal mode number
            %   -n    : Toroidal mode number
            %   -(ivol): Optional argument, required for fourier harmonics
            %            of the initial guess
            
            
            % Check if poloidal mode is within the resolution
            if m>obj.Mpol
                warning('M is larger than the Fourier resolution')
                out = 0;
                return
            end 
            
            % Check if toroidal mode is within the resolution
            if abs(n)>obj.Ntor
                warning('N is larger than the Fourier resolution')
                out = 0;
                return
            end
           
            % Need to make the difference between initial guess harmonics
            % or harmonics from the physics list. rhoi, bin, R0ic and Z0is
            % are related to the Henneberg representation
            if     strcmp(field, 'Ric') || strcmp(field, 'Ris') ...
                || strcmp(field, 'Zic') || strcmp(field, 'Zis') ...
                || strcmp(field, 'rhoi') || strcmp(field, 'bin') ...
                || strcmp(field, 'R0ic') || strcmp(field, 'Z0is')
                
                if isempty(obj.initial_guess)
                    error('No initial guess available')
                end
                
                ivol = varargin{1};
            
                if ivol<1 || ivol>obj.Nvol
                    error('Invalid ivol')
                end                   
                
                cat = 'initial_guess';
                
            else
                cat = 'physicslist';
                ivol = 1;
                
            end
                
            % Read relevant harmonic.
            out = obj.(cat).(field)(n+obj.Ntor+1, m+1, ivol);
            
        end
        
        % =================================================================
        % Setter
        function obj = set_fourier_harmonics( obj, field, im, in, value, varargin )
            %
            % SET_FOURIER_HARMONICS( FIELD, M, N, (IVOL) )
            % ============================================
            %
            % Sets the fourier harmonics associated to the field,
            % poloidal mode number m and toroidal mode number n
            %
            % INPUTS
            % ------
            %   -field: 'Rbc', 'Zws', 'Ric', ...
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
            
            % Loop over the input harmonics
            for imn=1:mn
                m = im(imn);
                n = in(imn);
                
                % If some modes are greater than the actual resolution,
                % increase the resolution accordingly
                if m>obj.Mpol
                    warning('Poloidal resolution has to be increased ...')
                    obj = obj.change_fourier_resolution( m, obj.Ntor );
                end 

                if abs(n)>obj.Ntor
                    warning('Toroidal resolution has to be increased ...')
                    obj = obj.change_fourier_resolution( obj.Mpol, abs(n) );
                end

                % Check if the quantity is located in the initial guess
                % structure or in the physicslist structure.
                if     strcmp(field, 'Ric') || strcmp(field, 'Ris') ...
                    || strcmp(field, 'Zic') || strcmp(field, 'Zis') ...
                    || strcmp(field, 'rhoi') || strcmp(field, 'bin') ...
                    || strcmp(field, 'R0ic') || strcmp(field, 'Z0is')

                    % If the initial guess is empty (i.e. no initial guess
                    % are available), then we create one filled with zeros
                    if isempty(obj.initial_guess)
                        warning('No initial guess available... filling with zeros')

                        if obj.Lboundary == 0
                            obj.initial_guess.Ric = zeros( obj.array_size(1), obj.array_size(2), obj.Nvol );
                            obj.initial_guess.Ris = zeros( obj.array_size(1), obj.array_size(2), obj.Nvol );
                            obj.initial_guess.Zic = zeros( obj.array_size(1), obj.array_size(2), obj.Nvol );
                            obj.initial_guess.Zis = zeros( obj.array_size(1), obj.array_size(2), obj.Nvol );                    
                        else
                            error('Henneberg representation not yet implemented' )
                        end
                    end

                    ivol = varargin{1};

                    if ivol<1 || ivol>obj.Nvol
                        error('Invalid ivol')
                    end                   

                    cat = 'initial_guess';

                else
                    cat = 'physicslist';
                    ivol = 1;

                end

                % Set the relevant harmonics
                obj.(cat).(field)(n+obj.Ntor+1, m+1, ivol) = value(imn);
            end
            
        end
        
        function obj = truncate_fourier_series( obj, Mpol, Ntor )
            %
            % TRUNCATE_FOURIER_SERIES( MPOL, NTOR )
            % =====================================
            %
            % Truncates all spectral quantities to the requested poloidal
            % and toroidal resolution. This can also be used to increase
            % the Fourier resolution.
            %
            % INPUTS
            % ------
            %   -Mpol: Poloidal resolution
            %   -Ntor: Toroidal resolution
            %
            % OUTPUT
            % ------
            %   -OBJ: Truncated instance of SPEC_Namelist
            %
           
            obj = obj.change_fourier_resolution( Mpol, Ntor );
            
        end
        
        function obj = change_boundary_representation( obj, new_lboundary )
            % 
            % CHANGE_BOUNDARY_REPRESENTATION( NEW_LBOUNDARY )
            % ===============================================
            %
            % Changes from the hudson representation to the henneberg's one
            % and vice versa.
            %
            % INPUT
            % -----
            %   -new_lboundary: New value for Lboundary
            
            error('Henneberg representation not yet implemented')
            
        end
        
        % =================================================================
        % Plotters
        function plot_plasma_boundary( obj, nt, phi, newfig, varargin )
           %
           % PLOT_PLASMA_BOUNDARY
           % ====================
           %
           % Plot the boundary given by Rbc, Zbs, Rbs, Zbc
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
           
           if obj.Lboundary==0
               Rbc = obj.physicslist.Rbc;
               Rbs = obj.physicslist.Rbs;
               Zbc = obj.physicslist.Zbc;
               Zbs = obj.physicslist.Zbs;

               obj.plot_surface_Lb0( Rbc, Zbs, Rbs, Zbc, nt, phi, newfig, varargin{:} )
           else
               
               error('Henneberg representation is not implemented')
           end
        end
        
        function plot_computational_boundary( obj, nt, phi, newfig, varargin )
           %
           % PLOT_PLASMA_BOUNDARY
           % ====================
           %
           % Plot the boundary given by Rbc, Zbs, Rbs, Zbc
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
           
           if obj.Lboundary == 0
               Rwc = obj.physicslist.Rwc;
               Rws = obj.physicslist.Rws;
               Zwc = obj.physicslist.Zwc;
               Zws = obj.physicslist.Zws;

               obj.plot_surface_Lb0( Rwc, Zws, Rws, Zwc, nt, phi, newfig, varargin{:} )
           else
               error('Henneberg representation not yet implemented')
           end
        end
           
        function plot_initial_guess( obj, nt, phi, newfig, varargin )
           %
           % PLOT_PLASMA_BOUNDARY( NT, PHI, NEWFIG, VARARGIN )
           % =================================================
           %
           % Plot the boundary given by Rbc, Zbs, Rbs, Zbc
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
           
           if obj.Lboundary == 0
               for ivol=1:obj.Nvol
                   Ric = obj.initial_guess.Ric(:,:,ivol);
                   Ris = obj.initial_guess.Ris(:,:,ivol);
                   Zic = obj.initial_guess.Zic(:,:,ivol);
                   Zis = obj.initial_guess.Zis(:,:,ivol);

                   if ivol>=2
                       newfig=0;
                   end

                   obj.plot_surface_Lb0( Ric, Zis, Ris, Zic, nt, phi, newfig, varargin{:} )
               end
           else
              
               error('Henneberg representation not yet implemented' )
           end
        end
        
        % =================================================================
        % Write method
        function write_input_file(obj, filename )
            %
            % WRITE_INPUT_FILE( filename )
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
            if obj.Ntor==0
                obj = obj.change_fourier_resolution( obj.Mpol, 1 );
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
            % elements with negative indices.
            S.shift.Rbc = [obj.Ntor+1, 1];
            S.shift.Rbs = [obj.Ntor+1, 1];
            S.shift.Zbc = [obj.Ntor+1, 1];
            S.shift.Zbs = [obj.Ntor+1, 1];
            S.shift.Rwc = [obj.Ntor+1, 1];
            S.shift.Rws = [obj.Ntor+1, 1];
            S.shift.Zws = [obj.Ntor+1, 1];
            S.shift.Zwc = [obj.Ntor+1, 1];
            S.shift.Vnc = [obj.Ntor+1, 1];
            S.shift.Vns = [obj.Ntor+1, 1];
            S.shift.Bnc = [obj.Ntor+1, 1];
            S.shift.Bns = [obj.Ntor+1, 1];
            S.shift.rhomn = [obj.Ntor+1, 1];

            % Remove unnecessary fields.
            if obj.Lboundary == 0
                if isfield( S.physicslist, 'rhomn' )
                    S.physicslist = rmfield(S.physicslist, 'rhomn');
                end
                if isfield( S.physicslist, 'bn' )
                    S.physicslist = rmfield(S.physicslist, 'bn');
                end
                if isfield( S.physicslist, 'R0c' )
                    S.physicslist = rmfield(S.physicslist, 'R0c');
                end
                if isfield( S.physicslist, 'Z0s' )
                    S.physicslist = rmfield(S.physicslist, 'Z0s');
                end
                
            else % Lboundary==1
                error('Henneberg representation not yet implemented')            
                
            end
            
            if obj.physicslist.Lfreebound==0 % Then no need to freeboundary info
                if isfield( S.physicslist, 'Rwc' )
                    S.physicslist = rmfield(S.physicslist, 'Rwc');
                end
                if isfield( S.physicslist, 'Rws' )
                    S.physicslist = rmfield(S.physicslist, 'Rws');
                end
                if isfield( S.physicslist, 'Zwc' )
                    S.physicslist = rmfield(S.physicslist, 'Zwc');
                end
                if isfield( S.physicslist, 'Zws' )
                    S.physicslist = rmfield(S.physicslist, 'Zws');
                end 
                if isfield( S.physicslist, 'Vnc' )
                    S.physicslist = rmfield(S.physicslist, 'Vnc');
                end
                if isfield( S.physicslist, 'Vns' )
                    S.physicslist = rmfield(S.physicslist, 'Vns');
                end
                if isfield( S.physicslist, 'Bnc' )
                    S.physicslist = rmfield(S.physicslist, 'Bnc');
                end
                if isfield( S.physicslist, 'Bns' )
                    S.physicslist = rmfield(S.physicslist, 'Bns');
                end 
                
            end
            
            
            % Build initial guess strings
            if obj.physicslist.Lboundary==0        
                if ~isempty(obj.initial_guess)
                    s = size(obj.initial_guess.Ric);
                    initialguess = cell(1, s(1)*s(2));
                    iline=0;
                    for ii=1:s(1)
                        for jj=1:s(2)
                            iline = iline+1;
                            mm = jj-1;
                            nn = ii-obj.Ntor-1;
                            
                            initialguess{iline} = sprintf( '%i   %i   ', mm, nn );
                            for ivol=1:s(3)
                               initialguess{iline} = sprintf( '%s   %0.12E   %0.12E   %0.12E   %0.12E', initialguess{iline}, ...
                                                            obj.initial_guess.Ric(ii,jj,ivol), ...
                                                            obj.initial_guess.Zis(ii,jj,ivol), ...
                                                            obj.initial_guess.Ris(ii,jj,ivol), ...
                                                            obj.initial_guess.Zic(ii,jj,ivol)     );
                            end
                        end
                    end
                else
                    initialguess = cell(0);
                end

            else

                error('Henneberg representation not yet implemented')



            end

            % Call namelist writer
            write_namelist( S, filename, initialguess );
       

        end
    end
    
    methods (Access=private)
        
        function obj = initialize_structure( obj )
            %
            % INITIALIZE_STRUCTURE( OBJ )
            % ===========================
            %
            % Use to check that all required inputs are correctly set;
            % check that the sizes of arrays are correct. Raise errors in
            % case crucial informations is missing.
            %
            % 
                       
            % Fill some important inputs
            if ~isfield(obj.physicslist, 'Nvol')
                error('Missing Nvol')
            end
            if ~isfield(obj.physicslist, 'Lfreebound')
                error('Missing Lfreebound')
            end
            
            obj.Mvol = obj.physicslist.Nvol + obj.physicslist.Lfreebound;
            obj.Nvol = obj.physicslist.Nvol;
            
            if isfield(obj.physicslist, 'Lboundary')
                obj.Lboundary = obj.physicslist.Lboundary;
            else
                obj.Lboundary = 0;
            end
            
            % PHYSICSLIST
            % -----------
            % Lrad
            if ~isfield(obj.physicslist, 'Lrad') && obj.verbose
               warning('Missing Lrad. Filling with 4...')
               obj.physicslist.Lrad = ones(1,obj.Mvol) * 4;
            else
                % Fill potential missing elements with 4s
               obj.physicslist.Lrad(end+1:obj.Mvol) = 4;                
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
               warning('Missing Lrad. Filling with 0...')
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
            
            % Ladiabatic
            if ~isfield(obj.physicslist, 'Ladiabatic')
                warning('Missing Ladiabatic. Setting to zero...')
                obj.physicslist.Ladiabatic = 0.0;
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
            if obj.physicslist.Ladiabatic==1
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
            
            % Lconstraint
            if ~isfield(obj.physicslist, 'Lconstraint')
               warning('Missing Lconstraint. Setting to 0...')
               obj.physicslist.Lconstraint = 0;              
            end
            
            if ~any(obj.physicslist.Lconstraint==[0,1,2,3])
                error('Invalid Lconstraint')
            end
            
            
            % Ivolume, Isurf
            if obj.physicslist.Lconstraint==3
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
            
            if obj.physicslist.Lconstraint==1
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
            if ~isfield(obj.physicslist, 'Mpol')
                error('Missing Mpol information')
            end
            if ~isfield(obj.physicslist, 'Ntor')
                error('Missing Ntor information')
            end
            if ~isfield(obj.physicslist, 'Lboundary')
               warning('Missing Lboundary. Setting to zero')
               obj.physicslist.Lboundary = 0;
            end
            
            Mpol_in = obj.physicslist.Mpol;
            Ntor_in = obj.physicslist.Ntor;
            
            if obj.Lboundary==0
                if ~isfield(obj.physicslist, 'Rbc')
                    obj.physicslist.shift.Rbc = [Ntor_in+1, 1];
                    obj.physicslist.Rbc = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Rbs')
                    obj.physicslist.shift.Rbs = [Ntor_in+1, 1];
                    obj.physicslist.Rbs = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Zbc')
                    obj.physicslist.shift.Zbc = [Ntor_in+1, 1];
                    obj.physicslist.Zbc = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Zbs')
                    obj.physicslist.shift.Zbs = [Ntor_in+1, 1];
                    obj.physicslist.Zbs = zeros(2*Ntor_in+1, Mpol_in);
                end
                
                % Check that sizes are consistent with each others
                if any(size(obj.physicslist.Rbc)~=size(obj.physicslist.Rbs))
                    error('Size mismatch between Rbc and Rbs')
                end
                if any(size(obj.physicslist.Rbc)~=size(obj.physicslist.Zbc))
                    error('Size mismatch between Rbc and Zbc')
                end
                if any(size(obj.physicslist.Rbc)~=size(obj.physicslist.Zbs))
                    error('Size mismatch between Rbc and Zbs')
                end
                
                
            else %Lboundary==1
                
                error('Henneberg representation not yet implemented')
                
                
            end
            
             

            if ~isfield(obj.physicslist, 'Rwc')
                obj.physicslist.shift.Rwc = [Ntor_in+1, 1];
                obj.physicslist.Rwc = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Rws')
                obj.physicslist.shift.Rws = [Ntor_in+1, 1];
                obj.physicslist.Rws = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Zwc')
                obj.physicslist.shift.Zwc = [Ntor_in+1, 1];
                obj.physicslist.Zwc = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Zws')
                obj.physicslist.shift.Zws = [Ntor_in+1, 1];
                obj.physicslist.Zws = zeros(2*Ntor_in+1, Mpol_in);
            end

            % Check that sizes are consistent with each others
            if any(size(obj.physicslist.Rwc)~=size(obj.physicslist.Rws))
                error('Size mismatch between Rwc and Rws')
            end
            if any(size(obj.physicslist.Rwc)~=size(obj.physicslist.Zwc))
                error('Size mismatch between Rwc and Zwc')
            end
            if any(size(obj.physicslist.Rwc)~=size(obj.physicslist.Zws))
                error('Size mismatch between Rwc and Zws')
            end            
                
                
            
            if ~isfield(obj.physicslist, 'Vnc')
                obj.physicslist.shift.Vnc = [Ntor_in+1, 1];
                obj.physicslist.Vnc = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Vns')
                obj.physicslist.shift.Vns = [Ntor_in+1, 1];
                obj.physicslist.Vns = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Bnc')
                obj.physicslist.shift.Bnc = [Ntor_in+1, 1];
                obj.physicslist.Bnc = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Bns')
                obj.physicslist.shift.Bns = [Ntor_in+1, 1];
                obj.physicslist.Bns = zeros(2*Ntor_in+1, Mpol_in);
            end
            
            % Check that sizes are consistent with each others
            if any(size(obj.physicslist.Vnc)~=size(obj.physicslist.Vns))
                error('Size mismatch between Vnc and Vns')
            end
            if any(size(obj.physicslist.Bnc)~=size(obj.physicslist.Bns))
                error('Size mismatch between Bnc and Bns')
            end
            
            % DIAGNOSTICSLIST
            % ---------------
            if ~isfield(obj.diagnosticslist, 'nPpts')
                warning('Missing nPpts. Setting to zero')
                obj.diagnosticslist.nPpts = 0;
            end
            
            if ~isfield(obj.diagnosticslist, 'nPtrj')
                warning('Missing nPtrj. Setting to zero')
                obj.diagnosticslist.nPtrj = 0;
            end            
        end
        
        
        function obj = change_fourier_resolution( obj, Mpol_new, Ntor_new )
            %
            % CHANGE_FOURIER_RESOLUTION( MPOL_NEW, NTOR_NEW )
            % ===============================================
            %
            % Change inner Fourier resolution of an instance of
            % SPEC_Namelist.
            %
            % INPUTS
            % ------
            %   -Mpol_new: New poloidal resolution
            %   -Ntor_new: New toroidal resolution
            %
            % OUTPUT
            % ------
            %   -obj: Updated instance of SPEC_Namelist
            %
       
            if Mpol_new<1
                error('InputError, Mpol_new should be larger than 0')
            end
            if Ntor_new<0
                error('InputError, Ntor_new should be larger or equal to zero')
            end
            
            obj.Mpol = Mpol_new;
            obj.Ntor = Ntor_new;
            obj.array_size = [2*obj.Ntor+1, obj.Mpol+1];
            
            % Check that all arrays have the same size; otherwise, fill
            % with zeros the missing elements
            if obj.Lboundary == 0
                if any(size(obj.physicslist.Rbc)~=obj.array_size)
                    obj = obj.reshape_array( 'Rbc' );
                    obj = obj.reshape_array( 'Rbs' );
                    obj = obj.reshape_array( 'Zbc' );
                    obj = obj.reshape_array( 'Zbs' );
                end

                if any(size(obj.physicslist.Rwc)~=obj.array_size)
                    obj = obj.reshape_array( 'Rwc' );
                    obj = obj.reshape_array( 'Rws' );
                    obj = obj.reshape_array( 'Zwc' );
                    obj = obj.reshape_array( 'Zws' );
                end
            else
                error('Henneberg representation not yet implemented')
                
            end
                
            
            if any(size(obj.physicslist.Vnc)~=obj.array_size)
                obj = obj.reshape_array( 'Vnc' );
                obj = obj.reshape_array( 'Vns' );
            end
            
            if any(size(obj.physicslist.Bnc)~=obj.array_size)
                obj = obj.reshape_array( 'Bnc' );
                obj = obj.reshape_array( 'Bns' );           
            end
        end
        
        
        function obj = set_fourier_resolution( obj )
            %
            % SET_FOURIER_RESOLUTION( OBJ )
            % =============================
            %
            % Set the obj internal Fourier resolution to the maximal value
            % required to store the given data; then, reshapees all arrays
            % to have the size obj.array_size.
            %
            % Check beforehand that the resolution of Rbc is the same as
            % Rbs, Zbc and Zbs. Do something similar for Vnc, Vns, Bnc, Bns
            % and Rwc, Rws, Zws, Zwc.
            %
           
            
            % Find largest Fourier resolution in the input file
            Mpol_in = obj.physicslist.Mpol;
            Ntor_in = obj.physicslist.Ntor;
            
            if obj.Lboundary == 0
                s_bc = size(obj.physicslist.Rbc);
                shift = obj.physicslist.shift.Rbc(1);
                Mpol_bc = s_bc(2)-1       ;
                Ntor_bc = max([abs(1-shift), s_bc(1)-shift]);
            else
                error('Henneberg representation not yet implemented')
            end
            
            s_wc = size(obj.physicslist.Rwc);
            shift = obj.physicslist.shift.Rwc(1);
            Mpol_wc = s_wc(2)-1       ;
            Ntor_wc = max([abs(1-shift), s_wc(1)-shift]);
            
            s_vb = size(obj.physicslist.Vnc);
            shift = obj.physicslist.shift.Vnc(1);
            Mpol_vb = s_vb(2)-1       ;
            Ntor_vb = max([abs(1-shift), s_vb(1)-shift]);
            
            Mpol_new = max([Mpol_in, Mpol_bc, Mpol_wc, Mpol_vb]);
            Ntor_new = max([Ntor_in, Ntor_bc, Ntor_wc, Ntor_vb]);
            obj = obj.change_fourier_resolution( Mpol_new, Ntor_new );
        end
        
        
        function obj = reshape_array( obj, field )
            %
            % RESHAPE_ARRAY( OBJ, FIELD )
            % =======================================
            %
            % Reshape the input array in a an array of size obj.array_size
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
            Ntor_array = max([abs(1-shift(1)), s(1)-shift(1)]);   
            
            new_array = zeros(obj.array_size);
            for ii=1:s(1)
                for jj=1:s(2)
                    nn = ii-Ntor_array-1;
                    mm = jj-shift(2);
                    
                    if abs(nn)>obj.Ntor || mm>obj.Mpol
                        continue
                    end
                    
                    new_array(nn+obj.Ntor+1, mm+1) = array(ii, jj);
                end
            end
            
            obj.physicslist.(field) = new_array;
        end
        
        
        function plot_surface_Lb0( obj, Rmnc, Zmns, Rmns, Zmnc, nt, phi, ...
                               newfig, varargin )
            %
            % PLOT_SURFACE( RMNC, ZMNS, RMNS, ZMNC, NT, NEWFIG, VARARGIN )
            % ===========================================================
            %
            % Plot a surface parametrized by the standard representation
            %
            % INPUTS
            % ------
            %   -Rmnc: Even Fourier modes of R, format (2*Ntor+1, Mpol+1)
            %   -Zmns: Odd  Fourier modes of Z, format (2*Ntor+1, Mpol+1)
            %   -Rmns: Odd  Fourier modes of R, format (2*Ntor+1, Mpol+1)
            %   -Zmnc: Even Fourier modes of Z, format (2*Ntor+1, Mpol+1)
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
            
            s = size(Rmnc);
            
            if any(s~=size(Rmns))
                error('InputError: Rmnc has not the same size as Rmns')
            end
            if any(s~=size(Zmns))
                error('InputError: Rmnc has not the same size as Zmns')
            end
            if any(s~=size(Zmnc))
                error('InputError: Rmnc has not the same size as Zmnc')
            end    
            if nt<1
                error('InputError: nt should be larger than 1')
            end
            
            N = (s(1)-1) / 2.0;
            
            tarr = linspace( 0, 2*pi, nt );
            R    = zeros( 1, nt );
            Z    = zeros( 1, nt );
            Nfp  = double(obj.physicslist.Nfp);
            
            
            for in=1:s(1)
                nn = in-1-N;
                for im=1:s(2)
                    mm = im-1;
                    
                    arg = mm*tarr - nn*Nfp*phi;
                    
                    R = R + Rmnc(in,im) * cos(arg) + Rmns(in,im) * sin(arg);
                    Z = Z + Zmnc(in,im) * cos(arg) + Zmns(in,im) * sin(arg);
                end
            end
            
            plot( R, Z, varargin{:} )
            
            axis equal
        end
        
        
        function plot_surface_Lb1( obj, rhomn, bn, R0c, Z0s, nt, phi, newfig, varargin )
           %
           % PLOT_SURFACE_LB1( RHOMN, BN, R0C, Z0S, NT, PHI, NEWFIG, VARARGIN )
           % ==================================================================
           %
           % Plots a surface using the Henneberg representation
           %
           % INPUTS
           % ------
           %   -rhomn:  rho_mn harmonics, format (2*Ntor+1, Mpol+1)
           %   -bn:     b_n harmonics, format (Ntor+1, 1)
           %   -R0c:    R_0c harmonics, format (Ntor+1, 1)
           %   -Z0s:    Z_0s harmonics, format (Ntor+1, 1)
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
                error('InputError: bn has not the size Ntor+1')
            end
            if any(length(R0c)~=N+1)
                error('InputError: R0c has not the size Ntor+1')
            end
            if any(length(Z0s)~=N+1)
                error('InputError: Z0s has not the size Ntor+1')
            end
            if nt<1
                error('InputError: nt should be larger than 1')
            end
            
            tarr = linspace(0, 2*pi, nt);
            rho = zeros(1, nt);
            alpha = obj.physicslist.twoalpha / 2.0;
            Nfp = double(obj.physicslist.Nfp);

            for im=1:s(2)
               mm = im-1;
               for in=1:s(1)
                  nn = in-N-1;
                  rho = rho + rhomn(in,im) * cos(mm*tarr + nn*Nfp*phi - alpha*Nfp*phi);
               end
            end
            
            R0 = 0;
            Z0 = 0;
            b  = 0;
            for nn=0:N
                R0 = R0 + R0c(nn+1) * cos(nn*Nfp*phi);
                b  = b  + bn(nn+1 ) * cos(nn*Nfp*phi);
                Z0 = Z0 + Z0s(nn+1) * sin(nn*Nfp*phi);
            end
            
            zeta = b.*sin( tarr - alpha*Nfp*phi );
            
            R = R0 + rho * cos(alpha*Nfp*phi) - zeta * sin(alpha*Nfp*phi);
            Z = Z0 + rho * sin(alpha*Nfp*phi) + zeta * cos(alpha*Nfp*phi);
            
            plot(R, Z, varargin{:})
            axis equal
        end
    end
    
end


















