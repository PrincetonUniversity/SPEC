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
        
    end
    
    methods (Access=public)
        % Class constructor
        function obj = SPEC_Namelist( filename )
        
            work = read_namelist( filename );
            obj.lists = fields(work);
            obj.physicslist = work.physicslist;
            obj.numericlist = work.numericlist;
            obj.locallist   = work.locallist;
            obj.globallist  = work.globallist;
            obj.diagnosticslist = work.diagnosticslist;
            obj.screenlist  = work.screenlist;  
            obj.initial_guess = [];
            
            % Check that the size of arrays makes sense, fills with zeros
            % otherwise
            obj = obj.initialize_structure();
            
            % Find the largest Fourier resolution used in the input file;
            % reformat all spectral quantities to have the same resolution
            obj = obj.set_fourier_resolution();
            
            
            
        end
        
    end
    
    methods (Access=private)
        
        function obj = initialize_structure( obj )
                       
            % Fill some important inputs
            if ~isfield(obj.physicslist, 'Nvol')
                error('Missing Nvol')
            end
            if ~isfield(obj.physicslist, 'Lfreebound')
                error('Missing Lfreebound')
            end
            
            obj.Mvol = obj.physicslist.Nvol + obj.physicslist.Lfreebound;
            
            % PHYSICSLIST
            % -----------
            % Lrad
            if ~isfield(obj.physicslist, 'Lrad')
               warning('Missing Lrad. Filling with 4...')
               obj.physicslist.Lrad = ones(1,obj.Mvol) * 4;
            else
                % Fill potential missing elements with 4s
               obj.physicslist.Lrad(end+1:obj.Mvol) = 4;                
            end
            
            % tflux
            if ~isfield(obj.physicslist, 'tflux')
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
            
            if obj.physicslist.Lboundary==0
                if ~isfield(obj.physicslist, 'Rbc')
                    obj.physicslist.Rbc = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Rbs')
                    obj.physicslist.Rbs = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Zbc')
                    obj.physicslist.Zbc = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Zbs')
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
                
                

                if ~isfield(obj.physicslist, 'Rwc')
                    obj.physicslist.Rwc = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Rws')
                    obj.physicslist.Rws = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Zwc')
                    obj.physicslist.Zwc = zeros(2*Ntor_in+1, Mpol_in);
                end
                if ~isfield(obj.physicslist, 'Zws')
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
                
                
                
            else
                error('To complete')
            end
            
            if ~isfield(obj.physicslist, 'Vnc')
                obj.physicslist.Vnc = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Vns')
                obj.physicslist.Vns = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Bnc')
                obj.physicslist.Bnc = zeros(2*Ntor_in+1, Mpol_in);
            end
            if ~isfield(obj.physicslist, 'Bns')
                obj.physicslist.Bns = zeros(2*Ntor_in+1, Mpol_in);
            end
            
            % Check that sizes are consistent with each others
            if any(size(obj.physicslist.Vnc)~=size(obj.physicslist.Vns))
                error('Size mismatch between Rwc and Vns')
            end
            if any(size(obj.physicslist.Vnc)~=size(obj.physicslist.Bnc))
                error('Size mismatch between Rwc and Bnc')
            end
            if any(size(obj.physicslist.Vnc)~=size(obj.physicslist.Bns))
                error('Size mismatch between Rwc and Bns')
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
       
        
        
        
        
        function obj = set_fourier_resolution( obj )
            
            % Check beforehand that the resolution of Rbc is the same as
            % Rbs, Zbc and Zbs. Do something similar for Vnc, Vns, Bnc, Bns
            % and Rwc, Rws, Zws, Zwc.
           
            % Find largest Fourier resolution in the input file
            Mpol_in = obj.physicslist.Mpol;
            Ntor_in = obj.physicslist.Ntor;
            
            s_bc = size(obj.physicslist.Rbc);
            shift = obj.physicslist.shift.Rbc(1);
            Mpol_bc = s_bc(2)-1       ;
            Ntor_bc = max([abs(1-shift), s_bc(1)-shift]);
            
            s_wc = size(obj.physicslist.Rwc);
            shift = obj.physicslist.shift.Rwc(1);
            Mpol_wc = s_wc(2)-1       ;
            Ntor_wc = max([abs(1-shift), s_wc(1)-shift]);
            
            s_vb = size(obj.physicslist.Vnc);
            shift = obj.physicslist.shift.Vnc(1);
            Mpol_vb = s_vb(2)-1       ;
            Ntor_vb = max([abs(1-shift), s_vb(1)-shift]);
            
            obj.Mpol = max([Mpol_in, Mpol_bc, Mpol_wc, Mpol_vb]);
            obj.Ntor = max([Ntor_in, Ntor_bc, Ntor_wc, Ntor_vb]);
            obj.array_size = [2*obj.Ntor+1, obj.Mpol+1];
            
            % Check that all arrays have the same size; otherwise, fill
            % with zeros the missing elements
            if any(size(obj.physicslist.Rbc)~=obj.array_size)
                obj.physicslist.Rbc = obj.reshape_array( obj.physicslist.Rbc, Ntor_bc );
                obj.physicslist.Rbs = obj.reshape_array( obj.physicslist.Rbs, Ntor_bc );
                obj.physicslist.Zbc = obj.reshape_array( obj.physicslist.Zbc, Ntor_bc );
                obj.physicslist.Zbs = obj.reshape_array( obj.physicslist.Zbs, Ntor_bc );
            end
            
            if any(size(obj.physicslist.Rwc)~=obj.array_size)
                obj.physicslist.Rwc = obj.reshape_array( obj.physicslist.Rwc, Ntor_wc );
                obj.physicslist.Rws = obj.reshape_array( obj.physicslist.Rws, Ntor_wc );
                obj.physicslist.Zwc = obj.reshape_array( obj.physicslist.Zwc, Ntor_wc );
                obj.physicslist.Zws = obj.reshape_array( obj.physicslist.Zws, Ntor_wc );
            end
            
            if any(size(obj.physicslist.Vnc)~=obj.array_size)
                obj.physicslist.Vnc = obj.reshape_array( obj.physicslist.Vnc, Ntor_vb );
                obj.physicslist.Vns = obj.reshape_array( obj.physicslist.Vns, Ntor_vb );
                obj.physicslist.Bnc = obj.reshape_array( obj.physicslist.Bnc, Ntor_vb );
                obj.physicslist.Bns = obj.reshape_array( obj.physicslist.Bns, Ntor_vb );           
            end
        end
        
        
        function new_array = reshape_array( obj, array, Ntor_array )
        
            s = size(array);
            
            new_array = zeros(obj.array_size);
            for ii=1:s(1)
                for jj=1:s(2)
                    nn = ii-Ntor_array-1;
                    mm = jj-1;
                    
                    new_array(nn+obj.Ntor+1, mm+1) = array(ii, jj);
                end
            end
        end
    end
    
end