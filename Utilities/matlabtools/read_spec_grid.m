function gdata = read_spec_grid(filename)

% Reads coordinate grid data (grid points, jacobian, cartesian components of B, etc.) 
%
% INPUT
% - filename : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%
% OUTPUT
% - gdata    : contains all the coordinate grid data, which can be fed into several routines for analyzing and ploting
%
% Note: So far only works if Lrad is the same for all volumes
%
% written by J.Loizu (2017)


global machform;

machform = 's';

data = read_hdf5(filename);

nvol = double(data.Mvol);
data.Rij = cell(nvol,1);
data.Zij = cell(nvol,1);
data.sg  = cell(nvol,1);
data.BR  = cell(nvol,1);
data.Bp  = cell(nvol,1);
data.BZ  = cell(nvol,1);

% Read the grid files

machine_format = machform;   % needs to be 's' (intel) or 'a' (gnu); format is swapped if an error occurs 
triedallform   =  0;         % whether all formats have been tried (1) or not (0) 
success        =  0;
int_format     = 'int32';
float_format   = 'float64';
spacer_format  = 'int32';

while(success==0)
try

  if(machine_format ~= machform)
   machine_format =  machform;  % update value
  end
  [filepath,name,ext]=fileparts(filename);
  grid_file = [filepath filesep '.' name '.grid'];
  fid        = fopen(grid_file,'r',machine_format);
  if (fid > 0)
    % Get Nt, Nz, Ntz, Mvol, Igeometry, pi2nfp
    fread(fid,1,spacer_format);
    param    = fread(fid,5,int_format);
    fread(fid,1,float_format);
    fread(fid,1,spacer_format);
    Nt       = param(1);
    Nz       = param(2);
    Ntz      = param(3);
    data.Nt  = Nt;
    data.Nz  = Nz;
    data.Ntz = Ntz; 
  else
   disp([' - File "' grid_file '" does not exist']);
  end

  for i=1:nvol
   % Get radial resolution
   fread(fid,1,spacer_format);
   Lrad         = fread(fid,1,int_format);
   fread(fid,1,spacer_format);
   data.Lrad(i) = Lrad;
   data.Rij{i} = zeros(Ntz,Lrad+1);
   data.Zij{i} = zeros(Ntz,Lrad+1);
   data.sg{i} = zeros(Ntz,Lrad+1);
   data.BR{i} = zeros(Ntz,Lrad+1);
   data.Bp{i} = zeros(Ntz,Lrad+1);
   data.BZ{i} = zeros(Ntz,Lrad+1);
   for l=1:Lrad+1
    % Get grid points
    fread(fid,1,spacer_format);
    Rout = fread(fid,Ntz,float_format);
    fread(fid,1,spacer_format);
    fread(fid,1,spacer_format);
    Zout = fread(fid,Ntz,float_format);
    fread(fid,1,spacer_format);   
    data.Rij{i}(:,l) = Rout;
    data.Zij{i}(:,l) = Zout;
    
    fread(fid,1,spacer_format);
    sg = fread(fid,Ntz,float_format);
    fread(fid,1,spacer_format);
    fread(fid,1,spacer_format);
    ijreal = fread(fid,Ntz,float_format);
    fread(fid,1,spacer_format);
    fread(fid,1,spacer_format);
    ijimag = fread(fid,Ntz,float_format);
    fread(fid,1,spacer_format);
    fread(fid,1,spacer_format);
    jireal = fread(fid,Ntz,float_format);
    fread(fid,1,spacer_format);
    data.sg{i}(:,l) = sg;
    data.BR{i}(:,l) = ijreal;
    data.Bp{i}(:,l) = ijimag;
    data.BZ{i}(:,l) = jireal;
   end
  end
  
  success  = 1;
  
catch

  if(triedallform==1)      
   disp(' - Could not read poincare file'); break;
  end
  if(machform~='s')
   machform     ='s'; 
   triedallform = 1;
  else
   machform     ='a';
   triedallform = 1;
  end
  if (fid ~= -1)
   fclose(fid);
   fid = -1;
  end
end

end

if (fid ~= -1) 
    fclose(fid);
end
%return output
gdata = data;

