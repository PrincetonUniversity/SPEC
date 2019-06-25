function fdata = read_spec_field(filename)

% Reads coefficients of the vector potential using output from SPEC 
%
% INPUT
% - filename : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%
% OUTPUT
% - fdata    : contains all the field data, which can be fed into several routines for analyzing and ploting
%
% Note: So far only works if Lrad is the same for all volumes
%
% written by J.Loizu (2017)
% modified by J.Loizu (10.2017)
% modified by J.Loizu (02.2018)

global machform;

machform  = 's';

data      = read_hdf5(filename);

nvol      = data.Mvol;

mn        = data.mn;

data.Lrad = zeros(nvol,1);  % allocate data set for Lrad(1:nvol)

data.Ate = cell(nvol,1);    % create cells for each volume's field 

data.Aze = cell(nvol,1);

data.Ato = cell(nvol,1);

data.Azo = cell(nvol,1);


% Read the field files

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
  field_file = [filepath filesep '.' name '.A'];
  fid        = fopen(field_file,'r',machine_format);
  if (fid > 0)
    % Read through Mvol, Mpol, Ntor, mn, Nfp, im(1:mn), in(1:mn)
    fread(fid,1,spacer_format);
    fread(fid,5,int_format);  % Mvol, Mpol, Ntor, mn, Nfp
    fread(fid,2,spacer_format);
    fread(fid,mn,int_format); % im
    fread(fid,2,spacer_format);
    fread(fid,mn,int_format); % in
    fread(fid,1,spacer_format);
  else
   disp([' - File "' field_file '" does not exist']); break;
  end  

  for v=1:nvol
   % Get radial resolution
   fread(fid,1,spacer_format);
   Lrad = fread(fid,1,int_format);
   fread(fid,1,spacer_format);
   data.Lrad(v,1) = Lrad;
   % Allocate data for the fields
   data.Ate{v} = zeros(Lrad+1,mn);
   data.Aze{v} = zeros(Lrad+1,mn);
   data.Ato{v} = zeros(Lrad+1,mn);
   data.Azo{v} = zeros(Lrad+1,mn);
   for i=1:mn
    fread(fid,1,spacer_format);
    data.Ate{v}(:,i) = fread(fid,Lrad+1,float_format); % Ate
    fread(fid,2,spacer_format);
    data.Aze{v}(:,i) = fread(fid,Lrad+1,float_format); % Aze
    fread(fid,2,spacer_format);
    data.Ato{v}(:,i) = fread(fid,Lrad+1,float_format); % Ato
    fread(fid,2,spacer_format);
    data.Azo{v}(:,i) = fread(fid,Lrad+1,float_format); % Azo
    fread(fid,1,spacer_format);
   end
  end
  
  success=1;
    
catch

  if(triedallform==1)      
   disp(' - Could not read field file'); break;
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
fdata = data;

