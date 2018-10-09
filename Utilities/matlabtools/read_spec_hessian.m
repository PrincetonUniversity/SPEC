function hdata = read_spec_hessian(filename)

% Reads Hessian matrix using output from SPEC
%
% INPUT
% - filename : path to the hdf5 output file (e.g. 'testcase.h5')
%
% OUTPUT
% - hdata    : contains the Hessian matrix, which can be fed into several routines for analyzing and ploting
%
% written by J.Loizu (2017)


global machform;

machform = 's'; 

data     = read_hdf5(filename);


% Read the hessian file

machine_format =  machform;  % needs to be 's' (intel) or 'a' (gnu); format is swapped if an error occurs 
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
  hessian_file = ['.' filename(1:length(filename)-3) '.GF'];
  fid        = fopen(hessian_file,'r',machine_format);
  if (fid > 0)
    % Read through Nvol, Mpol, Ntor, NGdof
    fread(fid,1,spacer_format);
    fread(fid,3,int_format);         % Nvol, Mpol, Ntor
    ngdof = fread(fid,1,int_format); % NGdof 
    fread(fid,1,spacer_format);  
  else
   disp(' - File does not exist'); break;
  end  

  fread(fid,1,spacer_format);  
  Hmatrix = fread(fid,[ngdof ngdof],float_format);
  fread(fid,1,spacer_format);

  data.Hmatrix = Hmatrix;
 
  success=1;
    
catch

  if(triedallform==1)      
   disp(' - Could not read hessian file'); break;
  end
  if(machform~='s')
   machform     ='s'; 
   triedallform = 1;
  else
   machform     ='a';
   triedallform = 1;
  end
  
end

end


fclose(fid);


%return output
hdata = data;
