function Hmatrix = read_spec_hessian(filename)

%
% READ_SPEC_HESSIAN( FILENAME )
% =============================
%
% Reads Hessian matrix using output from SPEC
%
% INPUT
% -----
% - filename : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%
% OUTPUT
% ------
% - Hmatrix    : Hessian matrix, which can be fed into several routines for analyzing and ploting
%
% written by J.Loizu (2017)
% updated by J.Loizu (2020)


global machform;

machform = 's'; 


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
  hessian_file = ['.' filename(1:length(filename)-6) '.GF.ma'];
  fid        = fopen(hessian_file,'r',machine_format);
  if (fid > 0)
    % Read through NGdof
    fread(fid,1,spacer_format);
    ngdof = fread(fid,1,int_format); % NGdof 
    fread(fid,1,spacer_format);  
  else
   disp(' - File does not exist'); break;
  end  

  fread(fid,1,spacer_format);  
  Hmatrix = fread(fid,[ngdof ngdof],float_format);
  fread(fid,1,spacer_format);
 
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
  if (fid ~= -1)
   fclose(fid);
   fid = -1;
  end  
end

end

if (fid ~= -1) 
    fclose(fid);
end

