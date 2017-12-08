function idata = read_spec_iota(filename)

% Reads rotational transform data from field-line-tracing using output from SPEC
%
% INPUT
% - filename : path to the hdf5 output file (e.g. 'testcase.sp.h5')
%
% OUTPUT
% - idata    : contains all the transform data, which can be fed into several routines for analyzing and ploting
%
% written by J.Loizu (2017)


global machform;

machform = 's';

data     = read_hdf5(filename);

nvol     = double(data.Nvol);


% Read the transform files

machine_format =  machform;  % needs to be 's' (intel) or 'a' (gnu); format is swapped if an error occurs 
triedallform   =  0;         % whether all formats have been tried (1) or not (0) 
success        =  0;
int_format     = 'int32';
float_format   = 'float64';
spacer_format  = 'int32';

j=1;

for i=1:nvol
    success      = 0;
    triedallform = 0;
    while(success==0)
    try
        if(machine_format ~= machform)
         machine_format =  machform;  % update value
        end
        transform_file=['.' filename(1:length(filename)-3) '.t.' num2str(i,'%4.4i') '.dat'];
        fid = fopen(transform_file,'r',machine_format);
        if (fid > 0)
            % Get subgrid
            fread(fid,1,spacer_format);
            npts=fread(fid,1,int_format);
            fread(fid,1,spacer_format);
	    if(npts>1e4), npts=-1; end;  % dirty fix to force switching of binary format 
            % Get Iota edges
            fread(fid,1,spacer_format);
            iotae=fread(fid,2,float_format);
            fread(fid,1,spacer_format);
            % Get Profile iota(s)
            fread(fid,1,spacer_format);
            out=fread(fid,[2 npts],float_format);
	    nptraj=size(out,2);
            data.sarr(j:j+nptraj-1) = out(1,:);
            data.iota(j:j+nptraj-1) = out(2,:);
	    j=j+nptraj;
            success = 1;
	    fclose(fid);
        else
        disp(' - File does not exist'); break
        end        
    catch
        if(triedallform==1)      
         disp(' - Could not read transform file'); break;
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
end


%return output
idata = data;

