function data = read_hdf5(filename)

%
% OUTDATED ?
%
% READ_HDF5( FILENAME )
% =====================
%
%   The READ_HDF5 function reads an HDF5 file and returns the contents of
%   that file as the fields of a structure.  Groups are treated as elements
%   of their parent structure.  If the file file cannot be opened a -1 is
%   returned.
%
%   Example
%       data=read_hdf5('input.h5');
%
%   Version 1.0
%   Maintained by: Samuel Lazerson (lazerson@pppl.gov)
%   Date  05/02/2012


% Try to read the file first
try
    data_info = h5info(filename,'/');
catch h5info_error
    data=-1;
    disp(['ERROR: Opening HDF5 File: ' filename]);
    disp(['  -identifier: ' h5info_error.identifier]);
    disp(['  -message:    ' h5info_error.message]);
    disp('      For information type:  help read_hdf5');
    return
end

ngroups = length(data_info.Groups);
nvars = length(data_info.Datasets);
% Get root datasets
for i = 1: nvars
    data.([data_info.Datasets(i).Name]) = h5read(filename,['/' data_info.Datasets(i).Name]);
    natts = length(data_info.Datasets(i).Attributes);
    for j=1:natts
        data.([data_info.Datasets(i).Attributes(j).Name]) = data_info.Datasets(i).Attributes(j).Value{1};
    end
    % Fix by Weatherby,Gerard <gweatherby@uchc.edu> 10/17/13
    %for j=2:natts 
    %    data.([data_info.Datasets(i).Attributes(j).Name]) = data_info.Datasets(i).Attributes(j).Value{1};
    %end
    %if natts == 1
    %    data.([data_info.Datasets(i).Attributes(1).Name]) = data_info.Datasets(i).Attributes(1);
    %end  
end

% Get each subgroup
if ngroups > 0
    for i = 1 : ngroups
    data.([data_info.Groups(i).Name]) = getGroup(['/' data_info.Groups(i).Name]);
    end
end
return
end

function data = getGroup(root)
data_info = h5info(filename,root);
ngroups = length(data_info.Groups);
nvars = length(data_info.Datasets);
% Get root datasets
for i = 1: nvars
    data.([data_info.Datasets(i).Name]) = h5read(filename,[root data_info.Datasets(i).Name]);
    natts = length(data_info.Datasets(i).Attributes);
    for j=1:natts
        data.([data_info.Datasets(i).Attributes(j).Name]) = data_info.Datasets(i).Attributes(j).Value{1};
    end
end

% Get each subgroup
if ngroups > 0
    for i = 1 : ngroups
    data.([data_info.Groups(i).Name]) = getGroup(['/' data_info.Groups(i).Name]);
    end
end
return
end
