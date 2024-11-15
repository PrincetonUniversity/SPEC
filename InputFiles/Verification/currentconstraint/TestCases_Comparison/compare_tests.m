

local_array = {};
global_array = {};
master_array = {};
ii=1;
fid = fopen('input_local');
tline = fgetl(fid);
while ischar(tline)
    local_array{ii} = [tline(1:end-2), 'h5'];
    global_array{ii} = [tline(1:6), '3', tline(8:end-2), 'h5'];
    master_array{ii} = [tline(1:9), '.master', tline(10:end-2), 'sp.h5'];
    tline = fgetl(fid);
    ii = ii+1;
end
fclose(fid);

for ii = 1:length(local_array)
    
    disp(['COMPARISON OF ', global_array{ii}])
    disp('------------------------------')
    
    data1 = read_spec(master_array{ii});
    data2 = read_spec(global_array{ii});
    
    try
        compare_spec_outputs(master_array{ii}, global_array{ii});
        compare_spec_outputs2(data1, data2, 1E-12);
        disp([newline,' '])
    catch err
        disp(err.message)
    end
    disp(' ')
    disp(' ')
end


for ii = 1:length(local_array)
    
    disp(['COMPARISON OF ', local_array{ii}])
    disp('------------------------------')
    data1 = read_spec(master_array{ii});
    data2 = read_spec(local_array{ii});
    try
        compare_spec_outputs(local_array{ii}, master_array{ii});
        compare_spec_outputs2(data1, data2, 1E-12);
        disp([newline,' '])
    catch err
        disp(err.message)
    end
    disp(' ')
    disp(' ')
end
