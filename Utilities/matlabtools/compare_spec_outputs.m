function compare_spec_outputs(fname1,fname2)

% Compares the interface geometry of two spec outputs
% Two outputs are considered "the same" if df ~ 1e-15 (or less)
%
% INPUT
% - fname1 : path to the hdf5 output file #1
% - fname2 : path to the hdf5 output file #2
%
% OUTPUT
% - Dabs   : absolute maximum distance
% - Drel   : relative maximum distance
% - df     : expected change in force-balance if one output is used as input for the other run
%
% written by J.Loizu (08.2017)
% modified by J.Loizu (10.2017)


Rmn1        = h5read(fname1,'/Rbc');
Zmn1        = h5read(fname1,'/Zbs');

  
Rmn2        = h5read(fname2,'/Rbc');
Zmn2        = h5read(fname2,'/Zbs');

  
maxR        = max(max(abs(Rmn1-Rmn2)));

maxZ        = max(max(abs(Zmn1-Zmn2)));

maxrelR     = max(max(abs(Rmn1-Rmn2)./abs(Rmn1)));

if(Zmn1 == 0)
maxrelZ     = 0;
else
maxrelZ     = max(max(abs(Zmn1-Zmn2)./abs(Zmn1)));
end

%Absolute maximum distance
Dabs        = max(maxR,maxZ);

%Relative maximum distance
Drel        = max(maxrelR,maxrelZ);

%Estimate for df~R*dR~(Dabs/Drel)*Dabs
df          = (Dabs^2)/Drel;

%Output 
display(['Dabs:            ' num2str(Dabs)]);
display(['Drel:            ' num2str(Drel)]);
display(['Estimate for df: ' num2str(df)]);

if(Dabs==0)
display(['The two outputs are exaclty the same']);
elseif(df < 1e-14)
display(['The two outputs can be considered the same']);
else
display(['The two outputs are not the same'])
end

