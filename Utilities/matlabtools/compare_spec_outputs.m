function [Dabs, Drel, df] = compare_spec_outputs(data1,data2)

% 
% COMPARE_SPEC_OUTPUTS( DATA1, DATA2 )
% ====================================
%
% Compares the interface geometry of two spec outputs
% Two outputs are considered "the same" if df ~ 1e-15 (or less)
%
% INPUT
% -----
% - data1: data obtained via read_spec( fname1 )
% - data2: data obtained via read_spec( fname2 )
%
% OUTPUT
% ------
% - Dabs   : absolute maximum distance
% - Drel   : relative maximum distance
% - df     : expected change in force-balance if one output is used as input for the other run
%
% written by J.Loizu (08.2017)
% modified by J.Loizu (10.2017)
%
% TODO: UPDATE - IS THIS STILL CORRECT ?

% 
% Rmn1        = h5read(fname1,'/Rbc');
% Zmn1        = h5read(fname1,'/Zbs');
% 
%   
% Rmn2        = h5read(fname2,'/Rbc');
% Zmn2        = h5read(fname2,'/Zbs');

%data        = read_spec(fname1);
Rmn1        = data1.output.Rbc;
Zmn1        = data1.output.Zbs;

%data        = read_spec(fname2);
Rmn2        = data2.output.Rbc;
Zmn2        = data2.output.Zbs;


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
if Drel==0
    df = 0;
else
    df = (Dabs^2)/Drel;
end

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

