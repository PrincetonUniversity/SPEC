MN    = [4 6 8 10 12 14];

files = {'vacuum.M04N04L6.sp.h5' 'vacuum.M06N06L6.sp.h5' 'vacuum.M08N08L6.sp.h5' 'vacuum.M10N10L6.sp.h5' 'vacuum.M12N12L8.sp.h5' 'vacuum.M14N14L8.sp.h5'};

mn    = length(MN);

for i=1:mn
 i

 fname     = files{i};
 fdata     = read_spec_field(fname);

 Berr1s(i) = fdata.beltramierror(1,1);
 Berr1t(i) = fdata.beltramierror(1,2);
 Berr1z(i) = fdata.beltramierror(1,3);

 Berr2s(i) = fdata.beltramierror(2,1);
 Berr2t(i) = fdata.beltramierror(2,2);
 Berr2z(i) = fdata.beltramierror(2,3);
end

figure
hold on

plot(MN,log10(Berr1s),'r*')
plot(MN,log10(Berr1t),'k*')
plot(MN,log10(Berr1z),'b*')
plot(MN,log10(Berr2s),'ro')
plot(MN,log10(Berr2t),'ko')
plot(MN,log10(Berr2z),'bo')

