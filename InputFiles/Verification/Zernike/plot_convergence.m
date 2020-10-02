MN    = [11:2:35];


mn    = length(MN);
Berravg = zeros([1,mn]);
Berrmax = zeros([1,mn]);
cd vac
for i=1:mn
 
 fname     = strcat('w7x_vac_M',num2str(MN(i)),'N',num2str(MN(i)),'L',num2str(MN(i)+4),'.sp.h5');
 fdata     = read_spec_field(fname);

 Berravg(i) = sum(fdata.beltramierror(1,4:6))/3;
 Berrmax(i) = max(fdata.beltramierror(1,7:9));

end
cd ..
figure
semilogy(MN,Berravg,'r*-')
hold on
semilogy(MN,Berrmax,'kx-')


cd taylor
for i=1:mn
 
 fname     = strcat('w7x_tay_M',num2str(MN(i)),'N',num2str(MN(i)),'L',num2str(MN(i)+4),'.sp.h5');
 fdata     = read_spec_field(fname);

 Berravg(i) = sum(fdata.beltramierror(1,4:6))/3;
 Berrmax(i) = max(fdata.beltramierror(1,7:9));

end
cd ..

semilogy(MN,Berravg,'rs--')
semilogy(MN,Berrmax,'k+--')

xlabel('M=N=L-4');
ylabel('Error');
legend('Vac Avg', 'Vac Max', 'Taylor Avg', 'Taylor Max');