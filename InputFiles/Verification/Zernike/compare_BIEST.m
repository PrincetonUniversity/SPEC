% compare vac

[Bx,By,Bz]=read_BIEST_vector('tmp-B0.data');
[x,y,z]=read_BIEST_vector('Svec.data');

NtB = 14*32;
NzB = 14*5*32;

Bx = reshape(Bx(:),NtB,NzB);
By = reshape(By(:),NtB,NzB);
Bz = reshape(Bz(:),NtB,NzB);

BxB = Bx(:,1:NzB/5);
ByB = By(:,1:NzB/5);
BzB = Bz(:,1:NzB/5);


MN    = [11:2:35];

mn    = length(MN);
maxerr = zeros([1,mn]);

for i=1:mn
 cd vac
 fname     = strcat('w7x_vac_M',num2str(MN(i)),'N',num2str(MN(i)),'L',num2str(MN(i)+4),'.sp.h5');
 gdata     = read_spec_grid(fname);
 cd ..
 
 NtSPEC = gdata.Nt;
 NzSPEC = gdata.Nz;
 
 R = reshape(gdata.Rij(1,:,2),NtSPEC,NzSPEC);
 Z = reshape(gdata.Zij(1,:,2),NtSPEC,NzSPEC);

 
 BrS = reshape(gdata.BR(1,:,2),NtSPEC,NzSPEC);
 BpS = reshape(gdata.Bp(1,:,2),NtSPEC,NzSPEC) .* R;
 BzS = reshape(gdata.BZ(1,:,2),NtSPEC,NzSPEC);
 
 
 BrS = resample2D(BrS,NtB,NzB/5);
 BpS = resample2D(BpS,NtB,NzB/5);
 BzS = resample2D(BzS,NtB,NzB/5);
 
 zeta = linspace(0,2*pi/double(gdata.Nfp),NzB/5+1);
 zeta = repmat(zeta(1:end-1), [NtB,1]); 
 BxS = BrS .* cos(zeta) - BpS .* sin(zeta);
 ByS = BrS .* sin(zeta) + BpS .* cos(zeta);
 
 
 maxerr(i) = max(abs(BzS(:) - BzB(:)));
 maxerr(i) = max(max(abs(BxS(:) - BxB(:))), maxerr(i));
 maxerr(i) = max(max(abs(ByS(:) - ByB(:))), maxerr(i));

end

maxB = max([max(abs(Bx(:))), max(abs(By(:))),max(abs(Bz(:)))]);
maxerr = maxerr / maxB;

figure
semilogy(MN,maxerr,'r*-')
xlabel('M=N=L-4');
ylabel('e_{max}');

% Taylor

[Bx,By,Bz]=read_BIEST_vector('tmp-B1.data');

NtB = 14*28;
NzB = 14*5*28;

Bx = reshape(Bx(:),NtB,NzB);
By = reshape(By(:),NtB,NzB);
Bz = reshape(Bz(:),NtB,NzB);

BxB = Bx(:,1:NzB/5);
ByB = By(:,1:NzB/5);
BzB = Bz(:,1:NzB/5);


MN    = [11:2:35];

mn    = length(MN);
maxerr = zeros([1,mn]);

for i=1:mn
 cd taylor
 fname     = strcat('w7x_tay_M',num2str(MN(i)),'N',num2str(MN(i)),'L',num2str(MN(i)+4),'.sp.h5');
 gdata     = read_spec_grid(fname);
 cd ..
 
 NtSPEC = gdata.Nt;
 NzSPEC = gdata.Nz;
 
 R = reshape(gdata.Rij(1,:,2),NtSPEC,NzSPEC);
 Z = reshape(gdata.Zij(1,:,2),NtSPEC,NzSPEC);

 
 BrS = reshape(gdata.BR(1,:,2),NtSPEC,NzSPEC);
 BpS = reshape(gdata.Bp(1,:,2),NtSPEC,NzSPEC) .* R;
 BzS = reshape(gdata.BZ(1,:,2),NtSPEC,NzSPEC);
 
 
 BrS = resample2D(BrS,NtB,NzB/5);
 BpS = resample2D(BpS,NtB,NzB/5);
 BzS = resample2D(BzS,NtB,NzB/5);
 
 zeta = linspace(0,2*pi/double(gdata.Nfp),NzB/5+1);
 zeta = repmat(zeta(1:end-1), [NtB,1]); 
 BxS = BrS .* cos(zeta) - BpS .* sin(zeta);
 ByS = BrS .* sin(zeta) + BpS .* cos(zeta);
 
 
 maxerr(i) = max(abs(BzS(:) - BzB(:)));
 maxerr(i) = max(max(abs(BxS(:) - BxB(:))), maxerr(i));
 maxerr(i) = max(max(abs(ByS(:) - ByB(:))), maxerr(i));

end

maxB = max([max(abs(Bx(:))), max(abs(By(:))),max(abs(Bz(:)))]);
maxerr = maxerr / maxB;

hold on
semilogy(MN,maxerr,'kx--')
legend('\mu=0', '\mu=1')
