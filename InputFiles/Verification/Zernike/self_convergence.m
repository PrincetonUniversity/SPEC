% compare vac

MN    = [11:2:35];

higestMN = 39;

mn    = length(MN);
maxerr = zeros([1,mn]);


cd vac
fname     = strcat('w7x_vac_M',num2str(higestMN),'N',num2str(higestMN),'L',num2str(higestMN+4),'.sp.h5');
gdata     = read_spec_grid(fname);

NtB = gdata.Nt;
NzB = gdata.Nz;

BrB = reshape(gdata.BR(1,:,2),NtB,NzB);
BpB = reshape(gdata.Bp(1,:,2).* gdata.Rij(1,:,2),NtB,NzB) ;
BzB = reshape(gdata.BZ(1,:,2),NtB,NzB);

zeta = linspace(0,2*pi/double(gdata.Nfp),NzB+1);
zeta = repmat(zeta(1:end-1), [NtB,1]); 

BxB = BrB .* cos(zeta) - BpB .* sin(zeta);
ByB = BrB .* sin(zeta) + BpB .* cos(zeta);
cd ..


for i=1:mn
 cd vac
 fname     = strcat('w7x_vac_M',num2str(MN(i)),'N',num2str(MN(i)),'L',num2str(MN(i)+4),'.sp.h5');
 gdata     = read_spec_grid(fname);
 
 NtSPEC = gdata.Nt;
 NzSPEC = gdata.Nz;
 cd ..
 R = reshape(gdata.Rij(1,:,2),NtSPEC,NzSPEC);
 Z = reshape(gdata.Zij(1,:,2),NtSPEC,NzSPEC);

 
 BrS = reshape(gdata.BR(1,:,2),NtSPEC,NzSPEC);
 BpS = reshape(gdata.Bp(1,:,2),NtSPEC,NzSPEC) .* R;
 BzS = reshape(gdata.BZ(1,:,2),NtSPEC,NzSPEC);
 
 
 BrS = resample2D(BrS,NtB,NzB);
 BpS = resample2D(BpS,NtB,NzB);
 BzS = resample2D(BzS,NtB,NzB);
 
 zeta = linspace(0,2*pi/double(gdata.Nfp),NzB+1);
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
hold on


cd taylor
fname     = strcat('w7x_tay_M',num2str(higestMN),'N',num2str(higestMN),'L',num2str(higestMN+4),'.sp.h5');
gdata     = read_spec_grid(fname);

NtB = gdata.Nt;
NzB = gdata.Nz;

BrB = reshape(gdata.BR(1,:,2),NtB,NzB);
BpB = reshape(gdata.Bp(1,:,2).* gdata.Rij(1,:,2),NtB,NzB) ;
BzB = reshape(gdata.BZ(1,:,2),NtB,NzB);

zeta = linspace(0,2*pi/double(gdata.Nfp),NzB+1);
zeta = repmat(zeta(1:end-1), [NtB,1]); 

BxB = BrB .* cos(zeta) - BpB .* sin(zeta);
ByB = BrB .* sin(zeta) + BpB .* cos(zeta);
cd ..


for i=1:mn
 cd taylor
 fname     = strcat('w7x_tay_M',num2str(MN(i)),'N',num2str(MN(i)),'L',num2str(MN(i)+4),'.sp.h5');
 gdata     = read_spec_grid(fname);
 
 NtSPEC = gdata.Nt;
 NzSPEC = gdata.Nz;
 cd ..
 R = reshape(gdata.Rij(1,:,2),NtSPEC,NzSPEC);
 Z = reshape(gdata.Zij(1,:,2),NtSPEC,NzSPEC);

 
 BrS = reshape(gdata.BR(1,:,2),NtSPEC,NzSPEC);
 BpS = reshape(gdata.Bp(1,:,2),NtSPEC,NzSPEC) .* R;
 BzS = reshape(gdata.BZ(1,:,2),NtSPEC,NzSPEC);
 
 
 BrS = resample2D(BrS,NtB,NzB);
 BpS = resample2D(BpS,NtB,NzB);
 BzS = resample2D(BzS,NtB,NzB);
 
 zeta = linspace(0,2*pi/double(gdata.Nfp),NzB+1);
 zeta = repmat(zeta(1:end-1), [NtB,1]); 
 BxS = BrS .* cos(zeta) - BpS .* sin(zeta);
 ByS = BrS .* sin(zeta) + BpS .* cos(zeta);
 
 maxerr(i) = max(abs(BzS(:) - BzB(:)));
 maxerr(i) = max(max(abs(BxS(:) - BxB(:))), maxerr(i));
 maxerr(i) = max(max(abs(ByS(:) - ByB(:))), maxerr(i));

end



maxB = max([max(abs(Bx(:))), max(abs(By(:))),max(abs(Bz(:)))]);
maxerr = maxerr / maxB;
semilogy(MN,maxerr,'bx-')
