% Comparison between SPEC and Freidberg for the relation between vertical field and plasma position and size: BV=BV(R,a)

froot    = 'wiretest_phiedge1_BV0.00';

bvs      = [35 37 39 41 43 45 47 49 51 53 55];

BV_SPEC  = 0.0001*bvs;

mu       = 0.1;

muI      = 0.1;

for i=1:length(bvs)

 fname  = strcat(froot, num2str(bvs(i)), '.sp.h5');
 
 pdata  = read_spec_poincare(fname);

 R(i)   = pdata.Rbc(1,2);

 ar     = sum(pdata.Rbc(2:end,2));

 az     = abs(sum(pdata.Zbs(2:end,2)));

 a(i)   = sqrt(ar*az);

 aoR(i) = a(i)/R(i);

end

li     = 1 - besselj(0,mu*a).*besselj(2,mu*a)./besselj(1,mu*a).^2;

BV_FB  = (muI./(4*pi*R)).*((li-3)./2 + log(8.*R./a)); 


figure
hold on

curve1 = BV_SPEC.*(1 - aoR.^2);
curve2 = BV_SPEC.*(1 + aoR.^2);
plot(BV_SPEC, curve1, 'k');
plot(BV_SPEC, curve2, 'k');
x2 = [BV_SPEC, fliplr(BV_SPEC)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'y');

plot(BV_SPEC, BV_SPEC, 'k--','LineWidth',2);

plot(BV_SPEC, BV_FB, 'rx','LineWidth',4,'MarkerSize',6)

xlabel('B_V^{SPEC}');
ylabel('B_V^{theory}');






