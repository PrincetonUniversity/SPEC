% Calculates Hessian from MRxMHD analytical linear theory 
%
% INPUTS
%
% mu2
% tf2
% pf3
% R3
%
% OUTPUTS
%
% K1, K2, K3
% Delta1, Delta2, Delta3
% Hmatrix
% eval, evec, lambda

Lconstraint = 2; % 0 or 2

mu = [0 mu2 0];
tf = [(1-tf2)/2  tf2  (1-tf2)/2];
pf = [-pf3  0.  pf3];

bcont = 0; % if bcont=1 then poloidal flux pf is adjusted in order to ensure continuity of B

% Calculate Equilibrium geometry from force-balance

x0 = R3*[tf(1) tf(2)];

options = optimoptions('fsolve','Display','none','OptimalityTolerance',1e-15,'FunctionTolerance',1e-15);

x  = fsolve(@(x)force(x,R3,mu,tf,pf,bcont), x0, options);

Delta = [x(1) x(2) R3-x(1)-x(2)];

% Evaluate magnetic helicity

mub  = mu.*Delta/2;

K(1) = helicity(mub(1),tf(1),pf(1));
K(2) = helicity(mub(2),tf(2),pf(2));
K(3) = helicity(mub(3),tf(3),pf(3));

% Evaluate A, C and D terms

D2   = Dfac(tf(2),pf(2),mu(2),Delta(2));
C1p  = Cfac(1  ,tf(1),pf(1),mu(1),Delta(1));
C2p  = Cfac(1  ,tf(2),pf(2),mu(2),Delta(2));
C2m  = Cfac(-1 ,tf(2),pf(2),mu(2),Delta(2));
C3m  = Cfac(-1 ,tf(3),pf(3),mu(3),Delta(3));

if(Lconstraint == 0)     % will calculate hessian at fixed mu
 for i=1:3
  if(mu(i)==0)
   Afac(i) = (tf(i)^2 + pf(i)^2)/(4*pi^2*Delta(i)^3);
  else
   Afac(i) = mu(i)^2*K(i)/(8*pi^2*Delta(i)*tan(mub(i)));
  end
 end
elseif(Lconstraint == 2) % will calculate hessian at fixed helicity
 for i=1:3
  if(mu(i)==0)
   Afac(i) = (tf(i)^2 + pf(i)^2)/(4*pi^2*Delta(i)^3);
  else 
   Afac(i) = mu(i)*K(i)/(4*pi^2*Delta(i)^2);
  end
 end
end

% Evaluate Hessian around the equilibrium

Hmatrix = zeros(4,4);

Hmatrix(1,1) = Afac(1) + Afac(2);
Hmatrix(2,2) = C2m - C1p;
Hmatrix(3,3) = Afac(2) + Afac(3);
Hmatrix(4,4) = C3m - C2p;
Hmatrix(3,1) = -Afac(2);
Hmatrix(4,2) = -D2;

Hmatrix(1,3) = Hmatrix(3,1);
Hmatrix(2,4) = Hmatrix(4,2);

Hmatrix

% Evaluate eigenvalues of Hessian

[evec , eval] = eig(Hmatrix);

lambda = min(diag(eval))

