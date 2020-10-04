function dpfdR = get_dpfdR_N3(tflux, pflux, mu, R, ns)

    dpfdR = zeros(2,3);





    % First interface derivative
    ivol=1;
    struc1_ou = ScrewPinch_Solution_Beltrami(0   , R(1), ns, pflux(1), tflux(1), mu(1), 1, 1, ivol);
    struc2_in = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 0, ivol);
    struc2_ou = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 1, ivol);
    struc3_in = ScrewPinch_Solution_Beltrami(R(2), R(3), ns, pflux(3), tflux(3), mu(3), 3, 0, ivol);
    
    
    
    M = zeros(2,2);
    M(1,1) = struc2_in.dBtdpf;
    M(2,1) =-struc2_ou.dBtdpf;
    %M(2,1) = 0;
    M(2,2) = struc3_in.dBtdpf;
    
    
    rhs = zeros(2,1);
    rhs(1) = struc1_ou.dBtdR - struc2_in.dBtdR;
    rhs(2) = struc2_ou.dBtdR - 0.0;

    x = M\rhs;
    
    dpfdR(1,1  ) = 0;
    dpfdR(1,2:3) = x(1:2);

    %disp(newline)
    %for ii=1:3
    %    disp(['Derivative of poloidal flux in vol=',num2str(ii), ' w.r.t R(1): ', num2str(dpfdR(1,ii), 16)])
    %end

    % Second interface derivative
    ivol=2;
    struc1_ou = ScrewPinch_Solution_Beltrami(0   , R(1), ns, pflux(1), tflux(1), mu(1), 1, 1, ivol);
    struc2_in = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 0, ivol);
    struc2_ou = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 1, ivol);
    struc3_in = ScrewPinch_Solution_Beltrami(R(2), R(3), ns, pflux(3), tflux(3), mu(3), 3, 0, ivol);
    
    rhs = zeros(2,1);
    rhs(1) = struc1_ou.dBtdR - struc2_in.dBtdR;
    rhs(2) = struc2_ou.dBtdR - struc3_in.dBtdR;

    x = M\rhs;

    dpfdR(2,1  ) = 0;
    dpfdR(2,2:3) = x(1:2);

    %disp(newline)
    %for ii=1:3
    %    disp(['Derivative of poloidal flux in vol=',num2str(ii), ' w.r.t R(2): ', num2str(dpfdR(2,ii), 16)])
    %end
    
end
  