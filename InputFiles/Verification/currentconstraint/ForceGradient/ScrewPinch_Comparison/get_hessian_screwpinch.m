function hessian = get_hessian_screwpinch(Nvol, tflux, pflux, mu, R, ns)


switch Nvol
    case 2
        a     = R(3);
        R     = R(2);

        struc1 = ScrewPinch_Solution_Beltrami(0, R, ns, pflux(1), tflux(1), mu(1), 1, 1, 1);
        struc2 = ScrewPinch_Solution_Beltrami(R, a, ns, pflux(2), tflux(2), mu(2), 2, 0, 1);

        dpfluxdR = (struc1.dBtdR - struc2.dBtdR) / struc2.dBtdpf;

        A = struc2.BdBdR ;
        B = struc2.BdBdpf * dpfluxdR;
        C = - struc1.BdBdR ;

        hessian = A+B+C;
        hessian = hessian * R;
    case 3
        R = R(2:end);

        ivol=1;
        struc1_ou1 = ScrewPinch_Solution_Beltrami(0   , R(1), ns, pflux(1), tflux(1), mu(1), 1, 1, ivol);
        struc2_in1 = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 0, ivol);
        struc2_ou1 = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 1, ivol);
        struc3_in1 = ScrewPinch_Solution_Beltrami(R(2), R(3), ns, pflux(3), tflux(3), mu(3), 3, 0, ivol);

        dpfdR = get_dpfdR_N3(tflux, pflux, mu, R, ns);

        hessian = zeros(2,2);
        hessian(1,1) = struc2_in1.BdBdR + struc2_in1.BdBdpf * dpfdR(1,2) - struc1_ou1.BdBdR;
        hessian(2,1) = struc3_in1.BdBdR + struc3_in1.BdBdpf * dpfdR(1,3) - struc2_ou1.BdBdR - struc2_ou1.BdBdpf * dpfdR(1,2);
        hessian(:,1) = hessian(:,1) * R(1);

        ivol=2;
        struc1_ou2 = ScrewPinch_Solution_Beltrami(0   , R(1), ns, pflux(1), tflux(1), mu(1), 1, 1, ivol);
        struc2_in2 = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 0, ivol);
        struc2_ou2 = ScrewPinch_Solution_Beltrami(R(1), R(2), ns, pflux(2), tflux(2), mu(2), 2, 1, ivol);
        struc3_in2 = ScrewPinch_Solution_Beltrami(R(2), R(3), ns, pflux(3), tflux(3), mu(3), 3, 0, ivol);


        hessian(1,2) = struc2_in2.BdBdR + struc2_in2.BdBdpf * dpfdR(2,2) - struc1_ou2.BdBdR;
        hessian(2,2) = struc3_in2.BdBdR + struc3_in2.BdBdpf * dpfdR(2,3) - struc2_ou2.BdBdR - struc2_ou2.BdBdpf * dpfdR(2,2);
        hessian(:,2) = hessian(:,2) * R(2);

    otherwise
        
        disp('Unsupported number of volumes')
end


