! An analytic helical coordinate map with different plasma/vacuum radial gaps.
! Check the actual surfacecurrent routine, including both Fourier parities.
program test_surfacecurrent_first
    use constants, only: zero, one, half
    use inputlist, only: Igeometry, Nvol, Lrad
    use allglobal, only: Mvol, mn, im, in, YESstellsym, NOTstellsym, &
        first_free_bound, iRbc, iZbs, iRbs, iZbc, Ate, Aze, Ato, Azo, TT
    implicit none
    integer :: sym, initial, v, j
    real :: radius(0:4), eps, eta, theta, phi, q, c, s, dr, r
    real :: rt, zt, rp, zp, rs, zs, jac, gu, gv, bt, bz
    real :: et(3), ep(3), expected(3), got(3), position(3), expected_position(3)

    Igeometry = 3
    Nvol = 3
    Mvol = 4
    mn = 2
    Lrad = 2
    allocate(im(mn), in(mn), TT(0:2,0:1,0:3))
    allocate(iRbc(mn,0:Mvol), iZbs(mn,0:Mvol))
    allocate(iRbs(mn,0:Mvol), iZbc(mn,0:Mvol))
    allocate(Ate(Mvol,0:0,mn), Aze(Mvol,0:0,mn))
    allocate(Ato(Mvol,0:0,mn), Azo(Mvol,0:0,mn))
    do v = 1, Mvol
        do j = 1, mn
            allocate(Ate(v,0,j)%s(0:2), Aze(v,0,j)%s(0:2))
            allocate(Ato(v,0,j)%s(0:2), Azo(v,0,j)%s(0:2))
            Ate(v,0,j)%s = zero
            Aze(v,0,j)%s = zero
            Ato(v,0,j)%s = zero
            Azo(v,0,j)%s = zero
        enddo
    enddo
    im = [0, 1]
    in = [0, 10]
    TT = zero
    TT(1,:,1) = one ! T_1'(s)=1 at both annular faces.
    radius = [0.0, 0.2, 0.4, 0.7, 0.8]
    Ate(Nvol,0,1)%s(1) = one
    Aze(Nvol,0,1)%s(1) = -0.25
    Ate(Mvol,0,1)%s(1) = 4.0
    Aze(Mvol,0,1)%s(1) = -2.0
    Ato(Nvol,0,2)%s(0:1) = [-0.03, 0.03]
    Azo(Nvol,0,2)%s(0:1) = [-0.02, 0.02]
    Ato(Mvol,0,2)%s(0:1) = [0.4, 0.4]
    Azo(Mvol,0,2)%s(0:1) = [0.5, 0.5]
    theta = 0.37
    phi = 0.02
    q = theta - 10 * phi
    c = cos(q)
    s = sin(q)
    do sym = 0, 1
        YESstellsym = sym.eq.1
        NOTstellsym = .not.YESstellsym
        eps = zero
        eta = zero
        if (NOTstellsym) then
            eps = 0.1
            eta = 0.04
        endif
        iRbc = zero
        iZbs = zero
        iRbs = zero
        iZbc = zero
        iRbc(1,:) = 3.6
        iRbc(2,:) = radius
        iZbs(2,:) = -radius
        iRbs(2,:) = eps * radius
        iZbc(2,:) = eta * radius
        r = 3.6 + radius(Nvol) * (c + eps*s)
        rt = radius(Nvol) * (-s + eps*c)
        zt = radius(Nvol) * (-c - eta*s)
        rp = -10 * rt
        zp = -10 * zt
        et = [rt*cos(phi), rt*sin(phi), zt]
        ep = [rp*cos(phi)-r*sin(phi), rp*sin(phi)+r*cos(phi), zp]
        expected_position = [r*cos(phi), r*sin(phi), radius(Nvol)*(-s+eta*c)]
        do initial = 0, 1
            first_free_bound = initial.eq.1
            dr = half * (radius(Mvol)-radius(Nvol))
            gu = 2.0
            gv = 4.0
            if (first_free_bound) then
                dr = half * (radius(Nvol)-radius(Nvol-1))
                gu = 0.25
                gv = one
            endif
            if (NOTstellsym) then
                if (first_free_bound) then
                    gu = gu - 0.02*s
                    gv = gv + 0.03*s
                else
                    gu = gu - 0.5*s
                    gv = gv + 0.4*s
                endif
            endif
            rs = dr * (c + eps*s)
            zs = dr * (-s + eta*c)
            jac = r * (zs*rt-rs*zt)
            bt = (gu*dot_product(et,et)+gv*dot_product(et,ep))/jac
            bz = (gu*dot_product(et,ep)+gv*dot_product(ep,ep))/jac
            expected = bz*et-bt*ep
            call surfacecurrent(theta, phi, position, got)
            if (maxval(abs(got-expected)).gt.1.e-11) error stop 'surface current mismatch'
            if (maxval(abs(position-expected_position)).gt.1.e-13) error stop 'position mismatch'
        enddo
    enddo
end program test_surfacecurrent_first
