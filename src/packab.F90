!> \defgroup grp_packing "packing" of Beltrami field solution vector
!>
!> \latexonly
!> \definecolor{Orange}{rgb}{1.0,0.5,0.0}
!> \definecolor{Cerulean}{rgb}{0.0,0.5,1.0}
!> \endlatexonly
!>
!> \file
!> \brief Packs, and unpacks, Beltrami field solution vector; \f${\bf a}\equiv\{A_{\theta,e,i,l}, A_{\zeta,e,i,l}, \mathrm{etc.}\}\f$.

!> \brief Packs and unpacks Beltrami field solution vector.
!> \ingroup grp_packing
!>
!> **construction of "vector" of independent degrees of freedom**
!>
!> <ul>
!> <li> Numerical routines for solving linear equations typically require the unknown, independent degrees of freedom
!>      to be "packed" into a vector, \f${\bf x}\f$. </li>
!> <li> The magnetic field is defined by the independent degrees of freedom in
!>      the Chebyshev-Fourier representation of the vector potential, \f${\color{red} A_{\theta,e,i,l}}\f$ and \f${\color{blue}A_{\zeta, e,i,l}}\f$;
!>      and the non-stellarator-symmetric terms if relevant, \f${\color{Orange}  A_{\theta,o,i,l}}\f$ and \f${\color{Cerulean}A_{\zeta ,o,i,l}}\f$;
!>      and the Lagrange multipliers, \f$a_i\f$, \f$b_i\f$, \f$c_i\f$, \f$d_i\f$, \f$e_i\f$, etc. as required to enforce the constraints:
!>      \f{eqnarray}{ {\bf x} \equiv \{ {\color{red} A_{\theta,e,i,l}},{\color{blue}A_{\zeta, e,i,l}},{\color{Orange}  A_{\theta,o,i,l}},{\color{Cerulean}A_{\zeta ,o,i,l}},a_i,b_i,c_i,d_i,e_i,f_i,g_1,h_1\}.
!>      \f} </li>
!> <li> The "packing" index is assigned in preset() . </li>
!> </ul>
!>
!> @param packorunpack
!> @param lvol
!> @param NN
!> @param solution
!> @param ideriv
subroutine packab(packorunpack, lvol, NN, solution, ideriv)
    use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    use constants, only: zero

    use numerical, only: small

    use fileunits, only: ounit

    use inputlist, only: Wpackab, Lrad

    use cputiming, only: Tpackab

    use allglobal, only: myid, ncpu, cpus, MPI_COMM_SPEC, &
                         mn, im, in, Ate, Aze, Ato, Azo, YESstellsym, NOTstellsym, &
                         TT, YESMatrixFree, &
                         Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                         Lmavalue, Lmbvalue, Lmcvalue, Lmdvalue, Lmevalue, Lmfvalue, Lmgvalue, Lmhvalue

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OPENMP
    USE OMP_LIB
#endif
    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(wp) :: cput, cpui, cpuo = 0 ! cpu time; cpu initial; cpu old; 31 Jan 13;

    character, intent(in) :: packorunpack
    integer, intent(in) :: lvol, NN, ideriv
    real(wp) :: solution(1:NN)

    integer :: ii, ll, id, llrad

    cpui = MPI_WTIME()
    cpuo = cpui
#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    llrad = Lrad(lvol) ! shorthand;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    select case (packorunpack)

    case ('U')

#ifdef DEBUG

        if (Wpackab) then

            if (YESstellsym) then

                ; ii = 1
                do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                    if (id .lt. 1 .or. id .gt. NN) then
                        write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                    end if

                    ; ; id = Aze(lvol, 0, ii)%i(ll); 
                    if (id .lt. 1 .or. id .gt. NN) then
                        write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                    end if

                end do ! end of do ll;
                do ii = 2, mn
                    do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                        if (id .lt. 1 .or. id .gt. NN) then
                            write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                        end if

                        ; ; id = Aze(lvol, 0, ii)%i(ll); 
                        if (id .lt. 1 .or. id .gt. NN) then
                            write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                        end if

                    end do ! end of do ll;
                end do ! end of do ii;

            else ! NOTstellsym;

                ; ii = 1
                do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                    if (id .lt. 1 .or. id .gt. NN) then
                        write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                    end if

                    ; ; id = Aze(lvol, 0, ii)%i(ll); 
                    if (id .lt. 1 .or. id .gt. NN) then
                        write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                    end if

                end do
                do ii = 2, mn
                    do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                        if (id .lt. 1 .or. id .gt. NN) then
                            write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                        end if

                        ; ; id = Aze(lvol, 0, ii)%i(ll); 
                        if (id .lt. 1 .or. id .gt. NN) then
                            write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                        end if

                        ; ; id = Ato(lvol, 0, ii)%i(ll); 
                        if (id .lt. 1 .or. id .gt. NN) then
                            write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                        end if

                        ; ; id = Azo(lvol, 0, ii)%i(ll); 
                        if (id .lt. 1 .or. id .gt. NN) then
                            write (6, '("packab :      fatal : myid=",i3," ; id.lt.1 .or. id.gt.NN ; unpacking illegal subscript ;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "packab : id.lt.1 .or. id.gt.NN : unpacking illegal subscript  ;"
                        end if

                    end do ! end of do ll;
                end do ! end of do ii;

            end if ! end ofif( YESstellsym );

        end if

#endif

        if (YESstellsym) then

            do ii = 1, mn
                do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Ate(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Ate(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Aze(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Aze(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Aze(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; ; Ato(lvol, ideriv, ii)%s(ll) = zero
                    ; ; ; Azo(lvol, ideriv, ii)%s(ll) = zero
                end do ! end of do ll;
            end do ! end of do ii;

            if (YESMatrixFree) then
                do ii = 1, mn
                    ; ; id = Lma(lvol, ii)
                    if (id /= 0) then; Lmavalue(lvol, ii) = solution(id)
                    else; Lmavalue(lvol, ii) = zero
                    end if

                    ; ; id = Lmb(lvol, ii)
                    if (id /= 0) then; Lmbvalue(lvol, ii) = solution(id)
                    else; Lmbvalue(lvol, ii) = zero
                    end if

                    if (ii .gt. 1) then
                        ; ; id = Lme(lvol, ii)
                        if (id /= 0) then; Lmevalue(lvol, ii) = solution(id)
                        else; Lmevalue(lvol, ii) = zero
                        end if

                    else; id = Lmg(lvol, ii)
                        if (id /= 0) then; Lmgvalue(lvol, ii) = solution(id)
                        else; Lmgvalue(lvol, ii) = zero
                        end if

                        ; ; id = Lmh(lvol, ii)
                        if (id /= 0) then; Lmhvalue(lvol, ii) = solution(id)
                        else; Lmhvalue(lvol, ii) = zero
                        end if
                    end if
                end do ! ii
            end if ! YESMatrixFree
        else ! NOTstellsym;

            ; ii = 1
            do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                if (id /= 0) then; Ate(lvol, ideriv, ii)%s(ll) = solution(id)
                else; Ate(lvol, ideriv, ii)%s(ll) = zero
                end if
                ; ; id = Aze(lvol, 0, ii)%i(ll); 
                if (id /= 0) then; Aze(lvol, ideriv, ii)%s(ll) = solution(id)
                else; Aze(lvol, ideriv, ii)%s(ll) = zero
                end if
                ; ; ; Ato(lvol, ideriv, ii)%s(ll) = zero         ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
                ; ; ; Azo(lvol, ideriv, ii)%s(ll) = zero         ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
            end do
            do ii = 2, mn
                do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Ate(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Ate(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Aze(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Aze(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Aze(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Ato(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Ato(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Ato(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Azo(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Azo(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Azo(lvol, ideriv, ii)%s(ll) = zero
                    end if
                end do ! end of do ll;
            end do ! end of do ii;

            if (YESMatrixFree) then
                do ii = 1, mn
                    ; ; id = Lma(lvol, ii)
                    if (id /= 0) then; Lmavalue(lvol, ii) = solution(id)
                    else; Lmavalue(lvol, ii) = zero
                    end if

                    ; ; id = Lmb(lvol, ii)
                    if (id /= 0) then; Lmbvalue(lvol, ii) = solution(id)
                    else; Lmbvalue(lvol, ii) = zero
                    end if

                    if (ii .gt. 1) then; id = Lmc(lvol, ii)
                        if (id /= 0) then; Lmcvalue(lvol, ii) = solution(id)
                        else; Lmcvalue(lvol, ii) = zero
                        end if

                        ; ; id = Lmd(lvol, ii)
                        if (id /= 0) then; Lmdvalue(lvol, ii) = solution(id)
                        else; Lmdvalue(lvol, ii) = zero
                        end if

                        ; ; id = Lme(lvol, ii)
                        if (id /= 0) then; Lmevalue(lvol, ii) = solution(id)
                        else; Lmevalue(lvol, ii) = zero
                        end if

                        ; ; id = Lmf(lvol, ii)
                        if (id /= 0) then; Lmfvalue(lvol, ii) = solution(id)
                        else; Lmfvalue(lvol, ii) = zero
                        end if

                    else; id = Lmg(lvol, ii)
                        if (id /= 0) then; Lmgvalue(lvol, ii) = solution(id)
                        else; Lmgvalue(lvol, ii) = zero
                        end if

                        ; ; id = Lmh(lvol, ii)
                        if (id /= 0) then; Lmhvalue(lvol, ii) = solution(id)
                        else; Lmhvalue(lvol, ii) = zero
                        end if
                    end if
                end do ! ii
            end if ! YESMatrixFree
        end if ! end of if( YESstellsym );

    case ('P')

!   FATAL( packab, .true., a trivial revision of packab is required to remove the if from the loop )

        solution = zero

        do ii = 1, mn
            do ll = 0, llrad
                ; ; id = Ate(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Ate(lvol, ideriv, ii)%s(ll)
                ; ; id = Aze(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Aze(lvol, ideriv, ii)%s(ll)
                if (ii .gt. 1 .and. NOTstellsym) then; id = Ato(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Ato(lvol, ideriv, ii)%s(ll)
                    ; ; id = Azo(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Azo(lvol, ideriv, ii)%s(ll)
                end if
            end do ! end of do ll;
        end do ! end of do ii;

    end select ! end of select case( packorunpack );

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

    if (Wpackab) then

        if (ideriv .eq. 0) then

            do ii = 1, mn
                write (ounit, 1000) myid, lvol, im(ii), in(ii), "Ate", Ate(lvol, 0, ii)%s(0:llrad)
                write (ounit, 1000) myid, lvol, im(ii), in(ii), "Aze", Aze(lvol, 0, ii)%s(0:llrad)
                write (ounit, 1000) myid, lvol, im(ii), in(ii), "Ato", Ato(lvol, 0, ii)%s(0:llrad)
                write (ounit, 1000) myid, lvol, im(ii), in(ii), "Azo", Azo(lvol, 0, ii)%s(0:llrad)
            end do ! end of do ii;

            do ii = 1, mn
                write (ounit, 1000) myid, lvol, im(ii), in(ii), "Ste", sum(Ate(lvol, 0, ii)%s(0:llrad)*TT(0:llrad, 0, 0)), sum(Ate(lvol, 0, ii)%s(0:llrad)*TT(0:llrad, 1, 0))
                write (ounit, 1000) myid, lvol, im(ii), in(ii), "Sze", sum(Aze(lvol, 0, ii)%s(0:llrad)*TT(0:llrad, 0, 0)), sum(Aze(lvol, 0, ii)%s(0:llrad)*TT(0:llrad, 1, 0))
            end do ! end of do ii;

        end if ! end of if( ideriv.eq.0);

    end if ! end of if( Wpackab );

1000 format("packab : "10x " : myid="i3" ; lvol="i3" ; ("i3" ,"i3" ) : "a3"=["999(es23.15","))

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
    cput = MPI_WTIME()
    Tpackab = Tpackab + (cput - cpuo)
    return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine packab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
