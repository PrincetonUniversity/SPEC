!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (&ldquo;packing&rdquo;) ! Packs, and unpacks, Beltrami field solution vector; ${\bf a}\equiv\{A_{\theta,e,i,l}, A_{\zeta,e,i,l}, etc.\}$.

!latex \briefly{Packs and unpacks Beltrami field solution vector.}

!latex \calledby{\link{dforce}, \link{ma02aa} and \link{mp00ac}}

!latex \tableofcontents

!latex \subsection{construction of ``vector'' of independent degrees of freedom}

!latex \begin{enumerate}
!latex \item Numerical routines for solving linear equations typically require the unknown, independent degrees of freedom
!latex       to be ``packed'' into a vector, ${\bf x}$.
!latex \item The magnetic field is defined by the independent degrees of freedom in
!latex       the Chebyshev-Fourier representation of the vector potential, $\Ate{i,l}$ and $\Aze{i,l}$;
!latex       and the non-stellarator-symmetric terms if relevant, $\Ato{i,l}$ and $\Azo{i,l}$;
!latex       and the Lagrange multipliers, $a_i$, $b_i$, $c_i$, $d_i$, $e_i$, etc. as required to enforce the constraints:
!latex       \be {\bf x} \equiv \{ \Ate{i,l},\Aze{i,l},\Ato{i,l},\Azo{i,l},a_i,b_i,c_i,d_i,e_i,f_i,g_1,h_1\}.
!latex       \ee
!latex \item The ``packing'' index is assigned in \link{preset}.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine packab( packorunpack, lvol, NN, solution, ideriv )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero
  
  use numerical, only : small
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpackab, Lrad
  
  use cputiming, only : Tpackab

  use allglobal, only : myid, ncpu, cpus, mn, im, in, Ate, Aze, Ato, Azo, YESstellsym, NOTstellsym, &
                        TT

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
  
  CHARACTER, intent(in) :: packorunpack
  INTEGER  , intent(in) :: lvol, NN, ideriv
  REAL                  :: solution(1:NN)
  
  INTEGER               :: ii, ll, id, llrad
    
  BEGIN(packab)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  llrad = Lrad(lvol) ! shorthand;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( packorunpack ) 
   
  case( 'U' )
   
#ifdef DEBUG
   
   if( Wpackab ) then
    
    if( YESstellsym ) then
     
     ;  ii = 1
      do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
       ;               ; id = Aze(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
      enddo ! end of do ll;
     do ii = 2, mn
      do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
       ;               ; id = Aze(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
      enddo ! end of do ll;
     enddo ! end of do ii;
      
    else ! NOTstellsym;
      
     ;  ii = 1
      do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
       ;               ; id = Aze(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
      enddo
     do ii = 2, mn
      do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
       ;               ; id = Aze(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
       ;               ; id = Ato(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
       ;               ; id = Azo(lvol,0,ii)%i(ll) ; FATAL( packab, id.lt.1 .or. id.gt.NN, unpacking illegal subscript )
      enddo ! end of do ll;
     enddo ! end of do ii;
     
    endif ! end ofif( YESstellsym );
     
   endif

#endif
    
   if( YESstellsym ) then
    
    do ii = 1, mn
     do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ate(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ate(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Aze(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Aze(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Aze(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ;                           ; Ato(lvol,ideriv,ii)%s(ll) = zero
      ;               ;                           ; Azo(lvol,ideriv,ii)%s(ll) = zero
     enddo ! end of do ll;
    enddo ! end of do ii;
    
   else ! NOTstellsym;
    
    ;  ii = 1
     do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ate(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ate(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Aze(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Aze(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Aze(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ;                           ; Ato(lvol,ideriv,ii)%s(ll) = zero         ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
      ;               ;                           ; Azo(lvol,ideriv,ii)%s(ll) = zero         ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
     enddo
    do ii = 2, mn
     do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ate(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ate(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Aze(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Aze(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Aze(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Ato(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ato(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ato(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Azo(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Azo(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Azo(lvol,ideriv,ii)%s(ll) = zero
       endif
     enddo ! end of do ll;
    enddo ! end of do ii;

   endif ! end of if( YESstellsym );
   
  case( 'P' )
   
!   FATAL( packab, .true., a trivial revision of packab is required to remove the if from the loop )
   
   solution = zero

   do ii = 1, mn
    do ll = 0, llrad
     ;                                    ; id = Ate(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Ate(lvol,ideriv,ii)%s(ll)
     ;                                    ; id = Aze(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Aze(lvol,ideriv,ii)%s(ll)
     if( ii.gt.1 .and. NOTstellsym ) then ; id = Ato(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Ato(lvol,ideriv,ii)%s(ll)
      ;                                   ; id = Azo(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Azo(lvol,ideriv,ii)%s(ll)
     endif
    enddo ! end of do ll;
   enddo ! end of do ii;
   
  end select ! end of select case( packorunpack );
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

  if( Wpackab ) then

   if( ideriv.eq.0 ) then

    do ii = 1, mn
     write(ounit,1000) myid, lvol, im(ii), in(ii), "Ate", Ate(lvol,0,ii)%s(0:llrad)
     write(ounit,1000) myid, lvol, im(ii), in(ii), "Aze", Aze(lvol,0,ii)%s(0:llrad)
     write(ounit,1000) myid, lvol, im(ii), in(ii), "Ato", Ato(lvol,0,ii)%s(0:llrad)
     write(ounit,1000) myid, lvol, im(ii), in(ii), "Azo", Azo(lvol,0,ii)%s(0:llrad)
    enddo ! end of do ii; 
    
    do ii = 1, mn
     write(ounit,1000) myid, lvol, im(ii), in(ii), "Ste", sum(Ate(lvol,0,ii)%s(0:llrad)*TT(0:llrad,0,0)), sum(Ate(lvol,0,ii)%s(0:llrad)*TT(0:llrad,1,0))
     write(ounit,1000) myid, lvol, im(ii), in(ii), "Sze", sum(Aze(lvol,0,ii)%s(0:llrad)*TT(0:llrad,0,0)), sum(Aze(lvol,0,ii)%s(0:llrad)*TT(0:llrad,1,0))
    enddo ! end of do ii;
        
   endif ! end of if( ideriv.eq.0);

  endif ! end of if( Wpackab );

1000 format("packab : " 10x " : myid="i3" ; lvol="i3" ; ("i3" ,"i3" ) : "a3"=["999(es23.15","))

#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(packab)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine packab

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
