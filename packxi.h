!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (&ldquo;packing&rdquo;) ! Packs, and unpacks, geometrical degrees of freedom; ${\bf \xi} \leftrightarrow \{ R_{i,l},Z_{i,l} \}$.

!latex \briefly{Packs, and unpacks, geometrical degrees of freedom; and sets coordinate axis.}

!latex \calledby{\link{dforce}, \link{global}, \link{hesian}, \link{newton} and \link{xspech}}
!latex \calls{\link{rzaxis}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{geometrical degrees of freedom}

!latex \begin{enumerate}
!latex \item The geometrical degrees-of-freedom,
!latex       namely the $R_{j,v}$ and $Z_{j,v}$ where $v$ labels the interface and $j$ labels the Fourier harmonic,
!latex       must be ``packxi'', and ``unpackxi'', into a single vector, $\boldxi$, so that standard numerical routines can be
!latex       called to find solutions to force-balance, i.e. ${\bf F}[\boldxi]=0$.
!latex \item A coordinate ``pre-conditioning'' factor is included: 
!latex       \be \boldxi_k \equiv \frac{R_{j,v}}{\Psi_{j,v}},
!latex       \ee
!latex       where $\Psi_{j,v} \equiv $ \type{psifactor(j,v)}, which is defined in \link{global}.
!latex \end{enumerate} 

!latex \subsection{coordinate axis}

!latex \begin{enumerate}
!latex \item The coordinate axis is not an independent degree-of-freedom of the geometry.
!latex       It is constructed by extrapolating the geometry of the innermost interface down to a line.
!latex \item Note that if the coordinate axis depends only on the geometry of the innermost interface
!latex       then the block tridiagonal structure of the the force-derivative matrix is preserved.
!l tex \item If required, the coordinate axis, $R_0(\z)$ and $Z_0(\z)$, is set as the ``midpoint'' of the innermost interface, $R_1(\t,\z)$ and $Z_1(\t,\z)$,
!l tex       according to
!l tex       \be R_0(\z) \equiv \frac{1}{2}\left[ R_1(0,\z) + R_1(\pi,\z)  \right], \qquad
!l tex           Z_0(\z) \equiv \frac{1}{2}\left[ Z_1(0,\z) + Z_1(\pi,\z)  \right].
!l tex       \ee
!l tex \item This becomes
!l tex       \be R_{j,0} = \sum_{i} R_{i,1} |\beta_{i,j}|, \qquad
!l tex           Z_{j,0} = \sum_{i} Z_{i,1}  \beta_{i,j},
!l tex       \ee
!l tex       where the $\beta \equiv $ \type{bjk} are set in \link{global}.
!l tex \item Note that this construction of the coordinate axis depends on the angle, $\t$,
!l tex       so that if a different angle parameterization of the innermost interface is used then the geometry of the coordinate axis changes.
!l tex \end{enumerate} 
!latex \item Define the arc-length weighted averages,
!latex       \be R_0(\z) \equiv \frac{\int_{0}^{2\pi} R_1(\t,\z) dl}{L(\z)}, \qquad
!latex           Z_0(\z) \equiv \frac{\int_{0}^{2\pi} Z_1(\t,\z) dl}{L(\z)},
!latex       \ee
!latex       where $L(\z)\equiv \int_{0}^{2\pi} dl$ and $dl \equiv \sqrt{ \partial_\t R_1(\t,\z)^2 + \partial_\t Z_1(\t,\z)^2 } \, d\t$.
!latex \item Note that if $dl$ does not depend on $\t$, i.e. if $\t$ is the equal arc-length angle, then the expressions simplify.
!latex \item Note that the geometry of the coordinate axis thus constructed only depends on the geometry of the innermost interface, by which I 
!latex       mean that the geometry of the coordinate axis is independent of the angle parameterization.
!latex \end{enumerate}

!latex \subsection{some numerical comments}

!latex \begin{enumerate}
!latex \item First, the differential poloidal length, $dl \equiv \sqrt{ R_\t^2 + Z_\t^2 }$, is computed in real space using 
!latex       an inverse FFT the from Fourier harmonics of $R$ and $Z$.
!latex \item Second, the Fourier harmonics of the $dl$ are computed using an FFT.
!latex       The integration over $\t$ to construct $L\equiv \int dl$ is now trivial: just multiply the $m=0$ harmonics of $dl$ by $2\pi$.
!latex       The \internal{ajk(1:mn)} variable is used.
!latex \item Next, the weighted $R dl$ and $Z dl$ are computed in real space, and the poloidal integral is similarly taken.
!latex \item Lastly, the Fourier harmonics are constructed using an FFT after dividing in real space.
!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine packxi( NGdof, position, Mvol, mn, iRbc, iZbs, iRbs, iZbc, packorunpack )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  
  use numerical, only :
  
  use fileunits, only : ounit
  
  use inputlist, only : Wpackxi, Igeometry, Ntor, Nvol, Lfindzero
  
  use cputiming, only : Tpackxi
  
  use allglobal, only : ncpu, myid, cpus, im, in, &
                        YESstellsym, NOTstellsym, &
                        ajk, Nt, Nz, Ntz, iRij, iZij, tRij, tZij, &
                        ijreal, ijimag, jireal, jiimag, efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        psifactor, Rscale 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  INTEGER, intent(in)    :: NGdof, Mvol, mn
  REAL                   :: position(0:NGdof), iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol)
  CHARACTER              :: packorunpack
  
  INTEGER                :: lvol, jj, kk, irz, issym, idof, ifail, ivol
  
  BEGIN(packxi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  idof = 0 ! initialize counter; 14 Jan 13;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do lvol = 1, Mvol-1 ! loop over internal interfaces;
   
   do jj = 1, mn ! loop over Fourier harmonics;
    
    do irz = 0, 1 ! loop over R & Z;
     
     if( Igeometry.lt.3 .and. irz.eq.1 ) cycle ! no dependence on Z; 14 Jan 13;
     
     do issym = 0, 1 ! loop over even & odd;
      
      if( YESstellsym .and. issym.eq.1 ) cycle
      
      if( issym.eq.0 .and. irz.eq.1 .and. jj.eq.1 ) cycle ! no dependence on Zbs_{0,0}; 14 Jan 13;
      if( issym.eq.1 .and. irz.eq.0 .and. jj.eq.1 ) cycle ! no dependence on Rbs_{0,0}; 14 Jan 13;
      
      idof = idof + 1
      
#ifdef DEBUG
      FATAL( packxi, idof.le.0 .or. idof.gt.NGdof, out of bounds )
#endif
      
      select case( packorunpack )
       
      case( 'P' ) !   pack vector of unknowns;
       
       if( irz.eq.0 .and. issym.eq.0 ) position(idof) = iRbc(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) position(idof) = iZbs(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) position(idof) = iRbs(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) position(idof) = iZbc(jj,lvol) / psifactor(jj,lvol)
       
      case( 'U' ) ! unpack vector of unknowns;
       
       if( irz.eq.0 .and. issym.eq.0 ) iRbc(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) iZbs(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) iRbs(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) iZbc(jj,lvol) = position(idof) * psifactor(jj,lvol)
       
      end select
      
     enddo ! end of do issym;
     
    enddo ! end of do irz;
    
   enddo ! end of do jj;
   
  enddo ! end of do lvol;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
#ifdef DEBUG
  FATAL( packxi, idof.ne.NGdof, counting error )
#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( YESstellsym ) then ! iRbc(    ,0:Mvol) = zero
   ;                     ; iZbs(1   ,0:Mvol) = zero
   ;                     ; iRbs(1:mn,0:Mvol) = zero
   ;                     ; iZbc(1:mn,0:Mvol) = zero
  else                   ! iRbc(    ,0:Mvol) = zero
   ;                     ; iZbs(1   ,0:Mvol) = zero
   ;                     ; iRbs(1   ,0:Mvol) = zero
   ;                     ! iZbc(    ,0:Mvol) = zero
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  select case( packorunpack )

  case( 'P' )

  case( 'U' )

   ivol = 1 ! take care with ivol: this variable name might be a global variable, but here it is local; 19 Jul 16; 
 
   if( (Nvol .ne. 1) .and. (Lfindzero .ne. 0) ) then  
    WCALL( packxi, rzaxis, ( Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), ivol ) ) ! set coordinate axis; 19 Jul 16; 
   endif

  end select
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(packxi)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine packxi

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
