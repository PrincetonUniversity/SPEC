!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Reads/writes vector potential from .AtAz.xxxx files.

!latex \end{enumerate} \subsubsection{representation of vector potential} \begin{enumerate}

!latex \item The magnetic vector potential is 
!latex       \be {\bf A} = A_\t \nabla \t + A_\z \nabla \z,
!latex       \ee
!latex       where
!latex       \be A_\t(\s,\t,\z) &=& \sum_{j,l} A_{\t,e,j,l} \; T_l(s) \cos(m_j\t-n_j\z)+ \sum_{j,l} A_{\t,o,j,l} \; T_l(s) \sin(m_j\t-n_j\z), \\
!latex           A_\z(\s,\t,\z) &=& \sum_{j,l} A_{\z,e,j,l} \; T_l(s) \cos(m_j\t-n_j\z)+ \sum_{j,l} A_{\z,o,j,l} \; T_l(s) \sin(m_j\t-n_j\z).
!latex       \ee

!latex \item Note that in each volume $\s\in[-1,1]$, and $T_l(s)$ are the Chebyshev polynomials.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ra00aa( writeorread ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, one, pi2
  
  use numerical, only : vsmall
  
  use fileunits, only : ounit, aunit
  
  use inputlist, only : Wmacros, Wra00aa, ext, Mpol, Ntor, Nvol, Lrad
  
  use cputiming, only : Tra00aa
  
  use allglobal, only : myid, ncpu, cpus, mn, im, in, Ate, Aze, Ato, Azo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  CHARACTER, intent(in) :: writeorread
  
  LOGICAL               :: exist                                
  
  INTEGER               :: vvol, oldNvol, ii, jj, oldmn, oldlrad, minlrad, llmodnp, ideriv
  INTEGER, allocatable  :: oldim(:), oldin(:)
  REAL   , allocatable  :: oldAte(:), oldAze(:), oldAto(:), oldAzo(:)
  
  BEGIN(ra00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ideriv = 0

#ifdef DEBUG
  if( Wra00aa ) then
  write(ounit,'("ra00aa : ", 10x ," : myid=",i3," ; writeorread="a2" ;")') myid, writeorread
  endif
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( writeorread )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \end{enumerate} \subsubsection{file format} \begin{enumerate}

!latex \item The format of the files containing the vector potential is as follows: \begin{enumerate}
!latex open(aunit, file="."//trim(ext)//".AtAzmn", status="replace", form="unformatted" )
!latex 
!latex write(aunit) Nvol, mn ! integers;
!latex 
!latex write(aunit) im(1:mn) ! integers;
!latex 
!latex write(aunit) in(1:mn) ! integers;
!latex 
!latex do vvol = 1, Nvol ! integers;
!latex 
!latex write(aunit) Lrad(vvol) ! integers; the radial resolution in each volume may be different;
!latex 
!latex do ii = 1, mn
!latex 
!latex write(aunit) Ate(vvol,ii)\%s(0:Lrad(vvol)) ! reals;
!latex 
!latex write(aunit) Aze(vvol,ii)\%s(0:Lrad(vvol)) ! reals;
!latex 
!latex write(aunit) Ato(vvol,ii)\%s(0:Lrad(vvol)) ! reals;
!latex 
!latex write(aunit) Azo(vvol,ii)\%s(0:Lrad(vvol)) ! reals;
!latex 
!latex enddo
!latex 
!latex enddo
!latex 
!latex close(aunit)
!latex 
!latex       \end{enumerate}

  case( 'W' ) ! write vector potential harmonics to file;
   
   if( myid.eq.0 ) then
    
    open(aunit, file="."//trim(ext)//".AtAzmn", status="replace", form="unformatted" )
    
    write(aunit) Nvol, mn
    
    write(aunit) im(1:mn)
    write(aunit) in(1:mn)
    
    do vvol = 1, Nvol ! shouldn't this be 1, Mvol? ; 15 Sep 15;
     
     write(aunit) Lrad(vvol) ! radial resolution depends on volume; 26 Feb 13;
     
     do ii = 1, mn ! loop over Fourier harmonics;
      
#ifdef DEBUG
      FATALMESS(ra00aa, .not.allocated(Ate(vvol,ideriv,ii)%s), error)
      FATALMESS(ra00aa, .not.allocated(Aze(vvol,ideriv,ii)%s), error)
      FATALMESS(ra00aa, .not.allocated(Ato(vvol,ideriv,ii)%s), error)
      FATALMESS(ra00aa, .not.allocated(Azo(vvol,ideriv,ii)%s), error)
      FATALMESS(ra00aa, Lrad(vvol).le.0, error)
#endif

      write(aunit) Ate(vvol,ideriv,ii)%s(0:Lrad(vvol))
      write(aunit) Aze(vvol,ideriv,ii)%s(0:Lrad(vvol))
      write(aunit) Ato(vvol,ideriv,ii)%s(0:Lrad(vvol))
      write(aunit) Azo(vvol,ideriv,ii)%s(0:Lrad(vvol))
      
     enddo ! end of do ii;  6 Feb 13;
     
    enddo ! end of do vvol;  6 Feb 13;
    
    close(aunit)
    
   endif ! end of if( myid.eq.0 ) ;  6 Feb 13;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 'R' ) ! read potential from file; interpolate onto new radial grid;
   
   if( myid.eq.0 ) then
    
    inquire(file=".AtAzmn",exist=exist)   
    
    if( .not.exist ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; .ext.AtAzmn does not exist ;")') cput-cpus, myid ; goto 9998
    endif
    
    open(aunit,file=".AtAzmn",status="old",form="unformatted",iostat=ios) ! this will contain initial guess for vector potential;
    
    if( ios.ne.0 ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; opening .ext.AtAzmn ;")') cput-cpus, myid ; goto 9997
    endif
    
    read(aunit,iostat=ios) oldNvol, oldmn ! these are the "old" resolution parameters;
    
    if( ios.ne.0 ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; reading oldNvol, oldmn ;")') cput-cpus, myid ; goto 9997
    endif
    
    IALLOCATE(oldim,(1:oldmn))
    IALLOCATE(oldin,(1:oldmn))
    
    read(aunit,iostat=ios) oldim(1:oldmn)
    read(aunit,iostat=ios) oldin(1:oldmn)
    
    do vvol = 1, oldNvol
     
     read(aunit,iostat=ios) oldlrad
     
     minlrad = min(oldlrad,Lrad(vvol))
     
     RALLOCATE(oldAte,(0:oldlrad))
     RALLOCATE(oldAze,(0:oldlrad))
     RALLOCATE(oldAto,(0:oldlrad))
     RALLOCATE(oldAzo,(0:oldlrad))
     
     do jj = 1, oldmn
      
      read(aunit,iostat=ios) oldAte(0:oldlrad)
      read(aunit,iostat=ios) oldAze(0:oldlrad)
      read(aunit,iostat=ios) oldAto(0:oldlrad)
      read(aunit,iostat=ios) oldAzo(0:oldlrad)
      
      do ii = 1, mn ! compare Fourier harmonic with old; 26 Feb 13;
       if( im(ii).eq.oldim(jj) .and. in(ii).eq.oldin(jj) ) then ; Ate(vvol,ideriv,ii)%s(0:minlrad) = oldAte(0:minlrad)
        ;                                                       ; Aze(vvol,ideriv,ii)%s(0:minlrad) = oldAze(0:minlrad)
        ;                                                       ; Ato(vvol,ideriv,ii)%s(0:minlrad) = oldAto(0:minlrad)
        ;                                                       ; Azo(vvol,ideriv,ii)%s(0:minlrad) = oldAzo(0:minlrad)
       endif
      enddo ! end of do ii; 26 Feb 13;
      
     enddo ! end of do jj; 26 Feb 13;
     
     DEALLOCATE(oldAte)
     DEALLOCATE(oldAze)
     DEALLOCATE(oldAto)
     DEALLOCATE(oldAzo)
     
    enddo ! end of do vvol; 26 Feb 13;
    
    DEALLOCATE(oldim)
    DEALLOCATE(oldin)
    
9997 continue
    
    close(aunit)
    
9998 continue
    
   endif  ! end of if( myid.eq.0 ) ; 26 Feb 13;
   
   do vvol = 1, Nvol
    
    llmodnp = 0 ! this node contains the information that is to be broadcast; 26 Feb 13;
    
    do ii = 1, mn 
     RlBCAST( Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
     RlBCAST( Aze(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
    enddo
   !if( NOTstellsym ) then
    do ii = 1, mn 
     RlBCAST( Ato(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
     RlBCAST( Azo(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
    enddo
   !endif
    
   enddo ! end of do vvol; 26 Feb 13;

#ifdef DEBUG
   if( Wra00aa ) then
    do vvol = 1, Nvol
     do ii = 1, 1
      write(ounit,'("ra00aa : ", 10x ," : myid=",i3," ; vvol=",i3," ; (",i3," ," ,i3," ) : Ate="99f9.4)') &
   myid, vvol, im(ii), in(ii), Ate(vvol,ideriv,ii)%s(0:Lrad(vvol))
     enddo
    enddo
   endif
#endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  case default
   
   FATALMESS(ra00aa, .true., invalid writeorread flag supplied on input )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(ra00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ra00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
