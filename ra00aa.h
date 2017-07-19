!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (output) ! Writes vector potential to .ext.sp.A .

!latex \briefly{Writes vector potential to .ext.sp.A .}

!latex \calledby{\link{preset} and \link{xspech}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{representation of vector potential} \begin{enumerate}

!latex \item \Ais

!latex \end{enumerate} \subsection{file format} \begin{enumerate}

!latex \item The format of the files containing the vector potential is as follows:

!latex \begin{verbatim}
!latex open(aunit, file="."//trim(ext)//".sp.A", status="replace", form="unformatted" )
!latex write(aunit) Mvol, Mpol, Ntor, mn, Nfp ! integers;
!latex write(aunit) im(1:mn) ! integers; poloidal modes;
!latex write(aunit) in(1:mn) ! integers; toroidal modes;
!latex do vvol = 1, Mvol ! integers; loop over volumes;
!latex write(aunit) Lrad(vvol) ! integers; the radial resolution in each volume may be different;
!latex do ii = 1, mn
!latex write(aunit) Ate(vvol,ii)%s(0:Lrad(vvol)) ! reals;
!latex write(aunit) Aze(vvol,ii)%s(0:Lrad(vvol)) ! reals;
!latex write(aunit) Ato(vvol,ii)%s(0:Lrad(vvol)) ! reals;
!latex write(aunit) Azo(vvol,ii)%s(0:Lrad(vvol)) ! reals;
!latex enddo ! end of do ii;
!latex enddo ! end of do vvol;
!latex close(aunit)
!latex \end{verbatim}

!latex \end{enumerate}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ra00aa( writeorread ) 
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero
  
  use numerical, only : 
  
  use fileunits, only : ounit, aunit
  
  use inputlist, only : Wmacros, Wra00aa, ext, Nfp, Mpol, Ntor, Lrad
  
  use cputiming, only : Tra00aa
  
  use allglobal, only : myid, ncpu, cpus, Mvol, mn, im, in, Ate, Aze, Ato, Azo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOCALS
  
  CHARACTER, intent(in) :: writeorread
  
  LOGICAL               :: exist                                
  
  INTEGER               :: vvol, oldMvol, oldMpol, oldNtor, oldmn, oldNfp, oldLrad, ii, jj, minLrad, llmodnp, ideriv
  INTEGER, allocatable  :: oldim(:), oldin(:)
  REAL   , allocatable  :: oldAte(:), oldAze(:), oldAto(:), oldAzo(:)
  
  BEGIN(ra00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ideriv = 0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( writeorread )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 'W' ) ! write vector potential harmonics to file;
   
   if( myid.eq.0 ) then
    
    open(aunit, file="."//trim(ext)//".sp.A", status="replace", form="unformatted" )
    
    write(aunit) Mvol, Mpol, Ntor, mn, Nfp
    write(aunit) im(1:mn)
    write(aunit) in(1:mn)
    
    do vvol = 1, Mvol
     
     write(aunit) Lrad(vvol) ! radial resolution depends on volume; 26 Feb 13;
     
     do ii = 1, mn ! loop over Fourier harmonics;
      
#ifdef DEBUG
      FATAL( ra00aa, .not.allocated(Ate(vvol,ideriv,ii)%s), error )
      FATAL( ra00aa, .not.allocated(Aze(vvol,ideriv,ii)%s), error )
      FATAL( ra00aa, .not.allocated(Ato(vvol,ideriv,ii)%s), error )
      FATAL( ra00aa, .not.allocated(Azo(vvol,ideriv,ii)%s), error )
      FATAL( ra00aa, Lrad(vvol).le.0, error )
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
   
   FATAL( ra00aa, .true., under reconstruction )

   if( myid.eq.0 ) then
    
    inquire(file="."//trim(ext)//".sp.A",exist=exist)   
    
    if( .not.exist ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; .ext.sp.A does not exist ;")') cput-cpus, myid ; goto 9998
    endif
    
    open(aunit,file="."//trim(ext)//".sp.A",status="old",form="unformatted",iostat=ios) ! this will contain initial guess for vector potential;
    
    if( ios.ne.0 ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; opening .ext.sp.A ;")') cput-cpus, myid ; goto 9997
    endif
    
    read(aunit,iostat=ios) oldMvol, oldMpol, oldNtor, oldmn, oldNfp ! these are the "old" resolution parameters;
    
    if( ios.ne.0 ) then
     write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; reading oldMvol, oldMpol, oldNtor, oldmn, oldNfp;")') cput-cpus, myid
     goto 9997
    endif
    
    if( oldNfp .ne.Nfp  ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; inconsistent Nfp ; ")') cput-cpus, myid ; goto 9997
    endif
    if( oldMvol.ne.Mvol ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; inconsistent Mvol ;")') cput-cpus, myid ; goto 9997
    endif
    
    SALLOCATE( oldim, (1:oldmn), 0 )
    SALLOCATE( oldin, (1:oldmn), 0 )
    
    read(aunit,iostat=ios) oldim(1:oldmn)
    read(aunit,iostat=ios) oldin(1:oldmn)
    
    do vvol = 1, oldMvol
     
     read(aunit,iostat=ios) oldLrad
     
     minLrad = min(oldLrad,Lrad(vvol))
     
     SALLOCATE( oldAte, (0:oldLrad), zero )
     SALLOCATE( oldAze, (0:oldLrad), zero )
     SALLOCATE( oldAto, (0:oldLrad), zero )
     SALLOCATE( oldAzo, (0:oldLrad), zero )
     
     do jj = 1, oldmn
      
      read(aunit,iostat=ios) oldAte(0:oldLrad)
      read(aunit,iostat=ios) oldAze(0:oldLrad)
      read(aunit,iostat=ios) oldAto(0:oldLrad)
      read(aunit,iostat=ios) oldAzo(0:oldLrad)
      
      do ii = 1, mn ! compare Fourier harmonic with old; 26 Feb 13;
       if( im(ii).eq.oldim(jj) .and. in(ii).eq.oldin(jj) ) then ; Ate(vvol,ideriv,ii)%s(0:minLrad) = oldAte(0:minLrad)
        ;                                                       ; Aze(vvol,ideriv,ii)%s(0:minLrad) = oldAze(0:minLrad)
        ;                                                       ; Ato(vvol,ideriv,ii)%s(0:minLrad) = oldAto(0:minLrad)
        ;                                                       ; Azo(vvol,ideriv,ii)%s(0:minLrad) = oldAzo(0:minLrad)
       endif
      enddo ! end of do ii; 26 Feb 13;
      
     enddo ! end of do jj; 26 Feb 13;
     
     DALLOCATE(oldAte)
     DALLOCATE(oldAze)
     DALLOCATE(oldAto)
     DALLOCATE(oldAzo)
     
    enddo ! end of do vvol; 26 Feb 13;
    
    DALLOCATE(oldim)
    DALLOCATE(oldin)
    
9997 continue
    
    close(aunit)
    
9998 continue
    
   endif  ! end of if( myid.eq.0 ) ; 26 Feb 13;
   
   do vvol = 1, Mvol
    
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  case default
   
   FATAL(ra00aa, .true., invalid writeorread flag supplied on input )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(ra00aa)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine ra00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


