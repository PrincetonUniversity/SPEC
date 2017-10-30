!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Constructs matrix that represents the vacuum linear system.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{enumerate} \subsubsection{representation of scalar potential} \begin{enumerate}

!latex \item The (stellarator-symmetric) scalar potential may be written
!latex       \be \phi = I \t + G \z + \sum_{p,j} \phi_{p,j} T_p \sin\a_j.
!latex       \ee
!latex       where the currents, $I$ and $G$, are assumed given and the $\phi_{p,j}$ are to be determined.

!latex \item An arbitrary ``basis-function'' is $\varphi \equiv T_l \sin \a_i$.

!latex \item From the identity $\int_{{\cal V}} \nabla \cdot ( \varphi \nabla \phi ) dv = \int_{\partial {\cal V}} \varphi \nabla \phi . d{\bf s}$, 
!latex       and using the ``weak condition'', $\int \varphi \nabla \cdot \nabla \phi \; dv = 0$, we obtain
!latex       \be \int_{\partial {\cal V}} \varphi \; \nabla \phi . d{\bf s} = \int_{{\cal V}} \nabla \varphi \cdot \nabla \phi \; dv.
!latex       \ee

!latex \item This routine will assume that $\nabla \phi . d{\bf s} \equiv b(\t,\z)\; d\t d\z$ on the computational boundary is given.

!latex \end{enumerate} \subsubsection{linear equations} \begin{enumerate}

!latex \item A system of linear equations is derived:
!latex       \be \begin{array}{ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc}
!latex            \ooint T_l \sin\a_i \; b
!latex            & = & I &            &     &     & \int \! ds & T_l^\prime &            & \ooint \sin\a_i &          & \bar g^{\s\t} \\
!latex            & + & G &            &     &     & \int \! ds & T_l^\prime &            & \ooint \sin\a_i &          & \bar g^{\s\z} \\
!latex            & + &   & \phi_{p,j} &     &     & \int \! ds & T_l^\prime & T_p^\prime & \ooint \sin\a_i & \sin\a_j & \bar g^{\s\s} \\
!latex            & + &   & \phi_{p,j} &     & m_j & \int \! ds & T_l^\prime & T_p        & \ooint \sin\a_i & \cos\a_j & \bar g^{\s\t} \\
!latex            & - &   & \phi_{p,j} &     & n_j & \int \! ds & T_l^\prime & T_p        & \ooint \sin\a_i & \cos\a_j & \bar g^{\s\z} \\
!latex            & + & I &            & m_i &     & \int \! ds & T_l        &            & \ooint \cos\a_i &          & \bar g^{\t\t} \\
!latex            & + & G &            & m_i &     & \int \! ds & T_l        &            & \ooint \cos\a_i &          & \bar g^{\t\z} \\
!latex            & + &   & \phi_{p,j} & m_i &     & \int \! ds & T_l        & T_p^\prime & \ooint \cos\a_i & \sin\a_j & \bar g^{\s\t} \\
!latex            & + &   & \phi_{p,j} & m_i & m_j & \int \! ds & T_l        & T_p        & \ooint \cos\a_i & \cos\a_j & \bar g^{\t\t} \\
!latex            & - &   & \phi_{p,j} & m_i & n_j & \int \! ds & T_l        & T_p        & \ooint \cos\a_i & \cos\a_j & \bar g^{\t\z} \\
!latex            & - & I &            & n_i &     & \int \! ds & T_l        &            & \ooint \cos\a_i &          & \bar g^{\z\t} \\
!latex            & - & G &            & n_i &     & \int \! ds & T_l        &            & \ooint \cos\a_i &          & \bar g^{\z\z} \\
!latex            & - &   & \phi_{p,j} & n_i &     & \int \! ds & T_l        & T_p^\prime & \ooint \cos\a_i & \sin\a_j & \bar g^{\s\z} \\
!latex            & - &   & \phi_{p,j} & n_i & m_j & \int \! ds & T_l        & T_p        & \ooint \cos\a_i & \cos\a_j & \bar g^{\t\z} \\
!latex            & + &   & \phi_{p,j} & n_i & n_j & \int \! ds & T_l        & T_p        & \ooint \cos\a_i & \cos\a_j & \bar g^{\z\z}
!latex       \end{array} \ee
!latex       where $\bar g^{uv} \equiv g^{uv} \sqrt g$ and summation over $p$ and $j$ is implied.

!latex \item The quantities $\verb+bns(i)+\equiv\ooint \sin\a_i \nabla\phi\cdot\nabla s \sqrt g$ are provided as input.

!latex \end{enumerate} \subsubsection{volume integrated metric information} \begin{enumerate}

!latex \item The required information, including the non-stellarator-symmetric terms, is 
!latex        \be 
!latex        \verb+TTee(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \cos\a_i \cos\a_j \\
!latex        \verb+TTeo(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \cos\a_i \sin\a_j \\
!latex        \verb+TToe(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \sin\a_i \cos\a_j \\
!latex        \verb+TToo(1,l,p,i,j)+ & \equiv & \int ds \; T_l^\prime \; T_p^\prime \; \ooint \bar g^{\s\s} \sin\a_i \sin\a_j
!latex        \ee
!latex        \be 
!latex        \verb+TTee(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \cos\a_i \cos\a_j \\
!latex        \verb+TTeo(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \cos\a_i \sin\a_j \\
!latex        \verb+TToe(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \sin\a_i \cos\a_j \\
!latex        \verb+TToo(2,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\t} \sin\a_i \sin\a_j
!latex        \ee
!latex        \be 
!latex        \verb+TTee(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \cos\a_i \cos\a_j \\
!latex        \verb+TTeo(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \cos\a_i \sin\a_j \\
!latex        \verb+TToe(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \sin\a_i \cos\a_j \\
!latex        \verb+TToo(3,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p^\prime \; \ooint \bar g^{\s\z} \sin\a_i \sin\a_j
!latex        \ee
!latex        \be 
!latex        \verb+TTee(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \cos\a_i \cos\a_j \\
!latex        \verb+TTeo(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \cos\a_i \sin\a_j \\
!latex        \verb+TToe(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \sin\a_i \cos\a_j \\
!latex        \verb+TToo(4,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\t} \sin\a_i \sin\a_j
!latex        \ee
!latex        \be 
!latex        \verb+TTee(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \cos\a_i \cos\a_j \\
!latex        \verb+TTeo(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \cos\a_i \sin\a_j \\
!latex        \verb+TToe(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \sin\a_i \cos\a_j \\
!latex        \verb+TToo(5,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\t\z} \sin\a_i \sin\a_j
!latex        \ee
!latex        \be 
!latex        \verb+TTee(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \cos\a_i \cos\a_j \\
!latex        \verb+TTeo(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \cos\a_i \sin\a_j \\
!latex        \verb+TToe(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \sin\a_i \cos\a_j \\
!latex        \verb+TToo(6,l,p,i,j)+ & \equiv & \int ds \; T_l        \; T_p        \; \ooint \bar g^{\z\z} \sin\a_i \sin\a_j
!latex        \ee

!!atex        \be 
!!atex        \verb+Te(1,l,i)+ & \equiv & \int ds \; T_l^\prime               \; \ooint \bar g^{\s\s} \cos\a_i          \\
!!atex        \verb+To(1,l,i)+ & \equiv & \int ds \; T_l^\prime               \; \ooint \bar g^{\s\s} \sin\a_i
!!atex        \ee
!latex        \be 
!latex        \verb+Te(2,l,i)+ & \equiv & \int ds \; T_l^\prime               \; \ooint \bar g^{\s\t} \cos\a_i          \\
!latex        \verb+To(2,l,i)+ & \equiv & \int ds \; T_l^\prime               \; \ooint \bar g^{\s\t} \sin\a_i
!latex        \ee
!latex        \be 
!latex        \verb+Te(3,l,i)+ & \equiv & \int ds \; T_l^\prime               \; \ooint \bar g^{\s\z} \cos\a_i          \\
!latex        \verb+To(3,l,i)+ & \equiv & \int ds \; T_l^\prime               \; \ooint \bar g^{\s\z} \sin\a_i
!latex        \ee
!latex        \be 
!latex        \verb+Te(4,l,i)+ & \equiv & \int ds \; T_l                      \; \ooint \bar g^{\s\t} \cos\a_i          \\
!latex        \verb+To(4,l,i)+ & \equiv & \int ds \; T_l                      \; \ooint \bar g^{\s\t} \sin\a_i
!latex        \ee
!latex        \be 
!latex        \verb+Te(5,l,i)+ & \equiv & \int ds \; T_l                      \; \ooint \bar g^{\t\z} \cos\a_i          \\
!latex        \verb+To(5,l,i)+ & \equiv & \int ds \; T_l                      \; \ooint \bar g^{\t\z} \sin\a_i
!latex        \ee
!latex        \be 
!latex        \verb+Te(6,l,i)+ & \equiv & \int ds \; T_l                      \; \ooint \bar g^{\z\z} \cos\a_i          \\
!latex        \verb+To(6,l,i)+ & \equiv & \int ds \; T_l                      \; \ooint \bar g^{\z\z} \sin\a_i
!latex        \ee

!latex \end{enumerate} \subsubsection{useful identities} \begin{enumerate}

!latex \item The following identities may be useful:
!latex       \be 
!latex       \sqrt g \; g^{\s\s} & = & \left( g_{\t\t} g_{\z\z} - g_{\t\z} g_{\t\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\s\t} & = & \left( g_{\t\z} g_{\s\z} - g_{\s\t} g_{\z\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\s\z} & = & \left( g_{\s\t} g_{\t\z} - g_{\t\t} g_{\s\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\t\t} & = & \left( g_{\z\z} g_{\s\s} - g_{\s\z} g_{\s\z} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\t\z} & = & \left( g_{\s\z} g_{\s\t} - g_{\t\z} g_{\s\s} \right) / \sqrt g \\
!latex       \sqrt g \; g^{\z\z} & = & \left( g_{\s\s} g_{\t\t} - g_{\s\t} g_{\s\t} \right) / \sqrt g
!latex       \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine va00aa( mn, lrad )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use constants, only : zero, half, quart, one

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wva00aa

  use cputiming, only : Tva00aa

  use allglobal, only : myid, ncpu, cpus, &
			pi2pi2nfphalf, &
                        Mvol, &
			NOTstellsym, &
			im, in, &
			Nmagneticdof, Ate, Ato, dMA, dMB, &
			TTee, TTeo, TToe, TToo, Te, To, &
			iBns, iBnc

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOCALS
				  
  INTEGER, intent(in)  :: mn, lrad

  INTEGER              :: NN, lvol, id, jd, ii, jj, ll, pp, mi, ni, mj, nj, mimj, minj, nimj, ninj, ideriv!, iuv

  BEGIN(va00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  NN = Nmagneticdof(Mvol) ; lvol=Mvol ; ideriv = 0 ! shorthand; Nmagneticdof was computed in al00aa;
  
  dMA(1:NN,1:NN) = zero ! these "zeros" are probably not required
  dMB(1:NN,0: 2) = zero

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do ii = 2, mn ; mi = im(ii) ; ni = in(ii)
   
   do ll = 0, lrad
    
    do jj = 2, mn ; mj = im(jj) ; nj = in(jj)
     
     mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj
     
     do pp = 0, lrad
      
      
      id=Ate(lvol,ideriv,ii)%i(ll) ; jd=Ate(lvol,ideriv,jj)%i(pp) ; dMA(id,jd) =      TToo(1,ll,pp,ii,jj) + mimj * TTee(4,ll,pp,ii,jj) + ninj * TTee(6,ll,pp,ii,jj) &
                                                                               + mi * TTeo(2,ll,pp,ii,jj) + mj   * TToe(2,pp,ll,ii,jj) - ni   * TTeo(3,ll,pp,ii,jj) &
                                                                               - nj * TToe(3,pp,ll,ii,jj) - minj * TTee(5,ll,pp,ii,jj) - nimj * TTee(5,ll,pp,ii,jj)
      
      if( NOTstellsym ) then
      FATALMESS(va00aa, .true., non-symmetric terms under construction )
!! DO NOT DELETE THE FOLLOWING: 03 Apr 13;
!     MA(ll,ii,pp,jj,2) =      TToe(1,ll,pp,ii,jj) - mimj * TTeo(4,ll,pp,ii,jj) - ninj * TTeo(6,ll,pp,ii,jj) & !! DO NOT DELETE; 03 Apr 13;
!                       + mi * TTee(2,ll,pp,ii,jj) - mj   * TToo(2,pp,ll,ii,jj) - ni   * TTee(3,ll,pp,ii,jj) &
!                       + nj * TToo(3,pp,ll,ii,jj) + minj * TTeo(5,ll,pp,ii,jj) + nimj * TTeo(5,ll,pp,ii,jj)
!     MA(ll,ii,pp,jj,3) =      TTeo(1,ll,pp,ii,jj) - mimj * TToe(4,ll,pp,ii,jj) - ninj * TToe(6,ll,pp,ii,jj) & !! DO NOT DELETE; 03 Apr 13;
!                       - mi * TToo(2,ll,pp,ii,jj) + mj   * TTee(2,pp,ll,ii,jj) + ni   * TToo(3,ll,pp,ii,jj) &
!                       - ni * TTee(2,pp,ll,ii,jj) + minj * TToe(5,ll,pp,ii,jj) + nimj * TToe(5,ll,pp,ii,jj)
!     MA(ll,ii,pp,jj,4) =      TTee(1,ll,pp,ii,jj) + mimj * TToo(4,ll,pp,ii,jj) + ninj * TToo(6,ll,pp,ii,jj) & !! DO NOT DELETE; 03 Apr 13;
!                       - mi * TToe(2,ll,pp,ii,jj) - mj   * TTeo(2,pp,ll,ii,jj) + ni   * TToe(3,ll,pp,ii,jj) &
!                       + nj * TTeo(3,pp,ll,ii,jj) - minj * TToo(5,ll,pp,ii,jj) - nimj * TToo(5,ll,pp,ii,jj)
!! DO NOT DELETE THE ABOVE;     03 Apr 13;
      endif ! end of if( NOTstellsym ) 
      
     enddo ! end of do pp;

    enddo ! end of do jj !! still inside do ii and do ll; 03 Apr 13;

   !id=Ate(lvol,ideriv,ii)%i(ll) ; jd=0 ; dMB(id,jd) = iBns(ii)                 !! was this the last bug; 17 Apr 13;
    id=Ate(lvol,ideriv,ii)%i(ll) ; jd=0 ; dMB(id,jd) = iBns(ii) * pi2pi2nfphalf !! 17 Apr 13;

    id=Ate(lvol,ideriv,ii)%i(ll) ; jd=1 ; dMB(id,jd) = To(2,ll,ii) + mi * Te(4,ll,ii) - ni * Te(5,ll,ii)
    id=Ate(lvol,ideriv,ii)%i(ll) ; jd=2 ; dMB(id,jd) = To(3,ll,ii) + mi * Te(5,ll,ii) - ni * Te(6,ll,ii)

    if( NOTstellsym ) then
    FATALMESS(va00aa, .true., non-symmetric terms under construction )
!   MB(ll,ii,1,0) = MB(ll,ii,1,0) + b(jj,2) * Teo(ll,ii,jj)
!   MB(ll,ii,2,0) = MB(ll,ii,2,0) + b(jj,1) * Toe(ll,ii,jj) + b(jj,2)*Too(ll,ii,jj)
!   MB(ll,ii,2,1) = Itor * ( To(2,ll,ii) - mi * Te(4,ll,ii) + ni * Te(5,ll,ii) )
!   MB(ll,ii,2,2) = Ipol * ( To(3,ll,ii) - mi * Te(5,ll,ii) + ni * Te(6,ll,ii) )
    endif ! end of if( NOTstellsym )
    
   enddo ! end of do ll; 
   
  enddo ! end of do ii;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RETURN(va00aa)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine va00aa

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
