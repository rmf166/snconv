    program main

      implicit none

      integer(8)                :: sol
      integer(8)                :: n
      integer(8)                :: kmax
      integer(8)                :: nx
      integer(8)                :: src
      real(8)                   :: c(2)

      integer(8)                :: s
      integer(8)                :: xn

      do sol=0,3
        do src=1,4
          do s=1,3
            if (s == 1) then
              c(1)=0.99d0
              c(2)=0.98d0
            elseif (s == 2) then 
              c(1)=0.95d0
              c(2)=0.90d0
            elseif (s == 3) then 
              c(1)=0.75d0
              c(2)=0.50d0
            endif
            n=6
            kmax=10000
            do xn=1,12
              nx=5*2**(xn-1)
              call solve_slab(sol,src,s,c,n,kmax,xn,nx) 
            enddo
          enddo
        enddo
      enddo

      contains

      subroutine solve_slab(sol,src,s,c,n,kmax,xn,nx) 

        implicit none

        integer(8), intent(in)    :: sol
        integer(8), intent(in)    :: src
        integer(8), intent(in)    :: s
        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: kmax
        integer(8), intent(in)    :: xn
        integer(8), intent(in)    :: nx
        real(8),    intent(in)    :: c(2)

        integer(8)                :: j
        integer(8)                :: jmax
        real(8)                   :: eps
        real(8)                   :: h
        real(8)                   :: x
        real(8)                   :: xm
        real(8)                   :: xp
        real(8)                   :: xj
        real(8)                   :: xmax
        real(8)                   :: xref
        real(8)                   :: pi
        real(8)                   :: mu(n/2)
        real(8)                   :: w (n/2)
        real(8), allocatable      :: phi (:)
        real(8), allocatable      :: q   (:)
        real(8), allocatable      :: sigt(:)
        real(8), allocatable      :: sigs(:)
        real(8), allocatable      :: xmsh(:)
        character(1)              :: prbopt
        character(1)              :: srcopt
        character(1)              :: solopt
        character(2)              :: xnopt
        character(13)             :: output

      ! mesh slab geometry

        xmax=40.0d0
        xref=10.0d0
        h   =xref/dble(nx)
        jmax=nint(xmax/h)

      ! dynamic allocation of arrays

        allocate(phi(jmax))
        allocate(q(jmax))
        allocate(sigt(jmax))
        allocate(sigs(jmax))
        allocate(xmsh(jmax))
        phi=0.0d0
        q=0.0d0
        sigt=0.0d0
        sigs=0.0d0
        xmsh=0.0d0

      ! build source based on options

        if (src == 1) then      ! flat
          q=1.0d0
        elseif (src == 2) then  ! hat
          do j=1,jmax
            xm=h*dble(j-1)
            xp=h*dble(j)
            xj=0.5d0*(xp+xm)
            q(j)=0.005d0*xj+0.95d0
            if (j > jmax/2) q(j)=-0.005d0*xj+1.15d0
          enddo
        elseif (src == 3) then  ! cosine
          pi=4.0d0*atan(1.0d0)
          do j=1,jmax
            xm=h*dble(j-1)
            xp=h*dble(j)
            xj=0.5d0*(xp+xm)
            q(j)=20.0d0/h*(cos(pi*xm/40.0d0)-cos(pi*xp/40.0d0))
          enddo
        elseif (src == 4) then  ! linear
          do j=1,jmax
            xm=h*dble(j-1)
            xp=h*dble(j)
            xj=0.5d0*(xp+xm)
            q(j)=xj/20.0d0
          enddo
        else
          write(0,'(a)') ' Incorrect source option selected.'
          stop
        endif

      ! cross sections from Kord's problem

        x=0.0d0
        do j=1,jmax
          xmsh(j)=x+h
          x=x+h
          sigt(j)=1.0d0
          sigs(j)=c(1)*sigt(j)
          if (j > nx .and. j < 3*nx+1) sigs(j)=c(2)*sigt(j)
        enddo

      ! set quadrature

        call quad(n,mu,w)

      ! solve fixed-source problem

        eps=5.0d-16
        if (sol == 0) then
          call solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)
        elseif (sol == 1) then
          call solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)
        elseif (sol == 2) then
          call solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)
        elseif (sol == 3) then
          call solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)
        else
          write(0,'(a)') ' Incorrect solution scheme selected.'
          stop
        endif

      ! write flux solution into file

        write(prbopt,'(i1)') s
        write(solopt,'(i1)') sol
        write(srcopt,'(i1)') src
        write(xnopt, '(i2)') xn
        output='p'//prbopt//'-'//srcopt//'-'//solopt//'-'//trim(adjustl(xnopt))//'.out'
        open(unit=1,file=output,action='write',status='unknown')
        write(1,*) jmax
        do j=1,jmax
          write(1,*) phi(j)
        enddo
        close(1)

      ! clean up arrays

        deallocate(phi)
        deallocate(q)
        deallocate(sigt)
        deallocate(sigs)
        deallocate(xmsh)

      end subroutine solve_slab

      subroutine solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(8)                :: j
        integer(8)                :: k
        integer(8)                :: m
        real(8)                   :: tau
        real(8)                   :: merr
        real(8)                   :: psi
        real(8)                   :: psi_in
        real(8)                   :: psi_bc(n/2)
        real(8),    allocatable   :: c1(:,:)
        real(8),    allocatable   :: c2(:,:)
        real(8),    allocatable   :: phio(:)
        real(8),    allocatable   :: s(:)

      ! pre-compute coeffs

        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        c1=0.0d0
        c2=0.0d0

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            c1(j,m)=0.5d0*tau
            c2(j,m)=0.5d0*h/mu(m)
          enddo
        enddo

      ! solve problem

        allocate(s(jmax))
        do j=1,jmax
          phi(j)=1.0d0
          s(j)  =0.5d0*(sigs(j)*phi(j)+q(j))
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=2.0d0*psi-psi_in
            enddo
            do j=jmax,1,-1
              psi =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=2.0d0*psi-psi_in
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/(phio(j)+1.0d-20)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s(j)=0.5d0*(sigs(j)*phi(j)+q(j))
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

        deallocate(c1)
        deallocate(c2)
        deallocate(s)
        deallocate(phio)

      end subroutine solve_dd

      subroutine solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(8)                :: j
        integer(8)                :: k
        integer(8)                :: m
        real(8)                   :: tau
        real(8)                   :: tau3
        real(8)                   :: tau5
        real(8)                   :: tau7
        real(8)                   :: merr
        real(8)                   :: psi
        real(8)                   :: psi_in
        real(8)                   :: psi_bc(n/2)
        real(8),    allocatable   :: alpha(:,:)
        real(8),    allocatable   :: c1(:,:)
        real(8),    allocatable   :: c2(:,:)
        real(8),    allocatable   :: phio(:)
        real(8),    allocatable   :: s(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0d0
        c1   =0.0d0
        c2   =0.0d0

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            if (tau < 0.01d0) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0d0-tau3/360.0d0+tau5/15120.0d0-tau7/604800.0d0
            else
              alpha(j,m)=1.0d0/tanh(tau/2.0d0)-2.0d0/tau
            endif
            c1(j,m)=0.5d0*tau*(1.0d0+alpha(j,m))
            c2(j,m)=0.5d0*h  *(1.0d0+alpha(j,m))/mu(m)
          enddo
        enddo

      ! solve problem

        allocate(s(jmax))
        do j=1,jmax
          phi(j)=1.0d0
          s(j)  =0.5d0*(sigs(j)*phi(j)+q(j))
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi=(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=(2.0d0*psi-(1.0d0-alpha(j,m))*psi_in)/(1.0d0+alpha(j,m))
            enddo
            do j=jmax,1,-1
              psi=(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=(2.0d0*psi-(1.0d0-alpha(j,m))*psi_in)/(1.0d0+alpha(j,m))
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/(phio(j)+1.0d-20)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s(j)=0.5d0*(sigs(j)*phi(j)+q(j))
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

        deallocate(alpha)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(s)

      end subroutine solve_sc

      subroutine solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(8)                :: j
        integer(8)                :: k
        integer(8)                :: m
        real(8)                   :: tau
        real(8)                   :: merr
        real(8)                   :: psi
        real(8)                   :: psil
        real(8)                   :: psi_in
        real(8)                   :: psi_out
        real(8)                   :: psi_bc(n/2)
        real(8),    allocatable   :: alpha(:,:)
        real(8),    allocatable   :: c1(:,:)
        real(8),    allocatable   :: c2(:,:)
        real(8),    allocatable   :: phio(:)
        real(8),    allocatable   :: phil(:)
        real(8),    allocatable   :: s(:)
        real(8),    allocatable   :: sl(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0d0
        c1   =0.0d0
        c2   =0.0d0

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            alpha(j,m)=1.0d0/(1.0d0+6.0d0/tau)
            c1(j,m)   =      (2.0d0/tau+alpha(j,m)-1.0d0)
            c2(j,m)   =1.0d0/(2.0d0/tau+alpha(j,m)+1.0d0)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        do j=1,jmax
          phi(j) =1.0d0
          phil(j)=0.0d0
          s (j)  =0.5d0*(sigs(j)*phi(j)+q(j))
          sl(j)  =0.0d0
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0d0
          phil=0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi_out=c2(j,m)*(2.0d0*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi=((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0-alpha(j,m)*sl(j)/sigt(j)
              psil=psi_out-psi
              psi_in=psi_out
              phi(j)=phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            do j=jmax,1,-1
              psi_out=c2(j,m)*(2.0d0*(s(j)-alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi=((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0+alpha(j,m)*sl(j)/sigt(j)
              psil=psi-psi_out
              psi_in=psi_out
              phi(j)=phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/(phio(j)+1.0d-20)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5d0*(sigs(j)*phi (j)+q(j))
            sl(j)=0.5d0*(sigs(j)*phil(j)  )
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

        deallocate(alpha)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(s)
        deallocate(sl)

      end subroutine solve_ld

      subroutine solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(8)                :: j
        integer(8)                :: k
        integer(8)                :: m
        real(8)                   :: tau
        real(8)                   :: tau3
        real(8)                   :: tau5
        real(8)                   :: tau7
        real(8)                   :: merr
        real(8)                   :: psi
        real(8)                   :: psil
        real(8)                   :: psi_in
        real(8)                   :: psi_out
        real(8)                   :: psi_bc(n/2)
        real(8),    allocatable   :: alpha(:,:)
        real(8),    allocatable   :: rbeta(:,:)
        real(8),    allocatable   :: c1(:,:)
        real(8),    allocatable   :: c2(:,:)
        real(8),    allocatable   :: phio(:)
        real(8),    allocatable   :: phil(:)
        real(8),    allocatable   :: s(:)
        real(8),    allocatable   :: sl(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(rbeta(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0d0
        rbeta=0.0d0
        c1   =0.0d0
        c2   =0.0d0

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            if (tau < 0.01d0) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0d0-tau3/360.0d0+tau5/15120.0d0-tau7/604800.0d0
            else
              alpha(j,m)=1.0d0/tanh(tau/2.0d0)-2.0d0/tau
            endif
            rbeta(j,m)=1.0d0/alpha(j,m)-6.0d0/tau
            c1(j,m)=      (2.0d0/tau+alpha(j,m)-1.0d0)
            c2(j,m)=1.0d0/(2.0d0/tau+alpha(j,m)+1.0d0)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        do j=1,jmax
          phi(j) =1.0d0
          phil(j)=0.0d0
          s (j)  =0.5d0*(sigs(j)*phi(j)+q(j))
          sl(j)  =0.0d0
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0d0
          phil=0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi_out=c2(j,m)*(2.0d0*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi=((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0-alpha(j,m)*sl(j)/sigt(j)
              psil=((rbeta(j,m)+1.0d0)*psi_out+(rbeta(j,m)-1.0d0)*psi_in)/2.0d0-rbeta(j,m)*psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            do j=jmax,1,-1
              psi_out=c2(j,m)*(2.0d0*(s(j)-alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi=((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0+alpha(j,m)*sl(j)/sigt(j)
              psil=-((rbeta(j,m)+1.0d0)*psi_out+(rbeta(j,m)-1.0d0)*psi_in)/2.0d0+rbeta(j,m)*psi
              psi_in=psi_out
              phi(j)=phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5d0*(sigs(j)*phi (j)+q(j))
            sl(j)=0.5d0*(sigs(j)*phil(j)  )
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

        deallocate(alpha)
        deallocate(rbeta)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(s)
        deallocate(sl)

      end subroutine solve_lc

      subroutine quad(n,mu,w)

        implicit none

        integer(8), intent(in)    :: n
        real(8),    intent(out)   :: mu(n/2)
        real(8),    intent(out)   :: w (n/2)

        if (n == 2) then
          mu(1)=0.5773502691d0
          w (1)=1.0000000000d0
        elseif (n == 4) then
          mu(1)=0.3399810435d0
          w (1)=0.6521451549d0
          mu(2)=0.8611363115d0
          w (2)=0.3478548451d0
        elseif (n == 6) then
          mu(1)=0.2386191860d0
          w (1)=0.4679139346d0
          mu(2)=0.6612093864d0
          w (2)=0.3607615730d0
          mu(3)=0.9324695142d0
          w (3)=0.1713244924d0
        elseif (n == 8) then
          mu(1)=0.1834346424d0
          w (1)=0.3626837834d0
          mu(2)=0.5255324099d0
          w (2)=0.3137066459d0
          mu(3)=0.7966664774d0
          w (3)=0.2223810344d0
          mu(4)=0.9602898564d0
          w (4)=0.1012285363d0
        elseif (n == 16) then
          mu(1)=0.0950125098d0
          w (1)=0.1894506105d0
          mu(2)=0.2816035508d0
          w (2)=0.1826034150d0
          mu(3)=0.4580167777d0
          w (3)=0.1691565194d0
          mu(4)=0.6178762444d0
          w (4)=0.1495959888d0
          mu(5)=0.7554044084d0
          w (5)=0.1246289713d0
          mu(6)=0.8656312024d0
          w (6)=0.0951585117d0
          mu(7)=0.9445750231d0
          w (7)=0.0622535239d0
          mu(8)=0.9894009350d0
          w (8)=0.0271524594d0
        endif

      end subroutine quad

    end program main
