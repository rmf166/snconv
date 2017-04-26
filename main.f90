    program main

      implicit none

      integer(8)                  :: sol
      integer(8)                  :: n
      integer(8)                  :: kmax
      integer(8)                  :: nx
      real(8)                     :: q
      real(8)                     :: c(2)

      call read_input(sol,q,c,n,kmax,nx) 

      call solve_slab(sol,q,c,n,kmax,nx) 

      contains

      subroutine read_input(sol,q,c,n,kmax,nx) 

        implicit none

        integer(8), intent(out)   :: sol
        integer(8), intent(out)   :: n
        integer(8), intent(out)   :: nx
        integer(8), intent(out)   :: kmax
        real(8),    intent(out)   :: q
        real(8),    intent(out)   :: c(2)

        open(unit=1,file='input',action='read',status='unknown')
        read(1,*) sol
        read(1,*) q
        read(1,*) c(1)
        read(1,*) c(2)
        read(1,*) n
        read(1,*) kmax
        read(1,*) nx
        close(1)
        nx=5**nx

      end subroutine read_input

      subroutine solve_slab(sol,q,c,n,kmax,nx) 

        implicit none

        integer(8), intent(in)    :: kmax
        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: nx
        integer(8), intent(in)    :: sol
        real(8),    intent(in)    :: q
        real(8),    intent(in)    :: c(2)

        integer(8)                :: j
        integer(8)                :: jmax
        real(8)                   :: eps
        real(8)                   :: h
        real(8)                   :: x
        real(8)                   :: xmax
        real(8)                   :: xref
        real(8)                   :: mu(n/2)
        real(8)                   :: w (n/2)
        real(8), allocatable      :: phi (:)
        real(8), allocatable      :: sigt(:)
        real(8), allocatable      :: sigs(:)
        real(8), allocatable      :: xmsh(:)

      ! mesh slab geometry

        xmax=40.0d0
        xref=10.0d0
        h   =xref/dble(nx)
        jmax=nint(xmax/h)

      ! cross sections from Kord's problem

        allocate(phi(jmax))
        allocate(sigt(jmax))
        allocate(sigs(jmax))
        allocate(xmsh(jmax))
        phi=0.0d0
        sigt=0.0d0
        sigs=0.0d0
        xmsh=0.0d0

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

      ! write solution into file

        open(unit=1,file='output',action='write',status='unknown')
        write(1,*) jmax
        do j=1,jmax
          write(1,*) phi(j)
        enddo
        close(1)

      end subroutine solve_slab

      subroutine solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q
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
        real(8)                   :: c1   (jmax,n/2)
        real(8)                   :: c2   (jmax,n/2)
        real(8)                   :: psi_bc(n/2)
        real(8)                   :: phio(jmax)
        real(8)                   :: s(jmax)

      ! pre-compute coeffs

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

        do j=1,jmax
          phi(j)=1.0d0
          s(j)  =0.5d0*(sigs(j)*phi(j)+q)
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        do k=1,kmax
          phio=phi
          phi =0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi   =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=2.0d0*psi-psi_in
            enddo
            do j=jmax,1,-1
              psi   =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
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
            s(j)=0.5d0*(sigs(j)*phi(j)+q)
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

      end subroutine solve_dd

      subroutine solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q
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
        real(8)                   :: alpha(jmax,n/2)
        real(8)                   :: c1   (jmax,n/2)
        real(8)                   :: c2   (jmax,n/2)
        real(8)                   :: psi_bc(n/2)
        real(8)                   :: phio(jmax)
        real(8)                   :: s(jmax)

      ! pre-compute coeffs

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

        do j=1,jmax
          phi(j)=1.0d0
          s(j)  =0.5d0*(sigs(j)*phi(j)+q)
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        do k=1,kmax
          phio=phi
          phi =0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi   =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=(2.0d0*psi-(1.0d0-alpha(j,m))*psi_in)/(1.0d0+alpha(j,m))
            enddo
            do j=jmax,1,-1
              psi   =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
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
            s(j)=0.5d0*(sigs(j)*phi(j)+q)
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

      end subroutine solve_sc

      subroutine solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q
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
        real(8)                   :: alpha(jmax,n/2)
        real(8)                   :: c1   (jmax,n/2)
        real(8)                   :: c2   (jmax,n/2)
        real(8)                   :: psi_bc(n/2)
        real(8)                   :: phil(jmax)
        real(8)                   :: phio(jmax)
        real(8)                   :: s (jmax)
        real(8)                   :: sl(jmax)

      ! pre-compute coeffs

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

        do j=1,jmax
          phi(j) =1.0d0
          phil(j)=0.0d0
          s (j)  =0.5d0*(sigs(j)*phi(j)+q)
          sl(j)  =0.0d0
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        do k=1,kmax
          phio=phi
          phi =0.0d0
          phil=0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi_out=c2(j,m)*(2.0d0*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi    =((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0-alpha(j,m)*sl(j)/sigt(j)
              psil   =psi_out-psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            do j=jmax,1,-1
              psi_out   =c2(j,m)*(2.0d0*(s(j)-alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi       =((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0+alpha(j,m)*sl(j)/sigt(j)
              psil      =psi-psi_out
              psi_in    =psi_out
              phi(j)    =phi(j) +psi*w(m)
              phil(j)   =phil(j)+psil*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/(phio(j)+1.0d-20)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5d0*(sigs(j)*phi (j)+q)
            sl(j)=0.5d0*(sigs(j)*phil(j)  )
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

      end subroutine solve_ld

      subroutine solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,phi)

        implicit none

        integer(8), intent(in)    :: n
        integer(8), intent(in)    :: jmax
        integer(8), intent(in)    :: kmax
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q
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
        real(8)                   :: alpha(jmax,n/2)
        real(8)                   :: rbeta(jmax,n/2)
        real(8)                   :: c1   (jmax,n/2)
        real(8)                   :: c2   (jmax,n/2)
        real(8)                   :: psi_bc(n/2)
        real(8)                   :: phil(jmax)
        real(8)                   :: phio(jmax)
        real(8)                   :: s (jmax)
        real(8)                   :: sl(jmax)

      ! pre-compute coeffs

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

        do j=1,jmax
          phi(j) =1.0d0
          phil(j)=0.0d0
          s (j)  =0.5d0*(sigs(j)*phi(j)+q)
          sl(j)  =0.0d0
        enddo

        psi_in=0.0d0
        psi_bc=0.0d0

        do k=1,kmax
          phio=phi
          phi =0.0d0
          phil=0.0d0
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            do j=1,jmax
              psi_out=c2(j,m)*(2.0d0*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi    =((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0-alpha(j,m)*sl(j)/sigt(j)
              psil   =((rbeta(j,m)+1.0d0)*psi_out+(rbeta(j,m)-1.0d0)*psi_in)/2.0d0-rbeta(j,m)*psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            do j=jmax,1,-1
              psi_out   =c2(j,m)*(2.0d0*(s(j)-alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi       =((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0+alpha(j,m)*sl(j)/sigt(j)
              psil      =-((rbeta(j,m)+1.0d0)*psi_out+(rbeta(j,m)-1.0d0)*psi_in)/2.0d0+rbeta(j,m)*psi
              psi_in    =psi_out
              phi(j)    =phi(j) +psi*w(m)
              phil(j)   =phil(j)+psil*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5d0*(sigs(j)*phi (j)+q)
            sl(j)=0.5d0*(sigs(j)*phil(j)  )
          enddo
        enddo

        if (k >= kmax) then
          write(0,'(a,i8,a,es12.5)') ' Source iteration did not converge, k = ',k,' err = ',merr
          stop
        endif

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
