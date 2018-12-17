    program main

      implicit none

      call drive

      contains

      subroutine drive

        implicit none

        integer(4)                :: sol
        integer(4)                :: n
        integer(4)                :: kmax
        integer(4)                :: nx
        integer(4)                :: src
        real(8)                   :: c(2)

        integer(4)                :: prb
        integer(4)                :: xn
        integer(4)                :: xnmax
        integer(4), parameter     :: xnr=18

        do sol=0,3
          xnmax=12
          if (sol == 2) then
            xnmax=9
          elseif (sol == 3) then
            xnmax=7
          endif
          do src=1,2
            do prb=1,4
              if (prb == 4 .and. src == 2) cycle
              if (prb == 1) then
                c(1)=0.99d0
                c(2)=0.98d0
              elseif (prb == 2) then 
                c(1)=0.95d0
                c(2)=0.90d0
              elseif (prb == 3) then 
                c(1)=0.75d0
                c(2)=0.50d0
              elseif (prb == 4) then
                c(1)=0.90d0
                c(2)=0.60d0
              endif
              n=6
              if (prb == 4) n=10
              kmax=10000
              do xn=1,xnmax
                nx=5*2**(xn-1)
                if (prb == 4) nx=100*2**(xn-1)
                call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
              enddo
            enddo
          enddo
          do src=3,4
            prb=1
            c=0.98d0
            n=6
            kmax=10000
            xnmax=12
            if (sol == 2) then
              xnmax=9
            elseif (sol == 3) then
              xnmax=7
            endif
            do xn=1,xnmax
              nx=5*2**(xn-1)
              call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
            enddo
          enddo
        enddo

      ! reference case (nxr) using dd (sol=0)

        sol=0
        do src=1,4
          if (src == 1 .or. src == 2) then
            do prb=1,4
              if (prb == 4 .and. src == 2) cycle
              if (prb == 1) then
                c(1)=0.99d0
                c(2)=0.98d0
              elseif (prb == 2) then 
                c(1)=0.95d0
                c(2)=0.90d0
              elseif (prb == 3) then 
                c(1)=0.75d0
                c(2)=0.50d0
              elseif (prb == 4) then
                c(1)=0.90d0
                c(2)=0.60d0
              endif
              n=6
              if (prb == 4) n=10
              kmax=10000
              xn=xnr
              nx=5*2**(xn-1)
              if (prb == 4) nx=100*2**(xn-1)
              call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
            enddo
          elseif (src == 3 .or. src == 4) then
            prb=1
            c=0.98d0
            n=6
            kmax=10000
            xn=xnr
            nx=5*2**(xn-1)
            call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
          endif
        enddo

      end subroutine drive

      subroutine solve_slab(sol,src,prb,c,n,kmax,xn,nx) 

        implicit none

        integer(4), intent(in)    :: sol
        integer(4), intent(in)    :: src
        integer(4), intent(in)    :: prb
        integer(4), intent(in)    :: n
        integer(4), intent(in)    :: kmax
        integer(4), intent(in)    :: xn
        integer(4), intent(in)    :: nx
        real(8),    intent(in)    :: c(2)

        integer(4)                :: j
        integer(4)                :: jmax
        integer(4)                :: bc(2)
        real(8)                   :: eps
        real(8)                   :: h
        real(8)                   :: x
        real(8)                   :: xm
        real(8)                   :: xp
        real(8)                   :: xj
        real(8)                   :: xmax
        real(8)                   :: xref
        real(8)                   :: pi
        real(8)                   :: xnorm
        real(8)                   :: mu(n/2)
        real(8)                   :: w (n/2)
        real(8), allocatable      :: phi (:)
        real(8), allocatable      :: q   (:)
        real(8), allocatable      :: ql  (:)
        real(8), allocatable      :: sigt(:)
        real(8), allocatable      :: sigs(:)
        real(8), allocatable      :: xmsh(:)
        character(1)              :: prbopt
        character(1)              :: srcopt
        character(1)              :: solopt
        character(2)              :: xnopt
        character(13)             :: datafile

      ! mesh slab geometry

        xmax=40.0d0
        xref=10.0d0
        if (prb == 4) then
          xmax=8.0d0
          xref=2.0d0
        endif
        h   =xref/dble(nx)
        jmax=nint(xmax/h)

      ! dynamic allocation of arrays

        allocate(phi(jmax))
        allocate(q(jmax))
        allocate(ql(jmax))
        allocate(sigt(jmax))
        allocate(sigs(jmax))
        allocate(xmsh(jmax))
        phi=0.0d0
        q=0.0d0
        ql=0.0d0
        sigt=0.0d0
        sigs=0.0d0
        xmsh=0.0d0

      ! build source based on options

        if (src == 1) then      ! flat
          q=1.0d0
          bc=1
          if (prb == 4) bc=0
        elseif (src == 2) then  ! hat
          do j=1,jmax
            xm=h*dble(j-1)
            xp=h*dble(j)
            xj=0.5d0*(xp+xm)
            q (j)=0.005d0*xj+0.95d0
            ql(j)=0.0025d0*h
            if (j > jmax/2) then
              q (j)=-0.005d0*xj+1.15d0
              ql(j)=-0.0025d0*h
            endif
          enddo
          bc=1
        elseif (src == 3) then  ! cosine
          pi=4.0d0*atan(1.0d0)
          do j=1,jmax
            xm=h*dble(j-1)
            xp=h*dble(j)
            xj=0.5d0*(xp+xm)
            q (j)=40.0d0/h*sin(pi*h/80.0d0)*sin(pi*xj/40.0d0)
            ql(j)=-120.0d0/(pi*h*h)*cos(pi*xj/40.0d0)*&
                  (pi*h*cos(pi*h/80.0d0)-80.0d0*sin(pi*h/80.0d0))
          enddo
          bc=1
        elseif (src == 4) then  ! linear
          do j=1,jmax
            xm=h*dble(j-1)
            xp=h*dble(j)
            xj=0.5d0*(xp+xm)
            q (j)=xj/20.0d0
            ql(j)=h/40.0d0
          enddo
          bc=0
        else
          write(0,'(a)') ' Incorrect source option selected.'
          stop
        endif

      ! set cross sections from Kord's problem

        x=0.0d0
        do j=1,jmax
          xmsh(j)=x+0.5d0*h
          x=x+h
          sigt(j)=1.0d0
          sigs(j)=c(1)*sigt(j)
          if (prb == 4) then
            if (xmsh(j) > xref .and. xmsh(j) < xmax-xref) then
              sigt(j)=100.0d0
              sigs(j)=c(2)*sigt(j)
            endif
          else
            if (j > nx .and. j < 3*nx+1) sigs(j)=c(2)*sigt(j)
          endif
        enddo

      ! set quadrature

        call quad(n,mu,w)

      ! solve fixed-source problem

        eps=1.0d-15
        if (sol == 0) then
          call solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi)
        elseif (sol == 1) then
          call solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi)
        elseif (sol == 2) then
          call solve_ld(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi)
        elseif (sol == 3) then
          call solve_lc(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi)
        else
          write(0,'(a)') ' Incorrect solution scheme selected.'
          stop
        endif

      ! Renormalize results

        xnorm=0.01d0
        if (prb == 2) then
          xnorm=0.05d0
        elseif (prb == 3) then
          xnorm=0.25d0
        elseif (prb == 4) then
          xnorm=1.0d0
        endif

      ! write flux solution into file

        write(prbopt,'(i1)') prb
        write(solopt,'(i1)') sol
        write(srcopt,'(i1)') src
        write(xnopt, '(i2)') xn
        datafile='p'//prbopt//'-'//srcopt//'-'//solopt//'-'//trim(adjustl(xnopt))//'.dat'
        open(unit=1,file=datafile,action='write',status='unknown')
        do j=1,jmax
          write(1,'(3(es25.16))') xmsh(j),phi(j)*xnorm,q(j)
        enddo
        close(1)

      ! clean up arrays

        deallocate(phi)
        deallocate(q)
        deallocate(ql)
        deallocate(sigt)
        deallocate(sigs)
        deallocate(xmsh)

      end subroutine solve_slab

      subroutine solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi)

        implicit none

        integer(4), intent(in)    :: n
        integer(4), intent(in)    :: jmax
        integer(4), intent(in)    :: kmax
        integer(4), intent(in)    :: bc(2)
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(4)                :: j
        integer(4)                :: k
        integer(4)                :: m
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
            if (bc(1) == 0) psi_in=0.0d0
            do j=1,jmax
              psi =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=2.0d0*psi-psi_in
            enddo
            if (bc(2) == 0) psi_in=0.0d0
            do j=jmax,1,-1
              psi =(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=2.0d0*psi-psi_in
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/(phi(j)+1.0d-20)))
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

      subroutine solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi)

        implicit none

        integer(4), intent(in)    :: n
        integer(4), intent(in)    :: jmax
        integer(4), intent(in)    :: kmax
        integer(4), intent(in)    :: bc(2)
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(4)                :: j
        integer(4)                :: k
        integer(4)                :: m
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
              alpha(j,m)=1.0d0/dtanh(tau/2.0d0)-2.0d0/tau
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
            if (bc(1) == 0) psi_in=0.0d0
            do j=1,jmax
              psi=(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=(2.0d0*psi-(1.0d0-alpha(j,m))*psi_in)/(1.0d0+alpha(j,m))
            enddo
            if (bc(2) == 0) psi_in=0.0d0
            do j=jmax,1,-1
              psi=(s(j)*c2(j,m)+psi_in)/(1.0d0+c1(j,m))
              phi(j)=phi(j)+psi*w(m)
              psi_in=(2.0d0*psi-(1.0d0-alpha(j,m))*psi_in)/(1.0d0+alpha(j,m))
            enddo
            psi_bc(m)=psi_in
          enddo

          merr=0.0d0
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/(phi(j)+1.0d-20)))
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

      subroutine solve_ld(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi)

        implicit none

        integer(4), intent(in)    :: n
        integer(4), intent(in)    :: jmax
        integer(4), intent(in)    :: kmax
        integer(4), intent(in)    :: bc(2)
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: ql(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(4)                :: j
        integer(4)                :: k
        integer(4)                :: m
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
            if (bc(1) == 0) psi_in=0.0d0
            do j=1,jmax
              psi_out=c2(j,m)*(2.0d0*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi=((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0-alpha(j,m)*sl(j)/sigt(j)
              psil=psi_out-psi
              psi_in=psi_out
              phi(j)=phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            if (bc(2) == 0) psi_in=0.0d0
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
            merr=max(merr,abs((phi(j)-phio(j))/(phi(j)+1.0d-20)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5d0*(sigs(j)*phi (j)+q (j))
            sl(j)=0.5d0*(sigs(j)*phil(j)+ql(j))
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

      subroutine solve_lc(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi)

        implicit none

        integer(4), intent(in)    :: n
        integer(4), intent(in)    :: jmax
        integer(4), intent(in)    :: kmax
        integer(4), intent(in)    :: bc(2)
        real(8),    intent(in)    :: h
        real(8),    intent(in)    :: q(jmax)
        real(8),    intent(in)    :: ql(jmax)
        real(8),    intent(in)    :: eps
        real(8),    intent(in)    :: sigt(jmax)
        real(8),    intent(in)    :: sigs(jmax)
        real(8),    intent(in)    :: mu(n/2)
        real(8),    intent(in)    :: w (n/2)
        real(8),    intent(out)   :: phi (jmax)

        integer(4)                :: j
        integer(4)                :: k
        integer(4)                :: m
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
              alpha(j,m)=1.0d0/dtanh(tau/2.0d0)-2.0d0/tau
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
            if (bc(1) == 0) psi_in=0.0d0
            do j=1,jmax
              psi_out=c2(j,m)*(2.0d0*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi=((1.0d0+alpha(j,m))*psi_out+(1.0d0-alpha(j,m))*psi_in)/2.0d0-alpha(j,m)*sl(j)/sigt(j)
              psil=((rbeta(j,m)+1.0d0)*psi_out+(rbeta(j,m)-1.0d0)*psi_in)/2.0d0-rbeta(j,m)*psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            if (bc(2) == 0) psi_in=0.0d0
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
            merr=max(merr,abs((phi(j)-phio(j))/(phi(j)+1.0d-20)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5d0*(sigs(j)*phi (j)+q (j))
            sl(j)=0.5d0*(sigs(j)*phil(j)+ql(j))
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

        integer(4), intent(in)    :: n
        real(8),    intent(out)   :: mu(n/2)
        real(8),    intent(out)   :: w (n/2)

        integer(4)                :: j
        integer(4), parameter     :: nmaxp=300
        real(8)                   :: xnew(nmaxp)
        real(8)                   :: wnew(nmaxp)

        xnew=0.0d0
        wnew=0.0d0
        call gauleg(-1.0d0,1.0d0,xnew,wnew,n)

        do j=1,n/2
          mu(j)=xnew(j)
          w(j) =wnew(j)
        enddo

      end subroutine quad

      subroutine gauleg(x1,x2,x,w,n)
 
        implicit none

        integer(4), intent(in)    :: n
        real(8),    intent(in)    :: x1
        real(8),    intent(in)    :: x2
        real(8),    intent(inout) :: x(n)
        real(8),    intent(inout) :: w(n)

        integer(4)                :: i
        integer(4)                :: j
        integer(4)                :: m
        integer(4)                :: kount
        integer(4), parameter     :: nmax=300
        real(8)                   :: xm
        real(8)                   :: xl
        real(8)                   :: p1
        real(8)                   :: p2
        real(8)                   :: p3
        real(8)                   :: pi
        real(8)                   :: pp
        real(8)                   :: z
        real(8)                   :: z1
        real(8)                   :: xtmp(nmax)  ! full set of abscissas
        real(8)                   :: wtmp(nmax)  ! full set of weights
        real(8),    parameter     :: eps=1.0d-15

        pi=4.0d0*atan(1.0d0)

        if (n > nmax) then
          write(0,'(a,1i6)') 'Gauss-Leg. integration problem --Increase PARAMETER: NMAX to at least:',n
          stop
        endif

        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)

        do i=1,m
          z=cos(pi*(i-0.25d0)/(n+0.5d0))
      1   continue
          p1=1.d0
          p2=0.d0
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
          enddo
      !   p1 is now the desired Legendre polynomial. we next compute pp, its derivative,
      !   by a standard relation involving also p2, the polynomial of one lower order.
          pp=n*(z*p1-p2)/(z*z-1.0d0)
          z1=z
          z=z1-p1/pp
          if (abs(z-z1) > eps) go to 1
          xtmp(i)=    xm-xl*z
          xtmp(n+1-i)=xm+xl*z
      !   the (n+1-i) terms are the symmetric counterparts
          wtmp(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
          wtmp(n+1-i)=wtmp(i)
        enddo

      ! (half set and assumed symmetric)
        kount=0
        do i=1,n
          if (xtmp(i) >= 0.0) then
            kount=kount+1
            x(kount)=xtmp(i)   ! abscissas
            w(kount)=wtmp(i)   ! weights
          endif
        enddo

      end subroutine gauleg

    end program main
