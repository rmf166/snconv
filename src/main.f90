    module global

      implicit none

      integer(4),    parameter     ::                    &
          sp = kind(1.0),                                &
          dp = selected_real_kind(2*precision(1.0_sp)),  &
          qp = selected_real_kind(2*precision(1.0_dp))
      integer(4),    parameter     :: kr=qp

    end module global

    program main

      implicit none

      !$ call omp_set_num_threads(4)

      if (.true.) then
        call drive_kss
      else
        call drive_dean
      endif

      contains

      subroutine drive_kss

        use global
        !$ use omp_lib, only: omp_set_num_threads

        implicit none

        integer(4)                   :: sol
        integer(4),    parameter     :: n=6
        integer(4),    parameter     :: kmax=10000
        integer(4)                   :: nx
        integer(4)                   :: src
        real(kind=kr)                :: c(2)

        integer(4)                   :: prb
        integer(4)                   :: xn
        integer(4)                   :: xnmax
        integer(4),    parameter     :: xnr=18

        !$ call omp_set_num_threads(n/2)

        do sol=0,3
          xnmax=12
          if (sol == 2) then
            xnmax=9
          elseif (sol == 3) then
            xnmax=7
          endif
          do src=1,2
            do prb=1,3
              if (prb == 1) then
                c(1)=0.99_kr
                c(2)=0.98_kr
              elseif (prb == 2) then
                c(1)=0.95_kr
                c(2)=0.90_kr
              elseif (prb == 3) then
                c(1)=0.75_kr
                c(2)=0.50_kr
              endif
              do xn=1,xnmax
                nx=5*2**(xn-1)
                call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
              enddo
            enddo
          enddo
          do src=3,4
            prb=1
            c=0.98_kr
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
            do prb=1,3
              if (prb == 1) then
                c(1)=0.99_kr
                c(2)=0.98_kr
              elseif (prb == 2) then 
                c(1)=0.95_kr
                c(2)=0.90_kr
              elseif (prb == 3) then 
                c(1)=0.75_kr
                c(2)=0.50_kr
              endif
              xn=xnr
              nx=5*2**(xn-1)
              call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
            enddo
          elseif (src == 3 .or. src == 4) then
            prb=1
            c=0.98_kr
            xn=xnr
            nx=5*2**(xn-1)
            call solve_slab(sol,src,prb,c,n,kmax,xn,nx)
          endif
        enddo

      end subroutine drive_kss

      subroutine drive_dean

        use global
        !$ use omp_lib, only: omp_set_num_threads

        implicit none

        integer(4)                   :: sol
        integer(4),    parameter     :: n=12
        integer(4),    parameter     :: kmax=10000
        integer(4)                   :: nx
        real(kind=kr)                :: c(2)

        integer(4)                   :: xn
        integer(4)                   :: xnmax
        integer(4),    parameter     :: xnr=18

        ! some parameters

        c(1)=0.50_kr
        c(2)=0.00_kr

        ! numerical solutions

        !$ call omp_set_num_threads(n/2)

        do sol=0,3
          xnmax=14
          do xn=1,xnmax
            nx=5*2**(xn-1)
            call solve_slab(sol,1,4,c,n,kmax,xn,nx)
          enddo
        enddo

      ! reference case (nxr) using lc (sol=3)

        sol=3
        xn=xnr
        nx=5*2**(xn-1)
        call solve_slab(sol,1,4,c,n,kmax,xn,nx)

      end subroutine drive_dean

      subroutine solve_slab(sol,src,prb,c,n,kmax,xn,nx)

        use global

        implicit none

        integer(4),    intent(in)    :: sol
        integer(4),    intent(in)    :: src
        integer(4),    intent(in)    :: prb
        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: xn
        integer(4),    intent(in)    :: nx
        real(kind=kr), intent(in)    :: c(2)

        integer(4)                   :: j
        integer(4)                   :: jmax
        integer(4)                   :: bc(2)
        real(kind=kr)                :: eps
        real(kind=kr)                :: h
        real(kind=kr)                :: x
        real(kind=kr)                :: xm
        real(kind=kr)                :: xp
        real(kind=kr)                :: xj
        real(kind=kr)                :: xmax
        real(kind=kr)                :: xref
        real(kind=kr)                :: pi
        real(kind=kr)                :: xnorm
        real(kind=kr)                :: mu(n/2)
        real(kind=kr)                :: w (n/2)
        real(kind=kr), allocatable   :: phi (:)
        real(kind=kr), allocatable   :: jnet(:)
        real(kind=kr), allocatable   :: q   (:)
        real(kind=kr), allocatable   :: ql  (:)
        real(kind=kr), allocatable   :: sigt(:)
        real(kind=kr), allocatable   :: sigs(:)
        real(kind=kr), allocatable   :: xmsh(:)
        character(1)                 :: prbopt
        character(1)                 :: srcopt
        character(1)                 :: solopt
        character(2)                 :: xnopt
        character(13)                :: datafile

      ! mesh slab geometry

        xmax=40.0_kr
        xref=10.0_kr
        if (prb == 4) then
          xmax=8.0_kr
          xref=2.0_kr
        endif
        h   =xref/nx
        jmax=nint(xmax/h)

      ! dynamic allocation of arrays

        allocate(phi(jmax))
        allocate(jnet(jmax+1))
        allocate(q(jmax))
        allocate(ql(jmax))
        allocate(sigt(jmax))
        allocate(sigs(jmax))
        allocate(xmsh(jmax))
        phi=0.0_kr
        jnet=0.0_kr
        q=0.0_kr
        ql=0.0_kr
        sigt=0.0_kr
        sigs=0.0_kr
        xmsh=0.0_kr

      ! build source based on options

        if (src == 1) then      ! flat
          q=1.0_kr
          bc=1
          if (prb == 4) bc=0
        elseif (src == 2) then  ! hat
          do j=1,jmax
            xm=h*(j-1)
            xp=h*(j)
            xj=0.5_kr*(xp+xm)
            q (j)=0.005_kr*xj+0.95_kr
            ql(j)=0.0025_kr*h
            if (j > jmax/2) then
              q (j)=-0.005_kr*xj+1.15_kr
              ql(j)=-0.0025_kr*h
            endif
          enddo
          bc=1
        elseif (src == 3) then  ! cosine
          pi=4.0_kr*atan(1.0_kr)
          do j=1,jmax
            xm=h*(j-1)
            xp=h*(j)
            xj=0.5_kr*(xp+xm)
            q (j)=40.0_kr/h*sin(pi*h/80.0_kr)*sin(pi*xj/40.0_kr)
            ql(j)=-120.0_kr/(pi*h*h)*cos(pi*xj/40.0_kr)*&
                  (pi*h*cos(pi*h/80.0_kr)-80.0_kr*sin(pi*h/80.0_kr))
          enddo
          bc=1
        elseif (src == 4) then  ! linear
          do j=1,jmax
            xm=h*(j-1)
            xp=h*(j)
            xj=0.5_kr*(xp+xm)
            q (j)=xj/20.0_kr
            ql(j)=h/40.0_kr
          enddo
          bc=0
        else
          write(0,'(a)') ' Incorrect source option selected.'
          stop
        endif

      ! set cross sections from Kord's problem

        x=0.0_kr
        do j=1,jmax
          xmsh(j)=x+0.5_kr*h
          x=x+h
          sigt(j)=1.0_kr
          sigs(j)=c(1)*sigt(j)
          if (prb == 4) then
            if (xmsh(j) > xref .and. xmsh(j) < xmax-xref) then
              q(j)   =0.0_kr
              sigt(j)=100.0_kr
              sigs(j)=c(2)*sigt(j)
            endif
          else
            if (j > nx .and. j < 3*nx+1) sigs(j)=c(2)*sigt(j)
          endif
        enddo

      ! set quadrature

        call quad(n,mu,w)

      ! solve fixed-source problem

        eps=1.0e-30
        if (sol == 0) then
          call solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet)
        elseif (sol == 1) then
          call solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet)
        elseif (sol == 2) then
          call solve_ld(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi,jnet)
        elseif (sol == 3) then
          call solve_lc(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi,jnet)
        else
          write(0,'(a)') ' Incorrect solution scheme selected.'
          stop
        endif

      ! Renormalize results

        xnorm=0.01_kr
        if (prb == 2) then
          xnorm=0.05_kr
        elseif (prb == 3) then
          xnorm=0.25_kr
        elseif (prb == 4) then
          xnorm=1.0_kr
        endif

      ! write flux solution into file

        write(prbopt,'(i1)') prb
        write(solopt,'(i1)') sol
        write(srcopt,'(i1)') src
        write(xnopt, '(i2)') xn
        datafile='p'//prbopt//'-'//srcopt//'-'//solopt//'-'//trim(adjustl(xnopt))//'.dat'
        open(unit=1,file=datafile,action='write',status='unknown')
        do j=1,jmax
          write(1,*) xmsh(j),phi(j)*xnorm,q(j),jnet(j)
        enddo
        close(1)

      ! clean up arrays

        deallocate(phi)
        deallocate(jnet)
        deallocate(q)
        deallocate(ql)
        deallocate(sigt)
        deallocate(sigs)
        deallocate(xmsh)

      end subroutine solve_slab

      subroutine quad(n,mu,w)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(out)   :: mu(n/2)
        real(kind=kr), intent(out)   :: w (n/2)

        integer(4)                   :: j
        integer(4),    parameter     :: nmaxp=300
        real(kind=kr)                :: xnew(nmaxp)
        real(kind=kr)                :: wnew(nmaxp)

        xnew=0.0_kr
        wnew=0.0_kr
        call gauleg(-1.0_kr,1.0_kr,xnew,wnew,n)

        do j=1,n/2
          mu(j)=xnew(j)
          w(j) =wnew(j)
        enddo

      end subroutine quad

      subroutine gauleg(x1,x2,x,w,n)
 
        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(in)    :: x1
        real(kind=kr), intent(in)    :: x2
        real(kind=kr), intent(inout) :: x(n)
        real(kind=kr), intent(inout) :: w(n)

        integer(4)                   :: i
        integer(4)                   :: j
        integer(4)                   :: m
        integer(4)                   :: kount
        integer(4),    parameter     :: nmax=300
        real(kind=kr)                :: xm
        real(kind=kr)                :: xl
        real(kind=kr)                :: p1
        real(kind=kr)                :: p2
        real(kind=kr)                :: p3
        real(kind=kr)                :: pi
        real(kind=kr)                :: pp
        real(kind=kr)                :: z
        real(kind=kr)                :: z1
        real(kind=kr)                :: xtmp(nmax)  ! full set of abscissas
        real(kind=kr)                :: wtmp(nmax)  ! full set of weights
        real(kind=kr), parameter     :: eps=1.0e-30

        pi=4.0_kr*atan(1.0_kr)
        if (n > nmax) then
          write(0,'(a,1i6)') 'Gauss-Leg. integration problem --Increase PARAMETER: NMAX to at least:',n
          stop
        endif

        m=(n+1)/2
        xm=0.5_kr*(x2+x1)
        xl=0.5_kr*(x2-x1)
        do i=1,m
          z=cos(pi*(i-0.25_kr)/(n+0.5_kr))
      1   continue
          p1=1.0_kr
          p2=0.0_kr
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0_kr*j-1.0_kr)*z*p2-(j-1.0_kr)*p3)/j
          enddo
      !   p1 is now the desired Legendre polynomial. we next compute pp, its derivative,
      !   by a standard relation involving also p2, the polynomial of one lower order.
          pp=n*(z*p1-p2)/(z*z-1.0_kr)
          z1=z
          z=z1-p1/pp
          if (abs(z-z1) > eps) go to 1
          xtmp(i)=    xm-xl*z
          xtmp(n+1-i)=xm+xl*z
      !   the (n+1-i) terms are the symmetric counterparts
          wtmp(i)=2.0_kr*xl/((1.0_kr-z*z)*pp*pp)
          wtmp(n+1-i)=wtmp(i)
        enddo

      ! (half set and assumed symmetric)
        kount=0
        do i=1,n
          if (xtmp(i) >= 0.0_kr) then
            kount=kount+1
            x(kount)=xtmp(i)   ! abscissas
            w(kount)=wtmp(i)   ! weights
          endif
        enddo

      end subroutine gauleg

      subroutine solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q(jmax)
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt(jmax)
        real(kind=kr), intent(in)    :: sigs(jmax)
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(out)   :: phi (jmax)
        real(kind=kr), intent(out)   :: jnet(jmax+1)

        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: merr
        real(kind=kr)                :: psi
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: s(:)

      ! pre-compute coeffs

        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        c1=0.0_kr
        c2=0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            c1(j,m)=0.5_kr*tau
            c2(j,m)=0.5_kr*h/mu(m)
          enddo
        enddo

      ! solve problem

        allocate(s(jmax))
        do j=1,jmax
          phi(j)=1.0_kr
          s(j)  =0.5_kr*(sigs(j)*phi(j)+q(j))
        enddo

        psi_in=0.0_kr
        psi_bc=0.0_kr

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0_kr
          jnet=0.0_kr
          !$omp parallel do private(j, psi_in, psi) reduction(+: phi, jnet)
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            if (bc(1) == 0) psi_in=0.0_kr
            do j=1,jmax
              jnet(j)=jnet(j)+psi_in*mu(m)*w(m)
              psi    =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
              phi(j) =phi(j)+psi*w(m)
              psi_in =2.0_kr*psi-psi_in
            enddo
            jnet(jmax+1)=jnet(jmax+1)+psi_in*mu(m)*w(m)
            if (bc(2) == 0) psi_in=0.0_kr
            jnet(jmax+1)=jnet(jmax+1)-psi_in*mu(m)*w(m)
            do j=jmax,1,-1
              psi    =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
              phi(j) =phi(j)+psi*w(m)
              psi_in =2.0_kr*psi-psi_in
              jnet(j)=jnet(j)-psi_in*mu(m)*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo
          !$omp end parallel do

          merr=0.0_kr
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/phi(j)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s(j)=0.5_kr*(sigs(j)*phi(j)+q(j))
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

      subroutine solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q(jmax)
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt(jmax)
        real(kind=kr), intent(in)    :: sigs(jmax)
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(out)   :: phi (jmax)
        real(kind=kr), intent(out)   :: jnet(jmax+1)

        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: tau3
        real(kind=kr)                :: tau5
        real(kind=kr)                :: tau7
        real(kind=kr)                :: merr
        real(kind=kr)                :: psi
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: s(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            if (tau < 0.01_kr) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0_kr-tau3/360.0_kr+tau5/15120.0_kr-tau7/604800.0_kr
            else
              alpha(j,m)=1.0_kr/tanh(tau/2.0_kr)-2.0_kr/tau
            endif
            c1(j,m)=0.5_kr*tau*(1.0_kr+alpha(j,m))
            c2(j,m)=0.5_kr*h  *(1.0_kr+alpha(j,m))/mu(m)
          enddo
        enddo

      ! solve problem

        allocate(s(jmax))
        do j=1,jmax
          phi(j)=1.0_kr
          s(j)  =0.5_kr*(sigs(j)*phi(j)+q(j))
        enddo

        psi_in=0.0_kr
        psi_bc=0.0_kr

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0_kr
          jnet=0.0_kr
          !$omp parallel do private(j, psi_in, psi) reduction(+: phi, jnet)
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            if (bc(1) == 0) psi_in=0.0_kr
            do j=1,jmax
              jnet(j)=jnet(j)+psi_in*mu(m)*w(m)
              psi    =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
              phi(j) =phi(j)+psi*w(m)
              psi_in =(2.0_kr*psi-(1.0_kr-alpha(j,m))*psi_in)/(1.0_kr+alpha(j,m))
            enddo
            jnet(jmax+1)=jnet(jmax+1)+psi_in*mu(m)*w(m)
            if (bc(2) == 0) psi_in=0.0_kr
            jnet(jmax+1)=jnet(jmax+1)-psi_in*mu(m)*w(m)
            do j=jmax,1,-1
              psi    =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
              phi(j) =phi(j)+psi*w(m)
              psi_in =(2.0_kr*psi-(1.0_kr-alpha(j,m))*psi_in)/(1.0_kr+alpha(j,m))
              jnet(j)=jnet(j)-psi_in*mu(m)*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo
          !$omp end parallel do

          merr=0.0_kr
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/phi(j)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s(j)=0.5_kr*(sigs(j)*phi(j)+q(j))
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

      subroutine solve_ld(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi,jnet)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q(jmax)
        real(kind=kr), intent(in)    :: ql(jmax)
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt(jmax)
        real(kind=kr), intent(in)    :: sigs(jmax)
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(out)   :: phi (jmax)
        real(kind=kr), intent(out)   :: jnet(jmax+1)

        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: merr
        real(kind=kr)                :: psi
        real(kind=kr)                :: psil
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_out
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: sl(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            alpha(j,m)=1.0_kr/(1.0_kr+6.0_kr/tau)
            c1(j,m)   =       (2.0_kr/tau+alpha(j,m)-1.0_kr)
            c2(j,m)   =1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        do j=1,jmax
          phi(j) =1.0_kr
          phil(j)=0.0_kr
          s (j)  =0.5_kr*(sigs(j)*phi(j)+q(j))
          sl(j)  =0.0_kr
        enddo

        psi_in=0.0_kr
        psi_bc=0.0_kr

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0_kr
          phil=0.0_kr
          jnet=0.0_kr
          !$omp parallel do private(j, psi_in, psi_out, psi, psil) reduction(+: phi, phil, jnet)
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            if (bc(1) == 0) psi_in=0.0_kr
            do j=1,jmax
              jnet(j)=jnet(j)+psi_in*mu(m)*w(m)
              psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt(j)
              psil   =psi_out-psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            jnet(jmax+1)=jnet(jmax+1)+psi_in*mu(m)*w(m)
            if (bc(2) == 0) psi_in=0.0_kr
            jnet(jmax+1)=jnet(jmax+1)-psi_in*mu(m)*w(m)
            do j=jmax,1,-1
              psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt(j)
              psil   =psi-psi_out
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
              jnet(j)=jnet(j)-psi_in*mu(m)*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo
          !$omp end parallel do

          merr=0.0_kr
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/phi(j)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5_kr*(sigs(j)*phi (j)+q (j))
            sl(j)=0.5_kr*(sigs(j)*phil(j)+ql(j))
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

      subroutine solve_lc(n,jmax,kmax,h,q,ql,eps,sigt,sigs,mu,w,bc,phi,jnet)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q(jmax)
        real(kind=kr), intent(in)    :: ql(jmax)
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt(jmax)
        real(kind=kr), intent(in)    :: sigs(jmax)
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(out)   :: phi (jmax)
        real(kind=kr), intent(out)   :: jnet(jmax+1)

        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: tau3
        real(kind=kr)                :: tau5
        real(kind=kr)                :: tau7
        real(kind=kr)                :: merr
        real(kind=kr)                :: psi
        real(kind=kr)                :: psil
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_out
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: rbeta(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: sl(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(rbeta(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        rbeta=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt(j)*h/mu(m)
            if (tau < 0.01_kr) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0_kr-tau3/360.0_kr+tau5/15120.0_kr-tau7/604800.0_kr
            else
              alpha(j,m)=1.0_kr/tanh(tau/2.0_kr)-2.0_kr/tau
            endif
            rbeta(j,m)=1.0_kr/alpha(j,m)-6.0_kr/tau
            c1(j,m)=       (2.0_kr/tau+alpha(j,m)-1.0_kr)
            c2(j,m)=1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        do j=1,jmax
          phi(j) =1.0_kr
          phil(j)=0.0_kr
          s (j)  =0.5_kr*(sigs(j)*phi(j)+q(j))
          sl(j)  =0.0_kr
        enddo

        psi_in=0.0_kr
        psi_bc=0.0_kr

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          phi =0.0_kr
          phil=0.0_kr
          jnet=0.0_kr
          !$omp parallel do private(j, psi_in, psi_out, psi, psil) reduction(+: phi, phil, jnet)
          do m=1,n/2
            psi_in=psi_bc(m) ! left specular bc
            if (bc(1) == 0) psi_in=0.0_kr
            do j=1,jmax
              jnet(j)=jnet(j)+psi_in*mu(m)*w(m)
              psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt(j)
              psil   =((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr-rbeta(j,m)*psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
            enddo
            jnet(jmax+1)=jnet(jmax+1)+psi_in*mu(m)*w(m)
            if (bc(2) == 0) psi_in=0.0_kr
            jnet(jmax+1)=jnet(jmax+1)-psi_in*mu(m)*w(m)
            do j=jmax,1,-1
              psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt(j)+c1(j,m)*psi_in)
              psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt(j)
              psil   =-((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr+rbeta(j,m)*psi
              psi_in =psi_out
              phi(j) =phi(j) +psi*w(m)
              phil(j)=phil(j)+psil*w(m)
              jnet(j)=jnet(j)-psi_in*mu(m)*w(m)
            enddo
            psi_bc(m)=psi_in
          enddo
          !$omp end parallel do

          merr=0.0_kr
          do j=1,jmax
            merr=max(merr,abs((phi(j)-phio(j))/phi(j)))
          enddo
          if (merr < eps) exit

          do j=1,jmax
            s (j)=0.5_kr*(sigs(j)*phi (j)+q (j))
            sl(j)=0.5_kr*(sigs(j)*phil(j)+ql(j))
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

    end program main
