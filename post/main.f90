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

      if (.true.) then
        call post_kss
      else
        call post_dean
      endif

      contains

      subroutine post_kss

        use global

        implicit none

        integer(4)                   :: src
        integer(4)                   :: prb

        do src=1,2
          do prb=1,3
            call postp(src,prb)
          enddo
        enddo
        do src=3,4
          prb=1
          call postp(src,prb) 
        enddo

      end subroutine post_kss

      subroutine post_dean

        use global

        implicit none

        call postp(1,4) 

      end subroutine post_dean

      subroutine postp(src,prb) 

        use global

        implicit none

        integer(4),    intent(in)    :: src
        integer(4),    intent(in)    :: prb

        integer(4)                   :: j
        integer(4)                   :: jmax
        integer(4)                   :: jmaxr
        integer(4)                   :: sol
        integer(4)                   :: nx
        integer(4)                   :: xn
        integer(4)                   :: xnmax
        integer(4),    parameter     :: xnr=18
        real(kind=kr)                :: xmax
        real(kind=kr)                :: xmsh
        real(kind=kr)                :: xref
        real(kind=kr)                :: hr
        real(kind=kr)                :: h(xnr)
        real(kind=kr)                :: l2(xnr,0:3)
        real(kind=kr), allocatable   :: phi(:)
        real(kind=kr), allocatable   :: phir(:)
        character(1)                 :: prbopt
        character(1)                 :: srcopt
        character(1)                 :: solopt
        character(2)                 :: xnopt
        character(13)                :: output

      ! reference solution

        write(prbopt,'(i1)') prb
        write(srcopt,'(i1)') src
        write(xnopt ,'(i2)') xnr
        output='p'//prbopt//'-'//srcopt//'-0-'//xnopt//'.dat'
        if (prb == 4) output='p'//prbopt//'-'//srcopt//'-3-'//xnopt//'.dat'

        xmax=40.0_kr
        xref=10.0_kr
        if (prb == 4) then
          xmax=8.0_kr
          xref=2.0_kr
        endif

        xn   =xnr
        nx   =5*2**(xn-1)
        hr   =xref/nx
        jmaxr=nint(xmax/hr)
        allocate(phir(jmaxr))
        phir=0.0_kr
        open(unit=1,file=output,action='read',status='unknown')
          do j=1,jmaxr
            read(1,*) xmsh,phir(j)
            phir(j)=phir(j)
          enddo
        close(1)

      ! all other solutions

        h=0.0_kr
        l2=0.0_kr
        do sol=0,3
          xnmax=12
          if (sol == 2) then
            xnmax=9
          elseif (sol == 3) then
            xnmax=7
          endif
          if (prb == 4) xnmax=14
          do xn=1,xnmax
            write(prbopt,'(i1)') prb
            write(solopt,'(i1)') sol
            write(srcopt,'(i1)') src
            write(xnopt, '(i2)') xn
            nx   =5*2**(xn-1)
            h(xn)=xref/nx
            jmax=nint(xmax/h(xn))
            allocate(phi(jmax))
            phi=0.0_kr
            output='p'//prbopt//'-'//srcopt//'-'//solopt//'-'//trim(adjustl(xnopt))//'.dat'
            open(unit=1,file=output,action='read',status='unknown')
              do j=1,jmax
                read(1,*) xmsh,phi(j)
                phi(j)=phi(j)
              enddo
              if (prb == 4) then
                call geterr(jmax,jmaxr,phi,phir,h(xn),hr,l2(xn,sol))
              else
                call geterr2(jmax,jmaxr,xmax,phi,phir,h(xn),hr,l2(xn,sol))
              endif
            close(1)
            deallocate(phi)
          enddo
        enddo

        deallocate(phir)

      ! print results to text file

        call results(xnr,src,prb,h,l2)

      end subroutine postp

      subroutine geterr(jmax,jmaxr,phi,phir,h,hr,l2)

        use global

        implicit none

        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: jmaxr
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: hr
        real(kind=kr), intent(inout) :: l2
        real(kind=kr), intent(in)    :: phi(jmax)
        real(kind=kr), intent(in)    :: phir(jmaxr)

        integer(4)                   :: j
        integer(4)                   :: jj
        integer(4)                   :: jjj
        integer(4)                   :: jx
        real(kind=kr), allocatable   :: phirr(:)

      ! integrate reference solution

        jx=jmaxr/jmax
        allocate(phirr(jmax))
        phirr=0.0_kr

        l2=0.0_kr
        do jj=1,jmax
          do j=1,jx
            jjj=j+(jj-1)*jx
            phirr(jj)=phirr(jj)+phir(jjj)*hr
          enddo
          phirr(jj)=phirr(jj)/h
          l2=l2+h*(phi(jj)-phirr(jj))**2
        enddo
        l2=sqrt(l2)

        deallocate(phirr)

      end subroutine geterr

      subroutine geterr2(jmax,jmaxr,xmax,phi,phir,h,hr,l2)

        use global

        implicit none

        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: jmaxr
        real(kind=kr), intent(in)    :: xmax
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: hr
        real(kind=kr), intent(inout) :: l2
        real(kind=kr), intent(in)    :: phi(jmax)
        real(kind=kr), intent(in)    :: phir(jmaxr)

        integer(4)                   :: j
        integer(4)                   :: jj
        integer(4)                   :: jjj
        integer(4)                   :: jmax2
        integer(4)                   :: jx
        real(kind=kr), allocatable   :: phi2(:)
        real(kind=kr), allocatable   :: phir2(:)

      ! compare over 2.0 cm mesh

        jmax2=nint(xmax/2.0_kr)
        allocate(phi2(jmax2))
        allocate(phir2(jmax2))
        phi2=0.0_kr
        phir2=0.0_kr

      ! integrate approximate solution

        jx=jmax/jmax2
        do jj=1,jmax2
          do j=1,jx
            jjj=j+(jj-1)*jx
            phi2(jj)=phi2(jj)+phi(jjj)*h
          enddo
          phi2(jj)=phi2(jj)/2.0_kr
        enddo

      ! integrate approx. and reference over 2.0 cm mesh

        jx=jmaxr/jmax2
        do jj=1,jmax2
          do j=1,jx
            jjj=j+(jj-1)*jx
            phir2(jj)=phir2(jj)+phir(jjj)*hr
          enddo
          phir2(jj)=phir2(jj)/2.0_kr
        enddo

      ! compute l2 error norm

        l2=0.0_kr
        do jj=1,jmax2
          l2=l2+(phi2(jj)-phir2(jj))**2
        enddo
        l2=sqrt(l2)

        deallocate(phi2)
        deallocate(phir2)

      end subroutine geterr2

      subroutine results(xnr,src,prb,h,l2)

        use global

        implicit none

        integer(4),    intent(in)    :: xnr
        integer(4),    intent(in)    :: src
        integer(4),    intent(in)    :: prb
        real(kind=kr), intent(in)    :: h(xnr)
        real(kind=kr), intent(in)    :: l2(xnr,0:3)

        integer(4)                   :: sol
        integer(4)                   :: xn
        integer(4)                   :: p
        real(kind=kr)                :: cnst
        real(kind=kr)                :: cord(xnr,0:3)
        character(1)                 :: prbopt
        character(1)                 :: srcopt
        character(6)                 :: output

        cord=0.0_kr
        do sol=0,3
          xn=12
          if (sol == 2) then
            xn=9
          elseif (sol == 3) then
            xn=7
          endif
          if (prb == 4) xn=14
          p=2
          if (sol == 2) p=3
          if (sol == 3) p=4
          cnst=l2(xn,sol)/(h(xn)**p)
          do xn=1,xn
            cord(xn,sol)=cnst*(h(xn)**p)
          enddo
        enddo

        write(prbopt,'(i1)') prb
        write(srcopt,'(i1)') src
        output='p'//prbopt//'-'//srcopt//'.p'
        open(unit=1,file=output,action='write',status='unknown')
        do xn=1,xnr
          write(1,'(9(es25.16))') h(xn),(l2(xn,sol),sol=0,3),(cord(xn,sol),sol=0,3)
        enddo
        close(1)

      end subroutine results

    end program main
