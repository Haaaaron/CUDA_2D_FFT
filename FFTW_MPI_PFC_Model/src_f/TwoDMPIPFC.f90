! To transform a one-dimensional complex array
!             double complex in, out
!             dimension in(N), out(N)
!             integer*8 plan
!
!             call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
!             call dfftw_execute_dft(plan, in, out)
!             call dfftw_destroy_plan(plan)
!             double precision in
! To transform a one-dimensional real array in Fortran, you might do:
!             dimension in(N)
!             double complex out
!             dimension out(N/2 + 1)
!             integer*8 plan
!
!             call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
!             call dfftw_execute_dft_r2c(plan, in, out)
!             call dfftw_destroy_plan(plan)
! ftn -fastsse Mipa=fast -O2 -r8 -I/home/u22/cva/louhi/usr/local/include/ TwoDMPIPFC.f /home/u22/cva/louhi/usr/local/lib/libfftw3.a
! mpiifort -r8 -O2 -lfftw3 TwoDMPIPFC.f
! Load the modules !!!!!
! mpif90 -fdefault-real-8 -fdefault-double-8 -O2 TwoDMPIPFC.f  -lfftw3


program si
    implicit none
    include "mpif.h"
    include "fftw3.f"
    integer,parameter :: lx=256,lxhp=lx/2+1,lhx=lx/2
    integer,parameter :: ly=128,lyhp=ly/2+1,lhy=ly/2
    integer,parameter :: nend=100
    !integer,parameter :: nend=1 !test
    double precision, parameter :: pi=3.14159265358979
    double precision, parameter :: sq3=1.732050807568877
    double precision, parameter :: at=7.255197456936871
    double precision, parameter :: qqt=2*pi/at
    double precision, parameter :: dt=0.125,VV=1.d0/(lx*ly)
    !double precision, parameter :: dy=(at*sq3/2.)/8,dx=dy
    double precision, parameter :: dx=at/8.0,dy=(at*sq3/2.)/8.0
    double precision, parameter :: r=-0.25,pm=-0.25
    double precision :: ttime(2),rnd
    double precision :: ttimesfsfftw,ttimempicom,ttimesssfftw
    integer :: j,i,kk, provided,j_len
    integer*8 planpsi,iplanpsi,plannt,planbff,iplanbff
    integer*8 planpsic,iplanpsic
    integer*8 plany1,plany2
    integer*8 iplany1,iplany2
    integer :: status(MPI_STATUS_SIZE),ierr
    integer :: myid,numprocs
    integer :: jminmyid,jmaxmyid,iminmyid,imaxmyid
    double precision :: qq,kkx,kky
    integer :: rcounts,scounts
    integer,allocatable :: rdispls(:),sdispls(:)
    integer,allocatable :: imin(:),imax(:)
    integer,allocatable :: jmin(:),jmax(:)
    integer :: lxr
    integer :: lyr,np,n,sumscounts,sumrcounts
    double complex,allocatable:: ckpsi(:,:),ntc(:,:)
    double precision, allocatable :: cpsi(:,:),cntc(:,:),cene(:,:)
    double precision :: qcpsizero(lyhp),qcpsilxhp(lyhp)
    double precision :: qcntczero(lyhp),qcntclxhp(lyhp)
    double precision :: qcenezero(lyhp),qcenelxhp(lyhp)
    double precision, allocatable :: psi(:,:),nt(:,:),tmp(:,:)
    double complex :: cppp(lxhp),ckpppin(ly),ckpppout(ly)
    double complex :: cplxhp(lyhp),cpzero(lyhp)
    double precision :: ppp(1:lx),pzero(ly),plxhp(ly)
    double complex :: cpsizero(lyhp),cpsilxhp(lyhp)
    double complex :: cntzero(lyhp),cntlxhp(lyhp)
    real*8  :: ene, totene

    !call MPI_INIT(ierr)




    call MPI_Init_thread( MPI_THREAD_FUNNELED, provided, ierr )
    write(*,*) provided
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    print *, "Process ", myid, " of ", numprocs, " is alive"
    call fft_init()

    ! print*,"energy.lxhp"
    ! print'(E25.17)',qcenelxhp
    ! print*,""
    ! print*,"energy.zero"
    ! print'(E25.17)',qcenezero
    ! print*,""
    ! print*,"psi.lxhp"
    ! print'(E25.17)',qcpsilxhp
    ! print*,""
    ! print*,"psi.zero"
    ! print'(E25.17)',qcpsizero
    ! print*,""
    ! print*,"ntc.lxhp"
    ! print'(E25.17)',qcntclxhp
    ! print*,""
    ! print*,"ntc.zero"
    ! print'(E25.17)',qcntczero
    ! print*,""
    ! print*,"psi.rest"
    ! print'(E25.17)',cpsi
    ! print*,""
    ! print*,"ntc.rest"
    ! print'(E25.17)',cntc
    ! print*,""
    ! print*,"energy.rest"
    ! print'(E25.17)',cene
    ! print*,""
    ! stop
    psi=0
    !print*,(0*2-1)*0.15-0.13*(dcos(dsqrt(3.d0)*i*dx/2.d0)*cos(j*dy/2.d0)-0.5*dcos(j*dy))+pm
    do j=jminmyid,jmaxmyid
        do i=1,lx
            psi(i,j)=-0.13*(dcos(dsqrt(3.d0)*i*dx/2.d0)*cos(j*dy/2.d0)-0.5*dcos(j*dy))+pm
            !print*,-0.13*(dcos(dsqrt(3.d0)*i*dx/2.d0)*cos(j*dy/2.d0)-0.5*dcos(j*dy))+pm
            !if (j == 3 .and. i == 3) psi(i,j) = 1
            !psi(i,j)= cos((j)*2*pi/2)*0+cos((i)*2*pi/2)!cos((j)*2*pi/2)+cos((j)*2*pi/4)+cos((j)*2*pi/8)+cos((j)*2*pi/16)+cos((j)*2*pi/32)
        enddo
    enddo

    rnd=sum(psi)/(lx*lyr)-pm

    psi=psi-rnd


    ttimesfsfftw=0.0
    ttimempicom=0.0
    ttimesssfftw=0.0

    call cpu_time(ttime(1))
    !call energy()
    if(myid.eq.0) then
       write(*,*) totene,'Init'
    endif
    write(*,*) ene,myid,'Init'
    !call update()

    !Main loop
    do n=1,nend
        call energy()
        if(myid.eq.0) then
            write(*,*) totene, psi(1,1)
        endif
        call update()
        !call tests()

        if(mod(n,500).eq.0) then
            call energy()

            if(myid.eq.0) then
                write(*,*) totene, n/50
            endif
        endif
    enddo

    call energy()

    if(myid.eq.0) then
       write(*,*) totene,'Fine'
    endif
    ! write(*,*) ene,myid,'Fine'

    call cpu_time(ttime(2))

    write(*,*) 'Total time', ttime(2)-ttime(1), myid, ttimesfsfftw, ttimempicom, ttimesssfftw

    call MPI_FINALIZE(ierr)

contains

    subroutine energy()

        double precision :: bffr(1:lx,jminmyid:jmaxmyid)
        double complex :: bffc(1:ly,iminmyid:imaxmyid)

        bffc=0
        call forward_fft(psi,bffc)

        select case(myid)
        case(0)
            cpzero=qcenezero*cpzero
            cplxhp=qcenelxhp*cplxhp
            do i=2,imaxmyid
                do j=1,ly
                    bffc(j,i)=cene(j,i)*bffc(j,i)
                enddo
            enddo
        case default

            bffc=cene*bffc

        end select
        !if (myid==1) print'(2E16.7)',bffc

        call backward_fft(bffc,bffr)

        ene=sum((1.d0/4.d0)*psi**4+(1.d0/2.d0)*psi*bffr)*VV

        call MPI_REDUCE(ene, totene, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        return

    end subroutine energy

    subroutine update()
        
        nt=psi**3
        call forward_fft(nt,ntc)
        ntc=ntc
        cntzero=cpzero
        cntlxhp=cplxhp

        call forward_fft(psi,ckpsi)     !in,out

        select case(myid)
        case(0)

            cpzero=qcpsizero*cpzero+qcntczero*cntzero
            cplxhp=qcpsilxhp*cplxhp+qcntclxhp*cntlxhp
            do i=2,imaxmyid
                do j=1,ly
                    ckpsi(j,i)=cpsi(j,i)*ckpsi(j,i)+cntc(j,i)*ntc(j,i)
                    !print'(E16.5,I5)',cpsi(j,i),i

                enddo
            enddo

        case default
            ckpsi=cpsi*ckpsi+cntc*ntc
        end select

        call backward_fft(ckpsi,psi)

        return

    end subroutine update

    subroutine coeff()
        double precision :: ccc

        select case(myid)
        case(0)
            do j=1,lyhp
                kkx=0
                kky=(2*pi*(j-1)/(ly*dy))
                qq=-(kkx**2+kky**2)
                ccc=r+1+2*qq+qq*qq
                qcenezero(j)=VV*ccc
                qcpsizero(j)=1.d0*VV/(1-dt*qq*ccc)
                qcntczero(j)=dt*qq*VV/(1-dt*qq*ccc)
            enddo

            do j=1,lyhp
                kkx=pi/dx
                kky=(2*pi*(j-1)/(ly*dy))
                qq=-(kkx**2+kky**2)
                ccc=r+1+2*qq+qq*qq
                qcenelxhp(j)=VV*ccc
                qcpsilxhp(j)=1.d0*VV/(1-dt*qq*ccc)
                qcntclxhp(j)=dt*qq*VV/(1-dt*qq*ccc)
            enddo

            do i=2,imaxmyid
                kkx=2*pi*(i-1)/(lx*dx)

                do j=1,ly
                    if(j.le.ly/2+1) then
                        kky=2*pi*(j-1)/(ly*dy)
                    else
                        kky=2*pi*(j-1-ly)/(ly*dy)

                    endif
                    qq=-(kkx**2+kky**2)

                    ccc=r+1+2*qq+qq*qq
                    cene(j,i)=VV*ccc
                    cpsi(j,i)=1.d0*VV/(1-dt*qq*ccc)
                    !print*,cpsi(j,i),i
                    cntc(j,i)=qq*dt*VV/(1-dt*qq*ccc)
                enddo
            enddo
            !stop
        case default
            do i=iminmyid,imaxmyid
                kkx=2*pi*(i-1)/(lx*dx)
                do j=1,ly
                    if(j.le.ly/2+1) then
                        kky=2*pi*(j-1)/(ly*dy)
                    else
                        kky=2*pi*(j-1-ly)/(ly*dy)
                    endif
                    qq=-(kkx*kkx+kky*kky)
                    ccc=r+1+2*qq+qq*qq
                    cene(j,i)=VV*ccc
                    cpsi(j,i)=1.d0*VV/(1-dt*qq*ccc)
                    !print*,cpsi(j,i),i

                    cntc(j,i)=qq*dt*VV/(1-dt*qq*ccc)
                enddo
            enddo
            !stop
        end select

        ! qq=-kkx*kkx+kky*kky
        ! ccc=r+1+2*qq+qq*qq
        ! cene(j,i)=VV*ccc
        ! cpsi(j,i)=VV/(1.d0-dt*qq*ccc)
        ! cntc(j,i)=VV*qq*dt/(1.d0-dt*ccc*qq)

        ! ckpsi(j,i)=cpsi(j,i)*ckpsi(j,i)+cntc(j,i)*ntc(j,i)
        return
    end subroutine

    subroutine fft_init()

        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

        allocate(sdispls(0:numprocs-1))
        allocate(rdispls(0:numprocs-1))
        allocate(imin(0:numprocs-1))
        allocate(imax(0:numprocs-1))
        allocate(jmin(0:numprocs-1))
        allocate(jmax(0:numprocs-1))

        lyr=ly/numprocs

        jmin(0)=1
        jmax(0)=lyr

        do np=1,numprocs-1
            jmin(np)=jmax(np-1)+1
            jmax(np)=jmax(np-1)+lyr
        enddo

        lxr=lhx/numprocs
        imin(0)=1
        imax(0)=lxr

        do np=1,numprocs-1
            imin(np)=imax(np-1)+1
            imax(np)=imax(np-1)+lxr
        enddo

        ! rcounts (0:numprocs-1) specifies the number of elements to be received from each process for forward transform
        rcounts=lxr*lyr

        ! scounts (0:numprocs-1) specifies the number of elemnts to be sent to each proc
        scounts=lxr*lyr

        ! sdipls specifies the displacement from which to take the the outgoing data for process j

        sdispls=0

        do np=1,numprocs-1
            sdispls(np)=sdispls(np-1)+scounts
        enddo

        ! rdipls specifies the at which to place the incoming data from each process

        rdispls=0
        do np=1,numprocs-1
            rdispls(np)=rdispls(np-1)+rcounts
        enddo

        sumscounts=scounts*numprocs
        sumrcounts=rcounts*numprocs

        call dfftw_plan_dft_r2c_1d(planpsi, lx, ppp, cppp, FFTW_ESTIMATE)  ! create the plan for forward ft of psi
        call dfftw_plan_dft_r2c_1d(plany1, ly, pzero, cpzero, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_1d(plany2, ly, plxhp, cplxhp, FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(planpsic, ly, ckpppin, ckpppout, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(iplanpsic, ly, ckpppout, ckpppin, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_1d(iplanpsi,lx,cppp,ppp, FFTW_ESTIMATE)  ! create the plan for forward ft of psi
        call dfftw_plan_dft_c2r_1d(iplany1,ly,cpzero,pzero, FFTW_ESTIMATE)  ! create the plan for forward ft of
        call dfftw_plan_dft_c2r_1d(iplany2,ly,cplxhp,plxhp, FFTW_ESTIMATE) ! create the plan for forward ft of psi

        call MPI_BARRIER( MPI_COMM_WORLD, ierr )

        jminmyid=jmin(myid)
        jmaxmyid=jmax(myid)
        iminmyid=imin(myid)
        imaxmyid=imax(myid)

        allocate(psi(1:lx,jminmyid:jmaxmyid))
        allocate(tmp(1:lx,jminmyid:jmaxmyid))
        allocate(nt(1:lx,jminmyid:jmaxmyid))
        allocate(ntc(1:ly,iminmyid:imaxmyid))
        allocate(ckpsi(1:ly,iminmyid:imaxmyid))
        allocate(cpsi(1:ly,iminmyid:imaxmyid))
        allocate(cntc(1:ly,iminmyid:imaxmyid))
        allocate(cene(1:ly,iminmyid:imaxmyid))
        !-----
        !test
        ! print *, jminmyid, jmaxmyid, iminmyid, imaxmyid
        ! print *, size(psi), shape(psi)
        ! j_len = jmaxmyid - jminmyid
        ! kk=0
        ! do i=1,lx
        !     do j=jminmyid,jmaxmyid
        !         psi(j,i) = kk
        !         kk=kk+1
        !     end do
        ! end do
        ! if (myid==0) then
        !     do i=1,lx
        !         do j=jminmyid,jmaxmyid
        !             print*,psi(j,i),myid,j_len*lx
        !         end do
        !     end do
        ! end if
        !print*,psi
        !----------
        call coeff()

        return
    end subroutine fft_init

    subroutine forward_fft(in,out)
        !------------------------------------------------------------------------------------------------------------------
        ! Forward Fourier transform of a real matrix in to a matrix out. The matrix out is transpose
        ! It uses variables defined in the subroutine fft_init. The size of the matrix must be defined in the main program.
        !------------------------------------------------------------------------------------------------------------------
        double complex :: ckintpsitpose(1:sumscounts),cktmp
        double complex :: work(1:sumrcounts)
        double precision :: in(1:lx,jminmyid:jmaxmyid)
        double complex :: out(1:ly,iminmyid:imaxmyid)
        double precision:: tt(2)
        integer :: ind

        call cpu_time(tt(1))
        do j=jminmyid,jmaxmyid
            ppp(1:lx)=in(1:lx,j)
            call dfftw_execute_dft_r2c(planpsi,ppp,cppp)
            cppp(1)=dcmplx(dreal(cppp(1)),dreal(cppp(lxhp)))
            !print*,cppp(1),j
            do i=1,lhx
                ckintpsitpose((i-1)*lyr+j-jminmyid+1)=cppp(i)
                !if (myid==1) print*,cppp(i),i,j,(i-1)*lyr+j-jminmyid+1
                !ckintpsitpose((i-1)*lyr+j-jminmyid+1)=cmplx(i,j)
            enddo
        enddo

        call cpu_time(tt(2))
        ttimesfsfftw=ttimesfsfftw+tt(2)-tt(1)

        call cpu_time(tt(1))
        call MPI_alltoall(ckintpsitpose, scounts, MPI_DOUBLE_COMPLEX, work, rcounts, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call cpu_time(tt(2))
        ttimempicom=ttimempicom+tt(2)-tt(1)

        call cpu_time(tt(1))
        select case(myid)
        case(0)
            do np=0,numprocs-1
                do j=jmin(np),jmax(np)
                    cktmp=work(rdispls(np)+j-jmin(np)+1)
                    pzero(j)=dreal(cktmp)
                    plxhp(j)=dimag(cktmp)
                enddo
            enddo

            call dfftw_execute_dft_r2c(plany1,pzero,cpzero)
            call dfftw_execute_dft_r2c(plany2,plxhp,cplxhp)

            do i=2,imaxmyid

                do np=0,numprocs-1
                    do j=jmin(np),jmax(np)
                        ckpppin(j)=work(rdispls(np)+(i-1)*lyr+j-jmin(np)+1)
                    enddo
                enddo

                call dfftw_execute_dft(planpsic, ckpppin, ckpppout)

                out(1:ly,i)=ckpppout(1:ly)
            enddo

        case default
            do i=iminmyid,imaxmyid
                do np=0,numprocs-1
                    do j=jmin(np),jmax(np)
                        ckpppin(j)=work(rdispls(np)+(i-iminmyid)*lyr+j-jmin(np)+1)
                    enddo
                enddo
                ! do j = 1, ly
                !     print*,ckpppin(j),j,i,ly
                ! end do
                call dfftw_execute_dft(planpsic, ckpppin, ckpppout)
                !print'(2F13.8)',ckpppout
                out(1:ly,i)=ckpppout(1:ly)
            enddo

        end select

        call cpu_time(tt(2))
        ttimesssfftw=ttimesssfftw+tt(2)-tt(1)

        return
    end subroutine forward_fft

    subroutine backward_fft(out,in)
        !-----------------------------------------------------------------------
        ! Inverse Fourier transform. Matrix out is transposed,complex, hermitian
        !-----------------------------------------------------------------------
        double complex :: ckintpsitpose(1:sumscounts),cktmp
        double complex :: work(1:sumrcounts)
        double precision :: in(1:lx,jminmyid:jmaxmyid)
        double complex :: out(1:ly,iminmyid:imaxmyid)
        double precision:: tt(2)

        call cpu_time(tt(1))
        select case(myid)
        case (0)

            call dfftw_execute_dft_c2r(iplany1,cpzero,pzero)
            call dfftw_execute_dft_c2r(iplany2,cplxhp,plxhp)

            do j=1,ly
                work((j-1)*lxr+1)=dcmplx(pzero(j),plxhp(j))

            enddo
            
            do i=2,imaxmyid
                ckpppout(1:ly)=out(1:ly,i)

                call dfftw_execute_dft(iplanpsic, ckpppout, ckpppin)

                do j=1,ly
                    work((j-1)*lxr+i)=ckpppin(j)

                enddo
            enddo
        case default
            do i=iminmyid,imaxmyid
                ckpppout(1:ly)=out(1:ly,i)

                call dfftw_execute_dft(iplanpsic, ckpppout, ckpppin)

                do j=1,ly
                    work((j-1)*lxr+i-iminmyid+1)=ckpppin(j)
                enddo
            enddo
        end select
        call cpu_time(tt(2))
        ttimesfsfftw=ttimesfsfftw+tt(2)-tt(1)

        call cpu_time(tt(1))
        call MPI_alltoall(work,rcounts,MPI_DOUBLE_COMPLEX,ckintpsitpose,scounts,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
        call cpu_time(tt(2))
        ttimempicom=ttimempicom+tt(2)-tt(1)

        call cpu_time(tt(1))
        do j = jminmyid,jmaxmyid
            do np = 0,numprocs-1
                do i = imin(np), imax(np)
                    cppp(i) = ckintpsitpose(sdispls(np)+(j-jminmyid)*lxr+i-imin(np)+1)
                enddo
            enddo

            cppp(lxhp)=dcmplx(dimag(cppp(1)),0.d0)
            cppp(1)=dcmplx(dreal(cppp(1)),0.d0)

            call dfftw_execute_dft_c2r(iplanpsi,cppp,ppp)
            in(1:lx,j)=ppp(1:lx)

        enddo
        call cpu_time(tt(2))
        ttimesssfftw=ttimesssfftw+tt(2)-tt(1)

        return
    end subroutine backward_fft

end program si
