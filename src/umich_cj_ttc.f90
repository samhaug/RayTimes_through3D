!******************************************************************
!umich_cj_tcc.f90
!Routine to calculate traveltime delays using ray theory
!Developed by Carlos Chaves and Jeroen Ritsema on 08/27/2016
!v1.0
!v.1.1 --> Reading station and event information
!v.1.2 --> Read any model in sph format / More information is
! necesssary in the input file.
! Don't forget to change the precison in .taup in latlon to
! be at least 3 digits
!******************************************************************

module model_par

    double precision, parameter :: R_EARTH = 6371000.d0
    double precision, parameter :: R_EARTH_KM = 6371.0
    double precision, parameter :: RCMB = 3480.d0
    double precision, parameter :: RICB = 1221.5
    double precision, parameter :: PI = 3.141592653589793d0
    double precision, parameter :: DEGREES_TO_RADIANS = PI/180.d0
    double precision, parameter :: RADIANS_TO_DEGREES = 180.0d0/PI

    double precision,dimension(:,:,:),allocatable :: V_dv_a,V_dv_b
    ! splines
    double precision,dimension(:),allocatable :: V_spknt
    double precision,dimension(:,:),allocatable :: V_qq0
    double precision,dimension(:,:,:),allocatable :: V_qq

end module model_par

!
!--------------------------------------------------------------------------------------------------
!
program umich_cj_tcc

    use model_par
    implicit none

    integer :: lines
    integer :: IIN=1,ier
    integer :: i,j,k
    integer :: NS !-->
    !integer, parameter :: NK = 20 !--> Number of splines (indeed are 21)
    integer :: NK !--> Number of splines (indeed are 21)
    integer :: region !--> Region where the stw105 1D model is calculated
    double precision :: dummy(4)
    double precision :: dtt !--> Travel time delay
    double precision :: ttheo !--> Theoretical travel time
    double precision :: ttomo !--> Travel time for the tomography model 
    double precision :: cospsi, delta
    double precision :: s1,s2,s3 !--> Slowness from the reference model
    double precision :: s1tomo,s2tomo,s3tomo !--> Slowness from a tomography model
    double precision :: r1,r2,r3 !--> radius in radians
    double precision :: theta1,theta2,theta3 !--> latitude in radians
    double precision :: ctheta1,ctheta2,ctheta3 !--> colatitude in radians
    double precision :: phi1,phi2,phi3 !--> longitude in radians
    double precision :: v1D !--> Velocity in km/s from a reference earth model
    double precision :: dv !--> Velocity perturbation from a tomography model (.sph file)
    double precision :: vv_1D_1,vh_1D_1,vv_1D_2,vh_1D_2,vv_1D_3,vh_1D_3 !-->
    ! Vertical and horizontal velocity (km/s)
    double precision :: vv_3D_1,vh_3D_1,vv_3D_2,vh_3D_2,vv_3D_3,vh_3D_3
    double precision,dimension(:,:),allocatable :: raypath
    real(kind=4) :: dvh,dvv !--> Vertical and horizontal velocity perturbation
    real(kind=4) :: xcolat,xlon,xrad !--> colatitude, longitude in degrees and
    ! radius in km
    character(len=300) :: frmt
    character(len=300) :: buffer,arg1,arg2,arg3,arg4
    character(len=100) :: INTEG !--> Method of numerical integration
    character(len=100) :: THREE_D_MODEL !--> 3D tomography model
    character(len=100) :: ONE_D_MODEL !--> 1D reference model
    character :: PHASE ! --> Compressional or shear-wave tomography model
    character(len=5) :: KFILE ! --> Kind of file: sph or hvd
    double precision ::  info(6) !--> Event and station information

    
    call getarg(1,buffer)
    read(buffer,*) lines
    call getarg(2,buffer)
    read(buffer,*) arg2
    call getarg(3,buffer)
    read(buffer,*) arg3
    !write(*,*) "ARG1: ",arg1
    !write(*,*) "ARG2: ",arg2
    !write(*,*) "ARG3: ",arg3
    dv = 0.d0
    ttheo = 0.d0
    ttomo = 0.d0
    v1D = 0.d0
    dvh = 0.0
    dvv = 0.0
    vv_1D_1=0.0
    vh_1D_1=0.0
    vv_1D_2=0.0
    vh_1D_2=0.0
    vv_1D_3=0.0
    vh_1D_3=0.0

    !open(unit=IIN,file='./evt',status='old',action='read',iostat=ier)
    !open(unit=IIN,file=arg1,status='old',action='read',iostat=ier)
    !if (ier /= 0) then
    !    write(*,*) 'Error opening file event information...', ier
    !    stop
    !else
    !    read(IIN,fmt=*,iostat=ier) (info(k),k=1,6)
    !endif
    !close(IIN)

    !open(unit=IIN,file='./noflines',status='old',action='read',iostat=ier)
    !open(unit=IIN,file=arg1,status='old',action='read',iostat=ier)
    !if (ier /= 0) then
    !    write(*,*) 'Error opening number of lines file...', ier
    !    stop
    !else
    !    read(IIN,fmt=*,iostat=ier) lines
    !endif
    !close(IIN)

    allocate(raypath(4,lines),stat=ier)
    !open(unit=IIN,file='./raypath',status='old',action='read',iostat=ier)
    open(unit=IIN,file=arg2,status='old',action='read',iostat=ier)
    if (ier /= 0) then
        write(*,*) 'Error opening file raypath', ier
        stop
    else
        do i=1,lines
            read(IIN,fmt=*,iostat=ier) (raypath(j,i),j=1,4)
        enddo
    endif
    close(IIN)

    !open(unit=IIN,file='./input',status='old',action='read',iostat=ier)
    open(unit=IIN,file=arg3,status='old',action='read',iostat=ier)
    if (ier /= 0) then
        write(*,*) 'Error opening file input', ier
        stop
    else
        read(IIN,fmt=*,iostat=ier) INTEG
        read(IIN,fmt=*,iostat=ier) ONE_D_MODEL
        read(IIN,fmt=*,iostat=ier) THREE_D_MODEL
        read(IIN,fmt=*,iostat=ier) PHASE
        read(IIN,fmt=*,iostat=ier) KFILE
    endif
    close(IIN)

    if (INTEG /= 'trapezoidal' .and. INTEG /= 'simpson') stop 'Unknown integration method. Check your input file'
    if (ONE_D_MODEL /= 'PREM' .and. ONE_D_MODEL /= 'STW105') stop 'Unknown 1D Model. Check your input file'
    if (PHASE /= 'P' .and. PHASE /= 'S') stop 'Unknown kind of Mantle tomography model'

    if (KFILE=='sph') then
        open(unit=IIN,file='/home/samhaug/RayTimes_through3D/models/'&
             //trim(THREE_D_MODEL)//'.sph',status='old',&
             action='read',iostat=ier)
        !open(unit=IIN,file=trim(THREE_D_MODEL)//'.sph',status='old',action='read',iostat=ier)
        if (ier /= 0) then
            write(*,*) 'Error opening "', trim(THREE_D_MODEL), '": ', ier
            stop 'Check your input file or the folder models because this routine did not find a tomography model.'
        endif
        read(IIN,*) (dummy(k),k=1,4)
        close(IIN)
        NS=int(dummy(1))
        NK=int(dummy(3)-4) !--> Crust is out
        call model_broadcast(NK,NS,THREE_D_MODEL)
    else if (KFILE=='hvd') then
        call model_1dref_broadcast()
        call model_s362ani_broadcast(THREE_D_MODEL)
    else
        stop 'Unknown tomography model format.'
    endif

!    if (THREE_D_MODEL  ==  'S40RTS' .or. THREE_D_MODEL  ==  'TX2015' .or. THREE_D_MODEL  == 'GYPSUM_S') then
!        NS=40
!        call model_broadcast(NK,NS,THREE_D_MODEL)
!        THREE_D_MODEL='sph'
!        PHASE = 'S'
!    else if (THREE_D_MODEL  ==  'S20RTS') then
!        NS=20
!        call model_broadcast(NK,NS,THREE_D_MODEL)
!        PHASE = 'S'
!        THREE_D_MODEL='sph'
!    else if (THREE_D_MODEL  ==  'SGLOBE_RANI') then
!        NS=35
!        call model_broadcast(NK,NS,THREE_D_MODEL)
!        THREE_D_MODEL='sph'
!        PHASE = 'S'
!    else if (THREE_D_MODEL  ==  'P12') then
!        NS=12
!        call model_broadcast(NK,NS,THREE_D_MODEL)
!        THREE_D_MODEL='sph'
!        PHASE = 'P'
!    else if (THREE_D_MODEL=='S362ANI_M' .or. THREE_D_MODEL=='S362WMANI_M' .and. ONE_D_MODEL == 'STW105') then
!        call model_1dref_broadcast()
!        call model_s362ani_broadcast(THREE_D_MODEL)
!        PHASE = 'S'
!    else
!        stop 'Unknown 3D model.'
!    endif

    do i=1,lines-1

        if(INTEG=='simpson') then
            r1=raypath(2,i)/R_EARTH_KM
            r2=raypath(2,i+1)/R_EARTH_KM
            r3=(r1+r2)/2.0
            theta1=raypath(3,i)*DEGREES_TO_RADIANS
            theta2=raypath(3,i+1)*DEGREES_TO_RADIANS
            theta3=(theta1+theta2)/2.0
            ctheta1=(90.0-raypath(3,i))*DEGREES_TO_RADIANS   !--> Colatitude
            ctheta2=(90.0-raypath(3,i+1))*DEGREES_TO_RADIANS !--> Colatitude
            ctheta3=(ctheta1+ctheta2)/2.0
            phi1=raypath(4,i)*DEGREES_TO_RADIANS
            phi2=raypath(4,i+1)*DEGREES_TO_RADIANS
            phi3=(phi1+phi2)/2.0
            cospsi=dsin(theta1)*dsin(theta2)+(dcos(theta1)*dcos(theta2)*dcos(phi2-phi1))
            delta=sqrt(r1**2+r2**2-(2.0*r1*r2*cospsi))*R_EARTH_KM
            if(ONE_D_MODEL=='PREM'.and.KFILE=='sph') then
                call model_prem_iso(r1,v1D,PHASE)
                call mantle_tomo(r1,ctheta1,phi1,dv,NK,NS)
                s1=1.0/v1D
                s1tomo=1.0/(v1D+v1D*dv)
                call model_prem_iso(r2,v1D,PHASE)
                call mantle_tomo(r2,ctheta2,phi2,dv,NK,NS)
                s2=1.0/v1D
                s2tomo=1.0/(v1D+v1D*dv)
                call model_prem_iso(r3,v1D,PHASE)
                call mantle_tomo(r3,ctheta3,phi3,dv,NK,NS)
                s3=1.0/v1D
                s3tomo=1.0/(v1D+v1D*dv)
            else if (ONE_D_MODEL=='STW105'.and.KFILE=='hvd') then
                if(r1 >= 0.d0 .and. r1 <= RICB/R_EARTH_KM) then
                    region = 3
                else if (r1 >= RICB/R_EARTH_KM .and. r1 < RCMB/R_EARTH_KM ) then
                    region = 2
                else
                    region = 1
                endif
                call model_1dref(r1,vv_1D_1,vh_1D_1,region,PHASE)
                xcolat = sngl(ctheta1*RADIANS_TO_DEGREES)
                xlon = sngl(phi1*RADIANS_TO_DEGREES)
                xrad = sngl(r1*R_EARTH_KM)
                call model_s362ani_subshsv(xcolat,xlon,xrad,dvh,dvv,PHASE)
                vv_3D_1=vv_1D_1*(1.0d0+dble(dvv))
                vh_3D_1=vh_1D_1*(1.0d0+dble(dvh))

                if(r2 >= 0.d0 .and. r2 <= RICB/R_EARTH_KM) then
                    region = 3
                else if (r2 >= RICB/R_EARTH_KM .and. r2 < RCMB/R_EARTH_KM ) then
                    region = 2
                else
                    region = 1
                endif

                call model_1dref(r2,vv_1D_2,vh_1D_2,region,PHASE)
                xcolat = sngl(ctheta2*RADIANS_TO_DEGREES)
                xlon = sngl(phi2*RADIANS_TO_DEGREES)
                xrad = sngl(r2*R_EARTH_KM)
                call model_s362ani_subshsv(xcolat,xlon,xrad,dvh,dvv,PHASE)
                vv_3D_2=vv_1D_2*(1.0d0+dble(dvv))
                vh_3D_2=vh_1D_2*(1.0d0+dble(dvh))

                if(r3 >= 0.d0 .and. r3 <= RICB/R_EARTH_KM) then
                    region = 3
                else if (r3 >= RICB/R_EARTH_KM .and. r3 < RCMB/R_EARTH_KM ) then
                    region = 2
                else
                    region = 1
                endif
                call model_1dref(r3,vv_1D_3,vh_1D_3,region,PHASE)
                xcolat = sngl(ctheta3*RADIANS_TO_DEGREES)
                xlon = sngl(phi3*RADIANS_TO_DEGREES)
                xrad = sngl(r3*R_EARTH_KM)
                call model_s362ani_subshsv(xcolat,xlon,xrad,dvh,dvv,PHASE)
                vv_3D_3=vv_1D_3*(1.0d0+dble(dvv))
                vh_3D_3=vh_1D_3*(1.0d0+dble(dvh))


                if(PHASE=='P') then
                    s1=1.0/sqrt((vv_1D_1**2+4.0*vh_1D_1**2)/5.0)
                    s1tomo=1.0/sqrt((vv_3D_1**2+4.0*vh_3D_1**2)/5.0)

                    s2=1.0/sqrt((vv_1D_2**2+4.0*vh_1D_2**2)/5.0)
                    s2tomo=1.0/sqrt((vv_3D_2**2+4.0*vh_3D_2**2)/5.0)

                    s3=1.0/sqrt((vv_1D_3**2+4.0*vh_1D_3**2)/5.0)
                    s3tomo=1.0/sqrt((vv_3D_3**2+4.0*vh_3D_3**2)/5.0)
                else if(PHASE=='S') then
                    s1=1.0/sqrt((2.0*(vv_1D_1**2)+vh_1D_1**2)/3.0)
                    s1tomo=1.0/sqrt((2.0*(vv_3D_1**2)+vh_3D_1**2)/3.0)

                    s2=1.0/sqrt((2.0*(vv_1D_2**2)+vh_1D_2**2)/3.0)
                    s2tomo=1.0/sqrt((2.0*(vv_3D_2**2)+vh_3D_2**2)/3.0)

                    s3=1.0/sqrt((2.0*(vv_1D_3**2)+vh_1D_3**2)/3.0)
                    s3tomo=1.0/sqrt((2.0*(vv_3D_3**2)+vh_3D_3**2)/3.0)
                else
                        stop 'PHASE problem'
                endif

            else
                stop'There is no model'
            endif

            ttheo=ttheo+((delta/6.0)*(s1+(4.0*s3)+s2)) !--> Simpson's method
            ttomo=ttomo+((delta/6.0)*(s1tomo+(4.0*s3tomo)+s2tomo)) !--> Simpson's method

        else if (INTEG=='trapezoidal') then
            r1=raypath(2,i)/R_EARTH_KM
            r2=raypath(2,i+1)/R_EARTH_KM

            theta1=raypath(3,i)*DEGREES_TO_RADIANS
            theta2=raypath(3,i+1)*DEGREES_TO_RADIANS

            ctheta1=(90.0-raypath(3,i))*DEGREES_TO_RADIANS !--> Colatitude
            ctheta2=(90.0-raypath(3,i+1))*DEGREES_TO_RADIANS !--> Colatitude

            phi1=raypath(4,i)*DEGREES_TO_RADIANS
            phi2=raypath(4,i+1)*DEGREES_TO_RADIANS

            cospsi=dsin(theta1)*dsin(theta2)+(dcos(theta1)*dcos(theta2)*dcos(phi2-phi1))
            delta=sqrt(r1**2+r2**2-(2.0*r1*r2*cospsi))*R_EARTH_KM
            if(ONE_D_MODEL=='PREM'.and.KFILE=='sph') then
                call model_prem_iso(r1,v1D,PHASE)
                call mantle_tomo(r1,ctheta1,phi1,dv,NK,NS)
                s1=1.0/v1D
                s1tomo=1.0/(v1D+v1D*dv)
                call model_prem_iso(r2,v1D,PHASE)
                call mantle_tomo(r2,ctheta2,phi2,dv,NK,NS)
                s2=1.0/v1D
                s2tomo=1.0/(v1D+v1D*dv)
            else if (ONE_D_MODEL=='STW105'.and.KFILE=='hvd') then
                call model_1dref_broadcast()
                call model_s362ani_broadcast(THREE_D_MODEL)
                if(r1 >= 0.d0 .and. r1 <= RICB/R_EARTH_KM) then
                    region = 3
                else if (r1 >= RICB/R_EARTH_KM .and. r1 < RCMB/R_EARTH_KM ) then
                    region = 2
                else
                    region = 1
                endif
                call model_1dref(r1,vv_1D_1,vh_1D_1,region,PHASE)
                xcolat = sngl(ctheta1*RADIANS_TO_DEGREES)
                xlon = sngl(phi1*RADIANS_TO_DEGREES)
                xrad = sngl(r1*R_EARTH_KM)
                call model_s362ani_subshsv(xcolat,xlon,xrad,dvh,dvv,PHASE)
                vv_3D_1=vv_1D_1*(1.0d0+dble(dvv))
                vh_3D_1=vh_1D_1*(1.0d0+dble(dvh))

                if(r2 >= 0.d0 .and. r2 <= RICB/R_EARTH_KM) then
                    region = 3
                else if (r2 >= RICB/R_EARTH_KM .and. r2 < RCMB/R_EARTH_KM ) then
                    region = 2
                else
                    region = 1
                endif

                call model_1dref(r2,vv_1D_2,vh_1D_2,region,PHASE)
                xcolat = sngl(ctheta1*RADIANS_TO_DEGREES)
                xlon = sngl(phi1*RADIANS_TO_DEGREES)
                xrad = sngl(r1*R_EARTH_KM)
                call model_s362ani_subshsv(xcolat,xlon,xrad,dvh,dvv,PHASE)
                vv_3D_2=vv_1D_2*(1.0d0+dble(dvv))
                vh_3D_2=vh_1D_2*(1.0d0+dble(dvh))

                if(PHASE=='P') then
                    s1=1.0/sqrt((vv_1D_1**2+4.0*vh_1D_1**2)/5.0)
                    s1tomo=1.0/sqrt((vv_3D_1**2+4.0*vh_3D_1**2)/5.0)

                    s2=1.0/sqrt((vv_1D_2**2+4.0*vh_1D_2**2)/5.0)
                    s2tomo=1.0/sqrt((vv_3D_2**2+4.0*vh_3D_2**2)/5.0)


                else if(PHASE == 'S') then
                    s1=1.0/sqrt((2.0*(vv_1D_1**2)+vh_1D_1**2)/3.0)
                    s1tomo=1.0/sqrt((2.0*(vv_3D_1**2)+vh_3D_1**2)/3.0)
                    s2=1.0/sqrt((2.0*(vv_1D_2**2)+vh_1D_2**2)/3.0)
                    s2tomo=1.0/sqrt((2.0*(vv_3D_2**2)+vh_3D_2**2)/3.0)


                else
                    stop 'PHASE problem'
                endif

            else
                stop 'There is no model'
            endif
            ttheo=ttheo+((delta/2.0)*(s1+s2)) !--> Trapezoidal method
            ttomo=ttomo+((delta/2.0)*(s1tomo+s2tomo)) !--> Trapezoidal method
        endif
        dtt=ttomo-ttheo !delay traveltime
enddo
!frmt="(i6.1, F10.3, F10.3, F10.3, F10.3, F10.3, F10.3, F10.3, F10.3)"
!write(*,frmt) int(info(1)), info(2), info(3), info(4), info(5), info(6), ttheo, ttomo, dtt
frmt="(F10.3, F10.3, F10.3)"
write(*,frmt) ttheo, ttomo, dtt

end program umich_cj_tcc

!
!--------------------------------------------------------------------------------------------------
!

subroutine model_broadcast(NK,NS,THREE_D_MODEL)

! standard routine to setup model

    use model_par

    implicit none

    integer ier
    integer :: NK,NS
    character(len=100) :: THREE_D_MODEL

    allocate(V_dv_a(0:NK,0:NS,0:NS), &
V_dv_b(0:NK,0:NS,0:NS), &
V_spknt(NK+1), &
V_qq0(NK+1,NK+1), &
V_qq(3,NK+1,NK+1), &
stat=ier)

    call read_model(NK,NS,THREE_D_MODEL)


end subroutine model_broadcast

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model(NK,NS,THREE_D_MODEL)

    use model_par

    implicit none
    integer :: NK,NS
    integer :: IIN = 1
    ! local parameters
    integer :: k,l,m,ier
    double precision :: dummy(4)
    character(len=100) :: THREE_D_MODEL

    open(unit=IIN,file='/home/samhaug/RayTimes_through3D/models/'&
        //trim(THREE_D_MODEL)//'.sph',status='old',action='read',iostat=ier)
    if (ier /= 0) then
        write(*,*) 'Error opening "', trim(THREE_D_MODEL), '": ', ier
        stop 'Unknown model'
    endif
    read(IIN,*) (dummy(k),k=1,4)
    do k = 0,NK
        do l = 0,NS
            read(IIN,*) V_dv_a(k,l,0),(V_dv_a(k,l,m),V_dv_b(k,l,m),m = 1,l)
        enddo
    enddo
    close(IIN)

    ! set up the splines used as radial basis functions by Ritsema
    call splhsetup(NK)

end subroutine read_model

!
!--------------------------------------------------------------------------------------------------
!
subroutine mantle_tomo(radius,theta,phi,dv,NK,NS)

    use model_par

    implicit none

    double precision :: radius,theta,phi,dv
    integer :: NK,NS
    ! local parameters
    double precision, parameter :: R_EARTH_ = R_EARTH ! 6371000.d0
    double precision, parameter :: RMOHO_ = R_EARTH - 24400.0 ! 6346600.d0
    double precision, parameter :: RCMB_ = 3480000.d0
    double precision, parameter :: ZERO_ = 0.d0

    integer :: l,m,k
    double precision :: r_moho,r_cmb,xr
    double precision :: dv_alm,dv_blm
    double precision :: rsple,radial_basis(0:NK)
    double precision :: sint,cost,x(2*NS+1),dx(2*NS+1)

    dv = ZERO_

    r_moho = RMOHO_ / R_EARTH_
    r_cmb = RCMB_ / R_EARTH_

    if (radius>=r_moho .or. radius <= r_cmb) return
    xr=-1.0d0+2.0d0*(radius-r_cmb)/(r_moho-r_cmb)
    if (xr > 1.0) print *,'xr > 1.0'
    if (xr < -1.0) print *,'xr < -1.0'
    do k = 0,NK
        radial_basis(k)=rsple(1,NK+1,V_spknt(1),V_qq0(1,NK+1-k),V_qq(1,1,NK+1-k),xr)
    enddo
    do l = 0,NS
        sint=dsin(theta)
        cost=dcos(theta)
        call lgndr(l,cost,sint,x,dx)
        dv_alm=0.0d0
        do k = 0,NK
            dv_alm=dv_alm+radial_basis(k)*V_dv_a(k,l,0)
        enddo
        dv=dv+dv_alm*x(1)
        do m = 1,l
            dv_alm=0.0d0
            dv_blm=0.0d0
            do k = 0,NK
                dv_alm=dv_alm+radial_basis(k)*V_dv_a(k,l,m)
                dv_blm=dv_blm+radial_basis(k)*V_dv_b(k,l,m)
            enddo
            dv=dv+(dv_alm*dcos(dble(m)*phi)+dv_blm*dsin(dble(m)*phi))*x(m+1)
        enddo
    enddo

end subroutine mantle_tomo

!----------------------------------

subroutine splhsetup(NK)

    use model_par

    implicit none


    integer :: NK
    ! local parameters
    integer :: i,j
    double precision :: qqwk(3,NK+1)

    V_spknt(1) = -1.00000d0
    V_spknt(2) = -0.78631d0
    V_spknt(3) = -0.59207d0
    V_spknt(4) = -0.41550d0
    V_spknt(5) = -0.25499d0
    V_spknt(6) = -0.10909d0
    V_spknt(7) = 0.02353d0
    V_spknt(8) = 0.14409d0
    V_spknt(9) = 0.25367d0
    V_spknt(10) = 0.35329d0
    V_spknt(11) = 0.44384d0
    V_spknt(12) = 0.52615d0
    V_spknt(13) = 0.60097d0
    V_spknt(14) = 0.66899d0
    V_spknt(15) = 0.73081d0
    V_spknt(16) = 0.78701d0
    V_spknt(17) = 0.83810d0
    V_spknt(18) = 0.88454d0
    V_spknt(19) = 0.92675d0
    V_spknt(20) = 0.96512d0
    V_spknt(21) = 1.00000d0

    do i = 1,NK+1
        do j = 1,NK+1
            if (i == j) then
                V_qq0(j,i)=1.0d0
            else
                V_qq0(j,i)=0.0d0
            endif
        enddo
    enddo
    do i = 1,NK+1
        call rspln(1,NK+1,V_spknt(1),V_qq0(1,i),V_qq(1,1,i),qqwk(1,1))
    enddo

end subroutine splhsetup

!----------------------------------

double precision function rsple(I1,I2,X,Y,Q,S)

    implicit none

! rsple returns the value of the function y(x) evaluated at point S
! using the cubic spline coefficients computed by rspln and saved in Q.
! If S is outside the interval (x(i1),x(i2)) rsple extrapolates
! using the first or last interpolation polynomial. The arrays must
! be dimensioned at least - x(i2), y(i2), and q(3,i2).

    integer i1,i2
    double precision  X(*),Y(*),Q(3,*),s

    integer i,ii
    double precision h

    i = 1
    II=I2-1

    !   GUARANTEE I WITHIN BOUNDS.
    I=MAX0(I,I1)
    I=MIN0(I,II)

    !   SEE IF X IS INCREASING OR DECREASING.
    if (X(I2)-X(I1) <  0) goto 1
    if (X(I2)-X(I1) >= 0) goto 2

    !   X IS DECREASING.  CHANGE I AS NECESSARY.
1   if (S-X(I) <= 0) goto 3
    if (S-X(I) >  0) goto 4

4   I=I-1

    if (I-I1 <  0) goto 11
    if (I-I1 == 0) goto 6
    if (I-I1 >  0) goto 1

3   if (S-X(I+1) <  0) goto 5
    if (S-X(I+1) >= 0) goto 6

5   I=I+1

    if (I-II <  0) goto 3
    if (I-II == 0) goto 6
    if (I-II >  0) goto 7

!   X IS INCREASING.  CHANGE I AS NECESSARY.
2   if (S-X(I+1) <= 0) goto 8
    if (S-X(I+1) >  0) goto 9

9   I=I+1

    if (I-II <  0) goto 2
    if (I-II == 0) goto 6
    if (I-II >  0) goto 7

8   if (S-X(I) <  0) goto 10
    if (S-X(I) >= 0) goto 6

10  I=I-1
    if (I-I1 <  0) goto 11
    if (I-I1 == 0) goto 6
    if (I-I1 >  0) goto 8

7   I=II
    goto 6
11  I=I1

!   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
6   H=S-X(I)
    RSPLE=Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))

end function rsple

!----------------------------------

subroutine rspln(I1,I2,X,Y,Q,F)

implicit none

! Subroutine rspln computes cubic spline interpolation coefficients
! for y(x) between grid points i1 and i2 saving them in q.The
! interpolation is continuous with continuous first and second
! derivatives. It agrees exactly with y at grid points and with the
! three point first derivatives at both end points (i1 and i2).
! X must be monotonic but if two successive values of x are equal
! a discontinuity is assumed and separate interpolation is done on
! each strictly monotonic segment. The arrays must be dimensioned at
! least - x(i2), y(i2), q(3,i2), and f(3,i2).
! F is working storage for rspln.

    integer i1,i2
    double precision X(*),Y(*),Q(3,*),F(3,*)

    integer i,j,k,j1,j2
    double precision y0,a0,b0,b1,h,h2,ha,h2a,h3a,h2b
    double precision YY(3),small
    equivalence (YY(1),Y0)
    data SMALL/1.0d-08/,YY/0.0d0,0.0d0,0.0d0/

    J1=I1+1
    Y0=0.0d0

    !   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL
    if (I2-I1  < 0) return
    if (I2-I1 == 0) goto 17
    if (I2-I1  > 0) goto 8

8   A0=X(J1-1)
!   SEARCH FOR DISCONTINUITIES.
    do  I=J1,I2
        B0=A0
        A0=X(I)
        if (DABS((A0-B0)/DMAX1(A0,B0)) < SMALL) goto 4
    enddo
17  J1=J1-1
    J2=I2-2
    goto 5
4    J1=J1-1
    J2=I-3
!   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
5   if (J2+1-J1 <  0) goto 9
    if (J2+1-J1 == 0) goto 10
    if (J2+1-J1 >  0) goto 11

!   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
10  J2=J2+2
    Y0=(Y(J2)-Y(J1))/(X(J2)-X(J1))
    do J=1,3
        Q(J,J1)=YY(J)
        Q(J,J2)=YY(J)
    enddo
    goto 12

!   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
11   A0 = 0.
    H=X(J1+1)-X(J1)
    H2=X(J1+2)-X(J1)
    Y0=H*H2*(H2-H)
    H=H*H
    H2=H2*H2
!   CALCULATE DERIVATIVE AT NEAR END.
    B0=(Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
    B1=B0

!   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
    do I=J1,J2
        H=X(I+1)-X(I)
        Y0=Y(I+1)-Y(I)
        H2=H*H
        HA=H-A0
        H2A=H-2.0d0*A0
        H3A=2.0d0*H-3.0d0*A0
        H2B=H2*B0
        Q(1,I)=H2/HA
        Q(2,I)=-HA/(H2A*H2)
        Q(3,I)=-H*H2A/H3A
        F(1,I)=(Y0-H*B0)/(H*HA)
        F(2,I)=(H2B-Y0*(2.0d0*H-A0))/(H*H2*H2A)
        F(3,I)=-(H2B-3.0d0*Y0*HA)/(H*H3A)
        A0=Q(3,I)
        B0=F(3,I)
    enddo

!   TAKE CARE OF LAST TWO ROWS.
    I=J2+1
    H=X(I+1)-X(I)
    Y0=Y(I+1)-Y(I)
    H2=H*H
    HA=H-A0
    H2A=H*HA
    H2B=H2*B0-Y0*(2.0d0*H-A0)
    Q(1,I)=H2/HA
    F(1,I)=(Y0-H*B0)/H2A
    HA=X(J2)-X(I+1)
    Y0=-H*HA*(HA+H)
    HA=HA*HA

!   CALCULATE DERIVATIVE AT FAR END.
    Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
    Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.0d0*A0))
    Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)

!   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
    do J=J1,J2
        K=I-1
        Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
        Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
        Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
        I=K
    enddo
    Q(1,I)=B1
!   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
9   J2=J2+2
    do J=1,3
        Q(J,J2)=YY(J)
    enddo

!   SEE IF THIS DISCONTINUITY IS THE LAST.
12  if (J2-I2 < 0) then
        goto 6
    else
        return
    endif

!   NO.  GO BACK FOR MORE.
6    J1=J2+2
    if (J1-I2 <= 0) goto 8
    if (J1-I2 >  0) goto 7

!   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
7    do J=1,3
        Q(J,I2)=YY(J)
    enddo

end subroutine rspln

!----------------------------------

subroutine lgndr(l,c,s,x,dx)

! computes Legendre function x(l,m,theta)
! theta=colatitude,c=cos(theta),s=sin(theta),l=angular order,
! sin(theta) restricted so that sin(theta) > 1.e-7
! x(1) contains m = 0, x(2) contains m = 1, x(k+1) contains m=k
! m=azimuthal (longitudinal) order 0 <= m <= l
! dx=dx/dtheta
!
! subroutine originally came from Physics Dept. Princeton through
! Peter Davis, modified by Jeffrey Park

    implicit none

! argument variables
    integer l
    double precision x(2*l+1),dx(2*l+1)
    double precision c,s

! local variables
    integer i,lp1,lpsafe,lsave
    integer m,maxsin,mmm,mp1

    double precision sqroot2over2,c1,c2,cot
    double precision ct,d,f1,f2
    double precision f3,fac,g1,g2
    double precision g3,rfpi,sqroot3,sos
    double precision ss,stom,t,tol
    double precision v,y

    tol = 1.d-05
    rfpi = 0.282094791773880d0
    sqroot3 = 1.73205080756890d0
    sqroot2over2 = 0.707106781186550d0

    if (s >= 1.0d0-tol) s=1.0d0-tol
    lsave=l
    if (l<0) l=-1-l
    if (l>0) goto 1
    x(1)=rfpi
    dx(1)=0.0d0
    l=lsave
    return
1 if (l /= 1) goto 2
    c1=sqroot3*rfpi
    c2=sqroot2over2*c1
    x(1)=c1*c
    x(2)=-c2*s
    dx(1)=-c1*s
    dx(2)=-c2*c
    l=lsave
    return
2 sos=s
    if (s<tol) s=tol
    cot=c/s
    ct=2.0d0*c
    ss=s*s
    lp1=l+1
    g3=0.0d0
    g2=1.0d0
    f3=0.0d0

! evaluate m=l value, sans (sin(theta))**l
    do i = 1,l
        g2=g2*(1.0d0-1.0d0/(2.0d0*i))
    enddo
    g2=rfpi*dsqrt((2*l+1)*g2)
    f2=l*cot*g2
    x(lp1)=g2
    dx(lp1)=f2
    v=1.0d0
    y=2.0d0*l
    d=dsqrt(v*y)
    t=0.0d0
    mp1=l
    m=l-1

! these recursions are similar to ordinary m-recursions, but since we
! have taken the s**m factor out of the xlm's, the recursion has the powers
! of sin(theta) instead
3 g1=-(ct*mp1*g2+ss*t*g3)/d
    f1=(mp1*(2.0d0*s*g2-ct*f2)-t*ss*(f3+cot*g3))/d-cot*g1
    x(mp1)=g1
    dx(mp1)=f1
    if (m == 0) goto 4
    mp1=m
    m=m-1
    v=v+1.0d0
    y=y-1.0d0
    t=d
    d=dsqrt(v*y)
    g3=g2
    g2=g1
    f3=f2
    f2=f1
    goto 3
! explicit conversion to integer added
4 maxsin=int(-72.0d0/log10(s))

! maxsin is the max exponent of sin(theta) without underflow
    lpsafe=min0(lp1,maxsin)
    stom=1.0d0
    fac=sign(1.0d0,dble((l/2)*2-l) + 0.50d0)

! multiply xlm by sin**m
    do m = 1,lpsafe
        x(m)=fac*x(m)*stom
        dx(m)=fac*dx(m)*stom
        stom=stom*s
    enddo

! set any remaining xlm to zero
    if (maxsin <= l) then
        mmm=maxsin+1
        do m=mmm,lp1
            x(m)=0.0d0
            dx(m)=0.0d0
        enddo
    endif

    s=sos
    l=lsave

end subroutine lgndr
