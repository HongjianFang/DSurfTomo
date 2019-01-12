! CODE FOR SURFACE WAVE TOMOGRAPHY USING DISPERSION MEASUREMENTS
! VERSION:
!      1.0 
! AUTHOR:
!      HONGJIAN FANG. fanghj@mail.ustc.edu.cn
! PURPOSE:
!      DIRECTLY INVERT SURFACE WAVE DISPERSION MEASUREMENTS FOR 3-D
! STUCTURE WITHOUT THE INTERMEDIATE STEP OF CONSTUCTION THE PHASE
! OR GROUP VELOCITY MAPS.
! REFERENCE:
! Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D.
! (2015). Direct inversion of surface wave dispersion for
! three-dimensional shallow crustal structure based on ray tracing:
! methodology and application. Geophysical Journal International,
! 201(3), 1251-1263.
! HISTORY:
!       2015/01/31 START TO REORGONIZE THE MESSY CODE   
!

program SurfTomo
  use lsmrModule, only:lsmr
  use  lsmrblasInterface, only : dnrm2
  implicit none

  ! VARIABLE DEFINE

  character inputfile*80
  character logfile*100
  character outmodel*100
  character outsyn*100
  logical ex
  character dummy*40
  character datafile*80

  integer nx,ny,nz
  real goxd,gozd
  real dvxd,dvzd
  integer nsrc,nrc
  real weight,weight0
  real damp
  real minthk
  integer kmax,kmaxRc,kmaxRg,kmaxLc,kmaxLg
  real*8,dimension(:),allocatable:: tRc,tRg,tLc,tLg
  real,dimension(:),allocatable:: depz
  integer itn
  integer nout
  integer localSize
  real mean,std_devs,balances,balanceb
  integer msurf
  real,dimension(:),allocatable:: obst,dsyn,cbst,wt,dtres,dist,datweight
  real,dimension(:),allocatable:: pvall,depRp,pvRp
  real sta1_lat,sta1_lon,sta2_lat,sta2_lon
  real dist1
  integer dall
  integer istep
  real,parameter :: pi=3.1415926535898
  integer checkstat
  integer ii,jj,kk
  real, dimension (:,:), allocatable :: scxf,sczf
  real, dimension (:,:,:), allocatable :: rcxf,rczf
  integer,dimension(:,:),allocatable::wavetype,igrt,nrc1
  integer,dimension(:),allocatable::nsrc1
  integer,dimension(:,:),allocatable::periods
  real,dimension(:),allocatable::rw
  integer,dimension(:),allocatable::iw,col
  real,dimension(:),allocatable::dv,norm
  real,dimension(:,:,:),allocatable::vsf
  real,dimension(:,:,:),allocatable::vsftrue
  character strf
  integer veltp,wavetp
  real velvalue
  integer knum,knumo,err
  integer istep1,istep2
  integer period
  integer knumi,srcnum,count1
  integer HorizonType,VerticalType
  character line*200
  integer iter,maxiter
  integer maxnar
  real acond
  real anorm
  real arnorm
  real rnorm
  real xnorm
  character str1
  real atol,btol
  real conlim
  integer istop
  integer itnlim
  integer lenrw,leniw
  integer nar,nar_tmp,nars
  integer count3,nvz,nvx
  integer m,maxvp,n
  integer i,j,k
  real Minvel,MaxVel
  real spfra
  real noiselevel
  integer ifsyn
  integer writepath
  real averdws
  real maxnorm
  real threshold0

! FOR MODEL VARIATION
  !------------------------------------------------
  integer idx
  integer counte
  real stdvs
!  integer numrand
!  real,allocatable,dimension(:,:)::modstat
!  real,allocatable,dimension(:)::modsig
  real gaussian
  external gaussian
  integer modest
  counte=0

  ! OPEN FILES FIRST TO OUTPUT THE PROCESS
  open(34,file='IterVel.out')
  nout=36
  open(nout,file='lsmr.txt')

  ! OUTPUT PROGRAM INFOMATION            
  write(*,*)
  write(*,*),'                             DSurfTomo (v1.3)'
  write(*,*),'PLEASE contact Hongjain Fang &
    (fanghj@mail.ustc.edu.cn) if you find any bug'
  write(*,*)

  ! READ INPUT FILE
  if (iargc() < 1) then
    write(*,*) 'input file [DSurfTomo.in(default)]:'
    read(*,'(a)') inputfile
    if (len_trim(inputfile) <=1 ) then
      inputfile = 'DSurfTomo.in'
    else
      inputfile = inputfile(1:len_trim(inputfile))
    endif
  else
    call getarg(1,inputfile)
  endif
  inquire(file = inputfile, exist = ex)
  if (.not. ex) stop 'unable to open the inputfile'

  open(10,file=inputfile,status='old')
  read(10,'(a30)')dummy
  read(10,'(a30)')dummy
  read(10,'(a30)')dummy
  read(10,*)datafile
  read(10,*) nx,ny,nz
  read(10,*) goxd,gozd
  read(10,*) dvxd,dvzd
  read(10,*) nsrc
  read(10,*) weight0,damp
  read(10,*) minthk
  read(10,*) Minvel,Maxvel
  read(10,*) maxiter
  read(10,*) spfra
  read(10,*) kmaxRc
  write(*,*)  'model origin:latitude,longitue'
  write(*,'(2f10.4)') goxd,gozd
  write(*,*) 'grid spacing:latitude,longitue'
  write(*,'(2f10.4)') dvxd,dvzd
  write(*,*) 'model dimension:nx,ny,nz'
  write(*,'(3i5)') nx,ny,nz
  write(logfile,'(a,a)')trim(inputfile),'.log' 
  open(66,file=logfile)
  write(66,*)
  write(66,*),'                         S U R F  T O M O'
  write(66,*),'PLEASE contact Hongjain Fang &
    (fanghj@mail.ustc.edu.cn) if you find any bug'
  write(66,*)
  write(66,*) 'model origin:latitude,longitue'
  write(66,'(2f10.4)') goxd,gozd
  write(66,*) 'grid spacing:latitude,longitue'
  write(66,'(2f10.4)') dvxd,dvzd
  write(66,*) 'model dimension:nx,ny,nz'
  write(66,'(3i5)') nx,ny,nz
  if(kmaxRc.gt.0)then
    allocate(tRc(kmaxRc),&
      stat=checkstat)
    if (checkstat > 0) stop 'error allocating RP'
    read(10,*)(tRc(i),i=1,kmaxRc)
    write(*,*)'Rayleigh wave phase velocity used,periods:(s)'
    write(*,'(50f7.1)')(tRc(i),i=1,kmaxRc)
    write(66,*)'Rayleigh wave phase velocity used,periods:(s)'
    write(66,'(50f7.1)')(tRc(i),i=1,kmaxRc)
  endif
  read(10,*)kmaxRg
  if(kmaxRg.gt.0)then
    allocate(tRg(kmaxRg), stat=checkstat)
    if (checkstat > 0) stop 'error allocating RP'
    read(10,*)(tRg(i),i=1,kmaxRg)
    write(*,*)'Rayleigh wave group velocity used,periods:(s)'
    write(*,'(50f7.1)')(tRg(i),i=1,kmaxRg)
    write(66,*)'Rayleigh wave group velocity used,periods:(s)'
    write(66,'(50f7.1)')(tRg(i),i=1,kmaxRg)
  endif
  read(10,*)kmaxLc
  if(kmaxLc.gt.0)then
    allocate(tLc(kmaxLc), stat=checkstat)
    if (checkstat > 0) stop 'error allocating RP'
    read(10,*)(tLc(i),i=1,kmaxLc)
    write(*,*)'Love wave phase velocity used,periods:(s)'
    write(*,'(50f7.1)')(tLc(i),i=1,kmaxLc)
    write(66,*)'Love wave phase velocity used,periods:(s)'
    write(66,'(50f7.1)')(tLc(i),i=1,kmaxLc)
  endif
  read(10,*)kmaxLg
  if(kmaxLg.gt.0)then
    allocate(tLg(kmaxLg), stat=checkstat)
    if (checkstat > 0) stop 'error allocating RP'
    read(10,*)(tLg(i),i=1,kmaxLg)
    write(*,*)'Love wave group velocity used,periods:(s)'
    write(*,'(50f7.1)')(tLg(i),i=1,kmaxLg)
    write(66,*)'Love wave group velocity used,periods:(s)'
    write(66,'(50f7.1)')(tLg(i),i=1,kmaxLg)
  endif
  read(10,*)ifsyn
  read(10,*)noiselevel
  read(10,*) threshold0
!  read(10,*) modest
!  read(10,*) numrand
  close(10)
  nrc=nsrc
  kmax=kmaxRc+kmaxRg+kmaxLc+kmaxLg

  ! READ MEASUREMENTS            
  open(unit=87,file=datafile,status='old')
  allocate(scxf(nsrc,kmax),sczf(nsrc,kmax),&
    rcxf(nrc,nsrc,kmax),rczf(nrc,nsrc,kmax),stat=checkstat)
  if(checkstat > 0)then
    write(6,*)'error with allocate'
  endif
  allocate(periods(nsrc,kmax),wavetype(nsrc,kmax),&
    nrc1(nsrc,kmax),nsrc1(kmax),&
    igrt(nsrc,kmax),stat=checkstat)
  if(checkstat > 0)then
    write(6,*)'error with allocate'
  endif
  allocate(obst(nrc*nsrc*kmax),dist(nrc*nsrc*kmax),&
    stat=checkstat)
  allocate(pvall(nrc*nsrc*kmax),depRp(nrc*nsrc*kmax),&
    pvRp(nrc*nsrc*kmax),stat=checkstat)
  IF(checkstat > 0)THEN
    write(6,*)'error with allocate'
  ENDIF
  istep=0
  istep2=0
  dall=0
  knumo=12345
  knum=0
  istep1=0
  do 
    read(87,'(a)',iostat=err) line
    if(err.eq.0) then
      if(line(1:1).eq.'#') then
        read(line,*) str1,sta1_lat,sta1_lon,period,wavetp,veltp
        if(wavetp.eq.2.and.veltp.eq.0) knum=period
        if(wavetp.eq.2.and.veltp.eq.1) knum=kmaxRc+period
        if(wavetp.eq.1.and.veltp.eq.0) knum=kmaxRg+kmaxRc+period
        if(wavetp.eq.1.and.veltp.eq.1) knum=kmaxLc+kmaxRg+&
          kmaxRc+period
        if(knum.ne.knumo) then
          istep=0
          istep2=istep2+1
        endif
        istep=istep+1
        istep1=0
        sta1_lat=(90.0-sta1_lat)*pi/180.0
        sta1_lon=sta1_lon*pi/180.0
        scxf(istep,knum)=sta1_lat
        sczf(istep,knum)=sta1_lon
        periods(istep,knum)=period
        wavetype(istep,knum)=wavetp
        igrt(istep,knum)=veltp
        nsrc1(knum)=istep
        knumo=knum
      else
        read(line,*) sta2_lat,sta2_lon,velvalue
        istep1=istep1+1
        dall=dall+1
        sta2_lat=(90.0-sta2_lat)*pi/180.0
        sta2_lon=sta2_lon*pi/180.0
        rcxf(istep1,istep,knum)=sta2_lat
        rczf(istep1,istep,knum)=sta2_lon
        call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist1)
        dist(dall)=dist1
        obst(dall)=dist1/velvalue
        pvall(dall)=velvalue
        nrc1(istep,knum)=istep1
      endif
    else
      exit
    endif
  enddo
  close(87)
  allocate(depz(nz), stat=checkstat)
  maxnar = spfra*dall*nx*ny*nz!sparsity fraction
  maxvp = (nx-2)*(ny-2)*(nz-1)
  allocate(dv(maxvp), stat=checkstat)
  allocate(norm(maxvp), stat=checkstat)
  allocate(vsf(nx,ny,nz), stat=checkstat)
  allocate(vsftrue(nx,ny,nz), stat=checkstat)
  ! FOR MODEL VARIATION
  !------------------------------------------------
!  allocate(modstat(numrand,maxvp))
!  allocate(modsig(maxvp))

  allocate(rw(maxnar), stat=checkstat)
  if(checkstat > 0)then
    write(6,*)'error with allocate: real rw'
  endif
  allocate(iw(2*maxnar+1), stat=checkstat)
  if(checkstat > 0)then
    write(6,*)'error with allocate: integer iw'
  endif
  allocate(col(maxnar), stat=checkstat)
  if(checkstat > 0)then
    write(6,*)'error with allocate:  integer iw'
  endif
  allocate(cbst(dall+maxvp),dsyn(dall),datweight(dall),wt(dall+maxvp),dtres(dall+maxvp),&
    stat=checkstat)

  ! MEASUREMENTS STATISTICS AND READ INITIAL MODEL           
  write(*,'(a,i7)') 'Number of all measurements',dall

  open(10,file='MOD',status='old')
  read(10,*) (depz(i),i=1,nz)
  do k = 1,nz
    do j = 1,ny
      read(10,*)(vsf(i,j,k),i=1,nx)
    enddo
  enddo
  close(10)
  write(*,*) 'grid points in depth direction:(km)'
  write(*,'(50f7.1)') depz



  ! CHECKERBOARD TEST
  if (ifsyn == 1) then
    write(*,*) 'Checkerboard Resolution Test Begin'
    vsftrue = vsf

    open(11,file='MOD.true',status='old')
    do k = 1,nz
      do j = 1,ny
        read(11,*) (vsftrue(i,j,k),i=1,nx)
      enddo
    enddo
    close(11)

    call synthetic(nx,ny,nz,maxvp,vsftrue,obst,&
      goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
      tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
      scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
      nsrc,nrc,noiselevel)
  endif



  ! ITERATE UNTILL CONVERGE
  writepath = 0
  do iter = 1,maxiter
    iw = 0
    rw = 0.0
    col = 0

    ! COMPUTE SENSITIVITY MATRIX
    if (iter == maxiter) then
      writepath = 1
      open(40,file='raypath.out')
    endif
    write(*,*) 'computing sensitivity matrix...'
    call CalSurfG(nx,ny,nz,maxvp,vsf,iw,rw,col,dsyn,&
      goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
      tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
      scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
      nsrc,nrc,nar,writepath)

    do i = 1,dall
      cbst(i) = obst(i) - dsyn(i)
    enddo

    ! write out rw, iw and cbst for testing
    !open(44,file='iw.dat')
    !write(44,*) nar
    !do i = 2,nar+1
    !write(44,*) iw(i)
    !enddo
    !do i = 1,nar
    !write(44,*) col(i)
    !enddo
    !close(44)

    !open(44,file='rw.dat')
    !do i = 1,nar
    !write(44,*) rw(i)
    !enddo
    !close(44)

    !open(44,file='residal.dat')
    !do i = 1,dall
    !write(44,*) cbst(i)
    !enddo
    !close(44)

    !threshold = threshold0+(maxiter/2-iter)/3*0.5
    do i = 1,dall
!      datweight(i) = 1.0
!      if(abs(cbst(i)) > threshold) then
!        datweight(i) = exp(-(abs(cbst(i))-threshold))
!      endif
      datweight(i) = 0.01+1.0/(1+0.05*exp(cbst(i)**2*threshold0))
      cbst(i) = cbst(i)*datweight(i)
    enddo

    do i = 1,nar
      rw(i) = rw(i)*datweight(iw(1+i))
    enddo

    norm=0
    do i=1,nar
      norm(col(i))=norm(col(i))+abs(rw(i))
    enddo
    averdws=0
    maxnorm=0
    do i=1,maxvp
      averdws =  averdws+norm(i)
      if(norm(i)>maxnorm) maxnorm=norm(i)
    enddo
    averdws=averdws/maxvp
    write(66,*)'Maximum and Average DWS values:',maxnorm,averdws
    !write(66,*)'Threshold is:',threshold

    ! WRITE OUT RESIDUAL FOR THE FIRST AND LAST ITERATION
    if(iter.eq.1) then
      open(88,file='residualFirst.dat')
      do i=1,dall
        write(88,*) dist(i),dsyn(i),obst(i), &
          dsyn(i)*datweight(i),obst(i)*datweight(i),datweight(i)
      enddo
      close(88)
    endif
    if(iter.eq.maxiter) then
      open(88,file='residualLast.dat')
      do i=1,dall
        write(88,*) dist(i),dsyn(i),obst(i), &
          dsyn(i)*datweight(i),obst(i)*datweight(i),datweight(i)
      enddo
      close(88)
    endif


    ! ADDING REGULARIZATION TERM
    !weight=dnrm2(dall,cbst,1)**2/dall*weight0
    weight=weight0
    nar_tmp=nar
    nars=0

    count3=0
    nvz=ny-2
    nvx=nx-2
    do k=1,nz-1
      do j=1,nvz
        do i=1,nvx
          if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
            count3=count3+1
            col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
            rw(nar+1)=2.0*weight
            iw(1+nar+1)=dall+count3
            cbst(dall+count3)=0
            nar=nar+1
          else
            count3=count3+1
            col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
            rw(nar+1)=6.0*weight
            iw(1+nar+1)=dall+count3
            rw(nar+2)=-1.0*weight
            iw(1+nar+2)=dall+count3
            col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1
            rw(nar+3)=-1.0*weight
            iw(1+nar+3)=dall+count3
            col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1
            rw(nar+4)=-1.0*weight
            iw(1+nar+4)=dall+count3
            col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i
            rw(nar+5)=-1.0*weight
            iw(1+nar+5)=dall+count3
            col(nar+5)=(k-1)*nvz*nvx+j*nvx+i
            rw(nar+6)=-1.0*weight
            iw(1+nar+6)=dall+count3
            col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i
            rw(nar+7)=-1.0*weight
            iw(1+nar+7)=dall+count3
            col(nar+7)=k*nvz*nvx+(j-1)*nvx+i
            cbst(dall+count3)=0
            nar=nar+7
          endif
        enddo
      enddo
    enddo
    m = dall + count3
    n = maxvp

    iw(1)=nar
    do i=1,nar
      iw(1+nar+i)=col(i)
    enddo
    if (nar > maxnar) stop 'increase sparsity fraction(spfra)'

    ! CALLING IRLS TO SOLVE THE PROBLEM

    leniw = 2*nar+1
    lenrw = nar
    dv = 0
    atol = 1e-4
    btol = 1e-4
    conlim = 100
    itnlim = 400
    istop = 0
    anorm = 0.0
    acond = 0.0
    arnorm = 0.0
    xnorm = 0.0
    localSize = 10

    call LSMR(m, n, leniw, lenrw,iw,rw,cbst, damp,&
      atol, btol, conlim, itnlim, localSize, nout,&
      dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
    !if(istop==3) print*,'istop = 3, large condition number'

    mean = sum(cbst(1:dall))/dall
    std_devs = sqrt(sum(cbst(1:dall)**2)/dall - mean**2)
    write(*,'(i2,a)'),iter,'th iteration...'
    write(*,'(a,f7.3)'),'weight is:',weight
    write(*,'(a,f8.1,a,f8.2,a,f8.3)'),'mean,std_devs and rms of &
      residual after weighting: ',mean*1000,'ms ',1000*std_devs,'ms ',&
      dnrm2(dall,cbst,1)/sqrt(real(dall))

    do i =1,dall
       cbst(i)=cbst(i)/datweight(i)
    enddo

    mean = sum(cbst(1:dall))/dall
    std_devs = sqrt(sum(cbst(1:dall)**2)/dall - mean**2)
    !write(*,'(i2,a)'),iter,'th iteration...'
    !write(*,'(a,f7.3)'),'weight is:',weight
    write(*,'(a,f8.1,a,f8.2,a,f8.3)'),'residual before weighting: ',mean*1000,'ms ',1000*std_devs,'ms ',&
      dnrm2(dall,cbst,1)/sqrt(real(dall))
    write(66,'(i2,a)'),iter,'th iteration...'
    write(66,'(a,f7.3)'),'weight is:',weight
    write(66,'(a,f8.1,a,f8.2,a,f8.3)'),'mean,std_devs and rms of &
      residual: ',mean*1000,'ms ',1000*std_devs,'ms ',&
      dnrm2(dall,cbst,1)/sqrt(real(dall))

    write(*,'(a,2f7.4)'),'min and max velocity variation ',&
      minval(dv),maxval(dv)
    write(66,'(a,2f7.4)'),'min and max velocity variation ',&
      minval(dv),maxval(dv)

    do k=1,nz-1
      do j=1,ny-2
        do i=1,nx-2
          if(dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i).ge.0.500) then
            dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)=0.500
          endif
          if(dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i).le.-0.500) then
            dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)=-0.500
          endif
          vsf(i+1,j+1,k)=vsf(i+1,j+1,k)+dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
          if(vsf(i+1,j+1,k).lt.Minvel) vsf(i+1,j+1,k)=Minvel
          if(vsf(i+1,j+1,k).gt.Maxvel) vsf(i+1,j+1,k)=Maxvel
        enddo
      enddo
    enddo
    write(34,*)'OUTPUT S VELOCITY AT ITERATION',iter
    do k=1,nz
      do j=1,ny
        write(34,'(100f7.3)') (vsf(i,j,k),i=1,nx)
      enddo
    enddo
    write(34,*)',OUTPUT DWS AT ITERATION',iter
    do k=1,nz-1
      do j=2,ny-1
        write(34,'(100f8.1)') (norm((k-1)*(ny-2)*(nx-2)+(j-2)*(nx-2)+i-1),i=2,nx-1)
      enddo
    enddo

    write(outmodel,'(a,a,i3.3)') trim(inputfile),'Measure.dat.iter',iter
    open(64,file=outmodel)
    do k=1,nz-1
      do j=1,ny-2
        do i=1,nx-2
          write(64,'(5f9.3)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsf(i+1,j+1,k)
        enddo
      enddo
    enddo
    close(64)


  enddo !end iteration

  ! OUTPUT THE VELOCITY MODEL

  write(*,*),'Program finishes successfully'
  write(66,*),'Program finishes successfully'

  if(ifsyn == 1) then
    open(65,file='Vs_model.real')
    write(outsyn,'(a,a)') trim(inputfile),'Syn.dat'
    open(63,file=outsyn)
    do k=1,nz-1
      do j=1,ny-2
        do i=1,nx-2
          write(65,'(5f9.3)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsftrue(i+1,j+1,k)
          write(63,'(5f9.3)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsf(i+1,j+1,k)
        enddo
      enddo
    enddo
    close(65)
    close(63)
    write(*,*),'Output True velocity model &
      to Vs_model.real'
    write(*,*),'Output inverted shear velocity model &
      to ',outsyn
    write(66,*),'Output True velocity model &
      to Vs_model.real'
    write(66,*),'Output inverted shear velocity model &
      to ',outsyn
  else
    write(outmodel,'(a,a)') trim(inputfile),'Measure.dat'
    open(64,file=outmodel)
    do k=1,nz-1
      do j=1,ny-2
        do i=1,nx-2
          write(64,'(5f9.3)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),vsf(i+1,j+1,k)
        enddo
      enddo
    enddo
    close(64)
    write(*,*),'Output inverted shear velocity model &
      to ',outmodel
    write(66,*),'Output inverted shear velocity model &
      to ',outmodel
  endif

  close(34)
  close(40)
  close(nout) !close lsmr.txt
  close(66) !close surf_tomo.log

!! USE RANDOM MODEL TO OBTAIN THE MODEL VARIATION
!  !modest = 1
!  if (modest ==1) then
!
!    write(*,*) 'model variation estimation begin...'
!    do iter = 1,numrand
!      call init_random_seed()
!      vsftrue=vsf
!      DO K=1,NZ-1
!        DO J=2,NY-1
!          DO I=2,NX-1
!            idx = (k-1)*(ny-2)*(nx-2)+(j-2)*(nx-2)+i-1
!            dv(idx) = 0.1/EXP(2*NORM(idx)/maxnorm)*gaussian()
!            VSFTRUE(I,J,K) = VSF(I,J,K)+dv(idx)
!          ENDDO
!        ENDDO
!      ENDDO
!      write(*,*),'maximum and minimum velocity variation',maxval(dv),minval(dv)
!
!      call synthetic(nx,ny,nz,maxvp,vsftrue,dsyn,&
!        goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
!        tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
!        scxf,sczf,rcxf,rczf,nrc1,nsrc1,kmax,&
!        nsrc,nrc,0.0)
!
!      do i = 1,dall
!        cbst(i) = obst(i) - dsyn(i)
!      enddo
!
!      write(*,*), dnrm2(dall,cbst,1)/sqrt(real(dall)), 1.05*std_devs
!      if (dnrm2(dall,cbst,1)/sqrt(real(dall)) < 1.05*std_devs) then
!        counte = counte + 1
!        modstat(counte,:) = dv
!      endif
!
!    enddo ! iteration for random models
!
!    write(*,*),'number of of models satisfy requirements',counte
!    modsig = 1.0
!    if (counte>0) then
!      do i=1,maxvp
!        !statis
!        !mean = sum(cbst(1:dall))/dall
!        !std_devs = sqrt(sum(cbst(1:dall)**2)/dall - mean**2)
!        mean = sum(modstat(1:counte,i))/counte
!        stdvs = sqrt(sum(modstat(1:counte,i)**2)/counte-mean**2) 
!        modsig(i) = stdvs
!      enddo
!    endif
!
!    write(*,*),'write model variation to "model_variation.dat"'
!    open (64,file='model_variation.dat')
!    do k=1,nz-1
!      do j=1,ny-2
!        do i=1,nx-2
!          idx = (k-1)*(ny-2)*(nx-2)+(j-1)*(nx-2)+i
!          write(64,'(5f8.4)') gozd+(j-1)*dvzd,goxd-(i-1)*dvxd,depz(k),modsig(idx)
!        enddo
!      enddo
!    enddo
!    close(64)
!    write(*,*) 'finishing model variation estimation'
!  endif



  deallocate(obst)
  deallocate(dsyn)
  deallocate(dist)
  deallocate(depz)
  deallocate(scxf,sczf)
  deallocate(rcxf,rczf)
  deallocate(wavetype,igrt,nrc1)
  deallocate(nsrc1,periods)
  deallocate(rw)
  deallocate(iw,col)
  deallocate(cbst,wt,dtres,datweight)
  deallocate(dv)
  deallocate(norm)
  deallocate(vsf)
  deallocate(vsftrue)
  if(kmaxRc.gt.0) then
    deallocate(tRc)
  endif
  if(kmaxRg.gt.0) then
    deallocate(tRg)
  endif
  if(kmaxLc.gt.0) then
    deallocate(tLc)
  endif
  if(kmaxLg.gt.0) then
    deallocate(tLg)
  endif

end program            



!-----------------------------------------------------------------------
!     Generate seed for random number generator of fortran
!     Note: only need to be called once inside one program
!-----------------------------------------------------------------------
subroutine init_random_seed()
  integer :: i,n,clock
  integer,dimension(:),allocatable :: seed
  call random_seed(size=n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed=clock+37*(/(i-1,i=1,n)/)
  call random_seed(PUT=seed) 
  deallocate(seed)
end subroutine
