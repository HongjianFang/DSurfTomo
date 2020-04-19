subroutine voronoiproj(leniw,lenrw,iw,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                       nx,ny,nz,nd,ncells,hvratio,damp,iproj,dv)
      use lsmrModule, only:lsmr

      implicit none
      integer leniw,lenrw
      integer nx,ny,nz
      integer iw(leniw)
      real depz(nz)
      real rw(lenrw)
      integer ncells
      real dv(*),dres(*)
      real goxd,gozd,dvxd,dvzd
      real damp
      real hvratio,cmb
      integer ndim,nd
      integer iproj

      real,parameter:: radius = 6371.0,ftol = 0.1,pi = 3.141592654
      integer ii,ix,iy,iz
      real,dimension(:),allocatable:: grow,gcol,subrow,dis,dws,xunknown
      real,dimension(:),allocatable:: lat,lon,rad,theta,phi,rrad,xpts,ypts,zpts
      real,dimension(:),allocatable :: rw_p,rwgp,norm
      integer,dimension(:),allocatable:: iw_p,row,col,iwgp,colgp
      integer idx
      integer maxnar,nzid
      integer iseed(4)
      real xs,ys,zs

      real atol,btol
      real conlim
      integer istop
      integer itnlim
      real acond
      real anorm
      real arnorm
      real rnorm
      real xnorm
      integer localSize,nout,itn
      integer leniw_p,lenrw_p,leniwgp,lenrwgp

      allocate(lat(nx-2),lon(ny-2),rad(nz-1))
      ndim = (nx-2)*(ny-2)*(nz-1)

      do ii = 1,nx-2
      lat(ii) = (goxd-(ii-1)*dvxd)*pi/180
      enddo

      do ii = 1,ny-2
      lon(ii) = (gozd+(ii-1)*dvzd)*pi/180
      enddo

      !cmb = radius - depz(nz-1)*hvratio
      do ii = 1,nz-1
      rad(ii) = radius-depz(ii)*hvratio
      !rad(ii) = cmb+depz(ii)*hvratio
      enddo

      allocate(theta(ncells),phi(ncells),rrad(ncells),dws(ncells),norm(ncells))
      allocate(xpts(ncells),ypts(ncells),zpts(ncells),dis(ncells),xunknown(ncells))
      allocate(rw_p(ndim))
      allocate(iw_p(2*ndim+1),row(ndim),col(ndim))

      iseed(1:3) = (/38,62,346/)
      iseed(4) = 2*iproj+1
      call slarnv(1,iseed,ncells,theta)
      theta = (gozd+theta*(ny-3)*dvzd)*pi/180
      call slarnv(1,iseed,ncells,phi)
      phi = pi/2-(goxd-phi*(nx-3)*dvxd)*pi/180
      call slarnv(1,iseed,ncells,rrad)
      rrad = radius-rrad*depz(nz-1)*hvratio

      xpts = rrad*sin(phi)*cos(theta)
      ypts = rrad*sin(phi)*sin(theta)
      zpts = rrad*cos(phi)

      idx = 0
      do iz = 1,nz-1
      do iy = 1,ny-2
      do ix = 1,nx-2
      xs = rrad(iz)*sin(pi/2-lat(ix))*cos(lon(iy))
      ys = rrad(iz)*sin(pi/2-lat(ix))*sin(lon(iy))
      zs = rrad(iz)*cos(pi/2-lat(ix))
      dis =  (xpts-xs)**2+(ypts-ys)**2+(zpts-zs)**2
      idx = idx+1
      col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
      row(idx) = minloc(dis,1)
      enddo
      enddo
      enddo
      rw_p = 1.0
      leniw_p = 2*ndim+1
      lenrw_p = ndim
      iw_p(1) = ndim
      iw_p(2:ndim+1) = row
      iw_p(ndim+2:2*ndim+1) = col

      allocate(grow(ndim),gcol(nd),subrow(ncells))
      maxnar = int(0.6*nd*ncells)
      allocate(iwgp(maxnar*2+1),colgp(maxnar),rwgp(maxnar))

      nzid = 0
      do ii = 1,nd
      grow = 0
      gcol = 0
      gcol(ii) = 1.0 
      call aprod(2,nd,ndim,grow,gcol,leniw,lenrw,iw,rw)
      subrow = 0
      call aprod(1,ncells,ndim,grow,subrow,leniw_p,lenrw_p,iw_p,rw_p)
      do ix = 1,ncells
      if(abs(subrow(ix))>ftol) then
          nzid = nzid+1
          rwgp(nzid) = subrow(ix)
          iwgp(1+nzid) = ii
          colgp(nzid) = ix
      endif
      enddo
      enddo
      leniwgp = nzid*2+1
      lenrwgp = nzid
      iwgp(1) = lenrwgp
      iwgp(nzid+2:nzid*2+1) = colgp(1:nzid)

      dws = 0
      !do ii=1,nzid
      !dws(iwgp(1+ii+nzid)) = dws(iwgp(1+ii+nzid))+abs(rwgp(ii))
      !enddo
      norm = 0
      do ii=1,nzid
      norm(iwgp(1+ii+nzid)) = norm(iwgp(1+ii+nzid))+rwgp(ii)**2
      enddo
    
      do ii =1,ncells
      norm(ii) = sqrt(norm(ii)/nd+0.01)
      enddo
    
      do ii =1,nzid
      rwgp(ii) = rwgp(ii)/norm(iwgp(1+ii+nzid))
      enddo

      conlim = 50
      itnlim = 100
      atol = 1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
      btol = 1e-3/(dvxd*nx*111.19/3.0)!1e-3
      istop = 0
      anorm = 0.0
      acond = 0.0
      arnorm = 0.0
      xnorm = 0.0
      localSize = int(ncells/4)
      !damp = dampvel
      ! using lsmr to solve for the projection coefficients
      !print*, 'LSMR beginning ...'

      !nout = -1
      nout = 36
      open(nout,file='lsmrout_sub.txt')

      call LSMR(nd, ncells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
      close(nout)
      do ii = 1,ncells
        xunknown(ii) = xunknown(ii)/norm(ii)
      enddo
      dv(1:ndim) = 0
      call aprod(2,ncells,ndim,dv,xunknown,leniw_p,lenrw_p,iw_p,rw_p)
      deallocate(grow,gcol,subrow)
      deallocate(theta,phi,rrad,dws,norm)
      deallocate(xpts,ypts,zpts,dis,xunknown)
      deallocate(iw_p,rw_p,row,col)
      deallocate(lat,lon,rad)
      deallocate(iwgp,colgp,rwgp)

      end subroutine
