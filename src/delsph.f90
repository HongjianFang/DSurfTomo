subroutine delsph(flat1,flon1,flat2,flon2,del)
implicit none
real,parameter:: R=6371.0
REAL,parameter:: pi=3.1415926535898
real flat1,flat2
real flon1,flon2
real del

real dlat
real dlon
real lat1
real lat2
real a
real c


!dlat=(flat2-flat1)*pi/180
!dlon=(flon2-flon1)*pi/180
!lat1=flat1*pi/180
!lat2=flat2*pi/180
dlat=flat2-flat1
dlon=flon2-flon1
lat1=pi/2-flat1
lat2=pi/2-flat2
a=sin(dlat/2)*sin(dlat/2)+sin(dlon/2)*sin(dlon/2)*cos(lat1)*cos(lat2)
c=2*atan2(sqrt(a),sqrt(1-a))
del=R*c
end subroutine
