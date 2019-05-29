module sourceph_mod

   implicit none

   contains
      subroutine sourceph(xcell, ycell, zcell, iseed)

         use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, TWOPI
         use photon_vars

         implicit none

         integer, intent(OUT)   :: xcell, ycell, zcell
         integer, intent(INOUT) :: iseed

         real :: ran2, theta, r, spotsizeradius

         spotsizeradius = 0.04d0 ! 0.4 mm

         zp = zmax - (1.d-5 * (2.d0*zmax/nzg))

         !http://mathworld.wolfram.com/DiskPointPicking.html
         !sample circle uniformly
         !sample radius between [0,r^2]
         r = ran2(iseed) * (spotsizeradius)**2
         theta = ran2(iseed) * twopi
         xp = sqrt(r) * cos(theta)
         yp = sqrt(r) * sin(theta)


         phi = TWOPI * ran2(iseed)
         cosp = cos(phi)
         sinp = sin(phi)          
         cost = -1.d0 
         sint = sqrt(1. - cost**2) 

         nxp = sint * cosp  
         nyp = sint * sinp
         nzp = cost

         !*************** Linear Grid *************************
         xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
         ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
         zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
         !*****************************************************
      end subroutine sourceph
end module sourceph_mod
