module generateFluence

    implicit none
        
    contains
    
        subroutine getFluenceDetectedPackets(delta, id)

            use constants,   only : fileplace, xmax, ymax, zmax
            use inttau2,     only : tauint1
            use iarray,      only : jmean
            use photon_vars, only : nxp, nyp, nzp, xp, yp, zp

            implicit none

            real,    intent(IN) :: delta
            integer, intent(IN) :: id

            integer            :: u, io, e, jseed, xcell, ycell, zcell, k
            character(len=256) :: line
            logical            :: tflag
            real               :: tau

            !read in photon position (xp, yp, zp), vector (nxp, nyp, nzp) and tau used
            open(newunit=u, file=trim(fileplace)//"photPos.dat", status="old")
            k = 1
            do
                read(u,"(a256)", iostat=io)line
                if(IS_IOSTAT_END(io))exit
                read(line, "(7(F10.7,1x),I1)")xp, yp, zp, nxp, nyp, nzp, tau,e

                !calculate initial cell
                xcell = int(250*(xp+xmax) / (2.*xmax)) + 1
                ycell = int(250*(yp+ymax) / (2.*ymax)) + 1
                zcell = int(500*(zp+zmax) / (2.*zmax)) + 1

                !set flags and seed
                tflag = .False.
                jseed = -231213132
                !rewind photon
                if(mod(k, 10000) == 0)print*,k,"9968789"
                if(e /= 0)call tauint1(xcell,ycell,zcell,tflag,jseed,delta, id, e, tau, .true.)
                k = k + 1
            end do
            close(u)

            !save rewound photons Fluence
            open(newunit=u,file=trim(fileplace)//"detc-flue-all.dat",access="stream",form="unformatted")
            write(u)jmean
            close(u)

        end subroutine getFluenceDetectedPackets

end module generateFluence