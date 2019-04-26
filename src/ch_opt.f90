module ch_opt

    implicit none

    contains

        subroutine opt_set(wave, f_array)

            use skin_layers
            use iarray,       only : albedo, rhokap
            use opt_prop,     only : strat_range, epi_range, pap_range, ret_range, hypo_range
            use fluorophores, only : fluro

            implicit none

            real,        intent(IN) :: wave
            type(fluro), intent(IN) :: f_array(:)

            real    :: alb_kap(2)

            !Strat corenum
            alb_kap = stratum(f_array, wave, strat_range(1))
            rhokap(strat_range(1):strat_range(2)) = alb_kap(2)
            albedo(strat_range(1):strat_range(2)) = alb_kap(1)

            !Living Epidermis
            alb_kap = Epidermis(f_array, wave, epi_range(1))
            rhokap(epi_range(1):epi_range(2)) = alb_kap(2)
            albedo(epi_range(1):epi_range(2)) = alb_kap(1)

            !PaPillary Dermis
            alb_kap = Pap_dermis(f_array, wave, pap_range(1))
            rhokap(pap_range(1):pap_range(2)) = alb_kap(2)
            albedo(pap_range(1):pap_range(2)) = alb_kap(1)

            !Reticular dermis
            alb_kap = Ret_Dermis(f_array, wave, ret_range(1))
            rhokap(ret_range(1):ret_range(2)) = alb_kap(2)
            albedo(ret_range(1):ret_range(2)) = alb_kap(1)  

            !Hypodermis
            alb_kap = Hypo_Dermis(f_array, wave, hypo_range(1))
            rhokap(hypo_range(1):hypo_range(2)) = alb_kap(2)
            albedo(hypo_range(1):hypo_range(2)) = alb_kap(1) 
        end subroutine opt_set


        subroutine set_default_opt(f_array)
        ! save the "default" optical properties to an array for faster program execution

            use fluorophores, only : fluro
            use opt_prop,     only : wave_incd, opt_array
            use skin_layers,  only : stratum, Epidermis, Pap_dermis, Ret_Dermis, Hypo_Dermis

            implicit none

            type(fluro), intent(IN) :: f_array(:)

            opt_array = 0.d0

            opt_array(1, :) = stratum(f_array, wave_incd, 500)
            opt_array(2, :) = Epidermis(f_array, wave_incd, 497)
            opt_array(3, :) = Pap_dermis(f_array, wave_incd, 489)
            opt_array(4, :) = Ret_Dermis(f_array, wave_incd, 471)
            opt_array(5, :) = Hypo_Dermis(f_array, wave_incd, 290)

        end subroutine set_default_opt


        subroutine opt_set_incident() 

            use skin_layers
            use iarray,    only : albedo, rhokap, zface
            use constants, only : nzg, zmax
            use opt_prop,  only : opt_array

            implicit none


            real    :: alb_kap(2), z
            integer :: i


            do i = 1, nzg
                z = zface(i) + zmax/nzg

                if(z <= 0.02)then
                    !Strat corenum
                    alb_kap = opt_array(1, :)
                    rhokap(i) = alb_kap(2)
                    albedo(i) = alb_kap(1)
                elseif(z <= 0.02 + 0.08)then
                    !Living Epidermis
                    alb_kap = opt_array(2, :)
                    rhokap(i) = alb_kap(2)
                    albedo(i) = alb_kap(1)
                elseif(z <= 0.02 + 0.08 + 0.18)then
                    !PaPillary Dermis
                    alb_kap = opt_array(3, :)
                    rhokap(i) = alb_kap(2)
                    albedo(i) = alb_kap(1)
                elseif(z <= 0.02 + 0.08 + 0.18 + 1.82)then
                    alb_kap = opt_array(4, :)
                    rhokap(i) = alb_kap(2)
                    albedo(i) = alb_kap(1)  
                else
                    !Hypodermis
                    alb_kap = opt_array(5, :)
                    rhokap(i) = alb_kap(2)
                    albedo(i) = alb_kap(1) 
                end if
            end do

        end subroutine opt_set_incident


        real function gethgg(wave)

            implicit none

            real, intent(IN) :: wave

            gethgg = 0.62d0 + (wave * 0.29d-3)

        end function gethgg

        subroutine writeOut(a,s,d, f_array)

            use skin_layers
            use constants,    only : fileplace
            use utils,        only : str
            use fluorophores, only : fluro

            implicit none

            integer, intent(IN) :: a, s, d
            type(fluro), intent(IN) :: f_array(:)

            integer :: i,q,w,e,r,t
            real    :: tmp(2), mua, mus

            print*,trim(fileplace)//'stratum-'//str(a)//"-"//str(s)//"-"//str(d)//'.dat'
            open(newunit=q,file=trim(fileplace)//'stratum-'//str(a)//"-"//str(s)//"-"//str(d)//'.dat',&
                 status="replace")
            open(newunit=w,file=trim(fileplace)//'epidermis-'//str(a)//"-"//str(s)//"-"//str(d)//'.dat',&
                 status="replace")
            open(newunit=e,file=trim(fileplace)//'paPillary-'//str(a)//"-"//str(s)//"-"//str(d)//'.dat',&
                 status="replace")
            open(newunit=r,file=trim(fileplace)//'reticular-'//str(a)//"-"//str(s)//"-"//str(d)//'.dat',&
                 status="replace")
            open(newunit=t,file=trim(fileplace)//'hypodermis-'//str(a)//"-"//str(s)//"-"//str(d)//'.dat',&
                 status="replace")

            !writes out wavelength, mua, mus

            do i = 300, 700
                tmp = stratum(f_array, real(i), 500)
                mua = tmp(2) - tmp(1)*tmp(2)
                mus =  tmp(1)*tmp(2)
                write(q,*)i,mua, mus

                tmp = epidermis(f_array, real(i), 497)
                mua = tmp(2) - tmp(1)*tmp(2)
                mus =  tmp(1)*tmp(2)
                write(w,*)i,mua, mus

                tmp = Pap_Dermis(f_array, real(i), 489)
                mua = tmp(2) - tmp(1)*tmp(2)
                mus =  tmp(1)*tmp(2)
                write(e,*)i,mua, mus

                tmp = Ret_Dermis(f_array, real(i), 471)
                mua = tmp(2) - tmp(1)*tmp(2)
                mus =  tmp(1)*tmp(2)
                write(r,*)i,mua, mus

                tmp = Hypo_Dermis(f_array, real(i), 290)
                mua = tmp(2) - tmp(1)*tmp(2)
                mus =  tmp(1)*tmp(2)
                write(t,*)i,mua, mus
            end do
            close(q)
            close(w)
            close(e)
            close(r)
            close(t)
        end subroutine writeOut
end module ch_opt