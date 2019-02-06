module skin_layers

    !blood values + mel value form louise/craigs paper
    !water, caro, bili from iglesias

    ! equations mix of my own and mm02 + iglesias + kb04

    implicit none

    contains
        function stratum(f_array, wave, zcell)

            use constants,    only : ln10
            use absorbers,    only : base, water, carotene, fat, bilirubin, Eumel, Pheomel, Deoxy_Hb, Oxy_Hb
            use opt_prop,     only : hgg
            use iarray,       only : conc
            use fluorophores, only : fluro!: nadh, fad!, tyrosine

            implicit none

            type(fluro), intent(IN) :: f_array(:)
            real,    intent(IN) :: wave
            integer, intent(IN) :: zcell

            integer :: i
            real :: stratum(2), a, bdash, fray
            real :: mua, mus_r, mus_m, mus, W, Baseline

            if(wave >= 100.0d0)then

                !set mua
                W        = 0.05d0      !water level
                Baseline = 0.d0 
                mua = ((0.1d0 - 0.3d-4*wave) + 0.125d0*wave/10.d0*base(wave))*(1.d0 - W) + W*water(wave)
                
                !set fluorophores
                do i = 1, size(f_array)
                    mua = mua + ln10 * (f_array(i)%getFluroWave(wave) * conc(zcell, i))
                end do

                !set mus
                ! formula from Jac13
                a     = 66.7d0
                fray  = .29d0
                bdash     = .689d0
                mus_r = fray*(wave/500.d0)**(-4.)
                mus_m = (1.d0 - fray)*(wave/500.d0)**(-bdash)
                mus   = a * (mus_r + mus_m)
                mus   = mus/(1.d0 - hgg)

                Stratum(2) = mua + mus
                Stratum(1) = mus / stratum(2)

            else
                print*,wave,'Wavelength out of range!'
                print*,'Setting Stratum mua = 0!!!'
                Stratum = 0.d0
            end if
        end function stratum


        function Epidermis(f_array, wave, zcell)

            use constants,    only : ln10
            use absorbers,    only : base, water, carotene, fat, bilirubin, Eumel, Pheomel, Deoxy_Hb, Oxy_Hb
            use opt_prop,     only : hgg
            use iarray,       only : conc
            use fluorophores, only : fluro

            implicit none

            type(fluro), intent(IN) :: f_array(:)
            real,    intent(IN) :: wave
            integer, intent(IN) :: zcell

            integer :: i
            real :: Epidermis(2), a, bdash, fray
            real :: mua, mus, mus_r, mus_m
            real :: nu_m , W, S, B, baseline
            real :: C_bili, C_caro

            if(wave .ge. 100.d0)then
                !set mua
                nu_m     = .01d0
                W        = 0.20d0     !water level
                S        = 0.75d0     !Blood oxygenation
                B        = 0.0d0      !blood fraction
                C_caro   = 2.1d-4
                C_bili   = 0.0d0
                Baseline = 0.3d0 

                mua = (nu_m * (Eumel(wave) + Pheomel(wave)) + (1.d0 - nu_m) * (ln10*C_caro*carotene(wave) + (1.d0-C_caro)*&
                       base(wave))) * (1.d0 - w) + w*water(wave)

                !set fluorophores
                do i = 1, size(f_array)
                    mua = mua + ln10 * (f_array(i)%getFluroWave(wave) * conc(zcell, i))
                end do


                !set mus
                ! formula from Jac13
                a     = 66.7d0
                fray  = 0.29d0
                bdash = 0.689d0
                mus_r = fray*(wave/500.d0)**(-4.)
                mus_m = (1.d0 - fray)*(wave/500.d0)**(-bdash)
                mus   = a * (mus_r + mus_m)
                mus   = mus/(1.d0 - hgg)

                Epidermis(2) = mua + mus
                Epidermis(1) = mus / Epidermis(2)

            else
                print*,wave,'Wavelength out of range!'
                print*,'Setting Epidermis mua = 0!!!'
                Epidermis = 0.d0
            end if
        end function Epidermis


        function Pap_dermis(f_array, wave, zcell)
            
            use constants,    only : ln10
            use absorbers,    only : base, water, carotene, fat, bilirubin, Eumel, Pheomel, Deoxy_Hb, Oxy_Hb
            use opt_prop,     only : hgg, blood_pap
            use iarray,       only : conc
            use fluorophores, only : fluro


            implicit none

            type(fluro), intent(IN) :: f_array(:)
            real,        intent(IN)    :: wave
            integer,     intent(IN) :: zcell

            real    :: Pap_dermis(2), C_bili, C_caro, nu_m, mus_m, mus_r
            real    :: mua, mus, a, fray, bdash,  W, S, B, baseline
            integer :: i

            if(wave .ge. 100.d0)then
                !set mua
                nu_m     = 0.0d0
                W        = 0.5d0      !water level
                S        = 0.75d0     !Blood oxygenation
                B        = blood_pap  !blood fraction
                C_caro   = 7.0d-5
                C_bili   = 0.05d0
                Baseline = 0.3d0


                mua = (B * (Oxy_Hb(wave) * S + Deoxy_Hb(wave) * (1.d0 - S) + ln10 * C_caro * carotene(wave) +&
                           ln10 * C_bili * bilirubin(wave)) + base(wave)*(1.d0 - B)) * (1.d0-W)+w*water(wave)

                
                !set fluorophores
                do i = 1, size(f_array)
                    mua = mua + ln10 * (f_array(i)%getFluroWave(wave) * conc(zcell, i))
                end do

                !set mus
                ! formula from Jac13
                a     = 43.6d0
                fray  = 0.41d0
                bdash = 0.35d0
                mus_r = fray*(wave/500.d0)**(-4.)
                mus_m = (1.d0 - fray)*(wave/500.d0)**(-bdash)
                mus   = a * (mus_r + mus_m)
                mus   = mus/(1.d0 - hgg)

                Pap_dermis(2) = mua + mus
                Pap_dermis(1) = mus / Pap_dermis(2)
            else
                print*,wave,'Wavelength out of range!'
                print*,'Setting Pap_dermis mua = 0!!!'
                Pap_dermis = 0.d0
            end if
        end function Pap_dermis


        function Ret_dermis(f_array, wave, zcell)

            use constants,    only : ln10
            use absorbers,    only : base, water, carotene, fat, bilirubin, Eumel, Pheomel, Deoxy_Hb, Oxy_Hb
            use opt_prop,     only : hgg, blood_ret
            use iarray,       only : conc
            use fluorophores, only : fluro

            implicit none

            type(fluro), intent(IN) :: f_array(:)
            integer,     intent(IN) :: zcell
            real,        intent(IN) :: wave

            integer :: i
            real    :: mua, mus, a, fray, bdash,  W, S, B, baseline
            real    :: C_bili, C_caro, nu_m, Ret_dermis(2), mus_r, mus_m

            if(wave .ge. 100.d0)then
                !set mua
                nu_m     = 0.0d0
                W        = 0.7d0      !water level
                S        = 0.75d0     !Blood oxygenation
                B        = blood_ret  !blood fraction
                C_caro   = 7.0d-5
                C_bili   = 0.05d0
                Baseline = 0.3d0

                mua = (B * (Oxy_Hb(wave) * S + Deoxy_Hb(wave) * (1.d0 - S) + ln10 * C_caro * carotene(wave) +&
                           ln10 * C_bili * bilirubin(wave)) + base(wave)*(1.d0 - B)) * (1.d0-W)+w*water(wave)

                !set fluorophores
                do i = 1, size(f_array)
                    mua = mua + ln10 * (f_array(i)%getFluroWave(wave) * conc(zcell, i))
                end do

                !set mus
                ! formula from Jac13
                a     = 43.6d0
                fray  = 0.41d0
                bdash = 0.35d0
                mus_r = fray*(wave/500.d0)**(-4.)
                mus_m = (1.d0 - fray)*(wave/500.d0)**(-bdash)
                mus   = a * (mus_r + mus_m)
                mus   = mus/(1.d0 - hgg)


                Ret_dermis(2) = mua + mus
                Ret_dermis(1) = mus / Ret_dermis(2)
            else
                print*,'Wavelength out of range!'
                print*,'Setting Ret_dermis mua = 0!!!'
                Ret_dermis = 0.d0
            end if
        end function Ret_dermis


        function Hypo_dermis(f_array, wave, zcell)   

            use constants,    only : ln10
            use absorbers,    only : base, water, carotene, fat, bilirubin, Eumel, Pheomel, Deoxy_Hb, Oxy_Hb
            use opt_prop,     only : hgg, blood_hypo
            use fluorophores, only : fluro
            use iarray,       only : conc

            implicit none

            type(fluro), intent(IN) :: f_array(:)
            real,        intent(IN) :: wave
            integer,     intent(IN) :: zcell

            integer :: i
            real    :: hypo_dermis(2), C_bili, C_caro, nu_m
            real    :: mua, mus, W, S, B, baseline

            if(wave .ge. 100.d0)then
                !set mua
                nu_m     = 0.0d0
                W        = 0.7d0      !water level
                S        = 0.75d0     !Blood oxygenation
                B        = blood_hypo !blood fraction
                C_caro   = 0.0d0
                C_bili   = 0.0d0
                Baseline = 0.3d0

                mua = (B * (Oxy_Hb(wave) * S + Deoxy_Hb(wave) * (1.d0 - S) + ln10 * C_caro * carotene(wave) +&
                           ln10 * C_bili * bilirubin(wave)) + base(wave)*(1.d0 - B)) * (1.d0-W) + w*water(wave) 

                !set mus
                mus = 1050.6d0 * wave**(-0.68) !in cm-1
                mus = mus / (1.0d0 - hgg)

                !set fluorophores
                do i = 1, size(f_array)
                    mua = mua + ln10 * (f_array(i)%getFluroWave(wave) * conc(zcell, i))
                end do

                Hypo_dermis(2) = mua + mus
                Hypo_dermis(1) = mus / hypo_dermis(2)
            else
                print*,'Wavelength out of range!'
                print*,'Setting Hypo_dermis mua = 0!!!'
                Hypo_dermis = 0.0d0
            end if
        end function hypo_dermis
end module skin_layers