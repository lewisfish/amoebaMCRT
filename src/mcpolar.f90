module monte

    implicit none

    contains

    subroutine mcpolar(concstmp, numproc, id, pflag, iq, iw, src, comm)

        use mpi_f08


        !shared data
        use constants
        use photon_vars
        use iarray
        use opt_prop

        !subroutines
        use subs
        use reader_mod
        use gridset_mod
        use sourceph_mod
        use inttau2
        use ch_opt
        use stokes_mod
        use fluorophores
        use writer_mod
        use utils, only : str

        implicit none

        real,           intent(IN)  :: concstmp(:)
        real,           intent(OUT) :: src(:)
        integer,        intent(IN)  :: numproc, id, iq, iw
        type(mpi_comm), intent(IN)  :: comm
        logical,        intent(IN)  :: pflag

        integer :: q, w, e, i
        integer :: nphotons, jseed, j, xcell, ycell, zcell, fluro_img(1000,3), fluroglobal(1000,3)
        logical :: tflag
        real    :: ran, delta, start, finish, wave_in, ran2
        character(len=256) :: paramsFile
        type(fluro), allocatable :: f_array(:)

        !mpi variables
        real    :: finish2, tmp, mu_tot, mua_skin

        !allocate and set arrays to 0
        call alloc_array(numproc)

        call zarray
        src = 0.d0
        !**** Read in parameters from the file input.params
        open(10,file=trim(resdir)//'input.params',status='old')
        read(10,*) nphotons
        read(10,*) xmax
        read(10,*) ymax
        read(10,*) zmax
        read(10,*) n1
        read(10,*) n2
        read(10,*) paramsFile
        close(10)

        !set optical properties and make cdfs.
        wave_in = 365.d0 
        hgg = gethgg(wave_in)

        g2 = hgg**2
        wave    = wave_in
        wave_incd = wave !set incdent wavelength

        call init_fluros(f_array, trim(paramsFile), id, comm)

        call gridset(f_array, id)

        !***** Set small distance for use in optical depth integration routines 
        !***** for roundoff effects when crossing cell walls
        delta = 1.e-8*(2.*zmax/nzg)

        call set_default_opt(f_array) !set up default optical properties
        call opt_set_incident()!set optical properties to default

        call mpi_barrier(comm)
        call cpu_time(start)


        q = 2
        w = 2
        e = 2
        blood_pap = .06d0
        blood_ret = .045d0
        blood_hypo = .05d0
        jseed = -123456789+id!-97463834+id

        call cpu_time(start)
        call set_default_opt(f_array) !set up default optical properties
        fluro_img = 0
        fluroglobal = 0
        do j = 1, nphotons

            e = 1
            wave = wave_in
            f_array(:)%bool = .false.
            hgg = gethgg(wave)
            g2 = hgg*hgg
            call opt_set_incident()

            
            tflag=.FALSE.

            if(mod(j,1000000) == 0 .and. .not. pflag)then
                print *, j,' scattered photons completed on core: ',id
            end if

            if(j == 1000000 .and. id == 0 .and. pflag)then
                call cpu_time(finish2)
                print*,' '
                tmp = (finish2-start)/1000000.*real(nphotons)
                if(tmp >= 60.)then
                    tmp = tmp / 60.
                    if(tmp > 60)then
                        tmp = tmp / 60.
                        print*,str(tmp),' hrs'
                    else
                        print*,str(tmp),' mins'
                    end if
                else
                    print*, str(tmp),' s'
                end if
            end if

            !***** Release photon from point source *******************************
            call sourceph(xcell,ycell,zcell,jseed)

            !****** Find scattering location
            call tauint1(xcell,ycell,zcell,tflag,jseed,delta, id, e)

            !******** Photon scatters in grid until it exits (tflag=TRUE) 
            do while(tflag.eqv..FALSE.)

                mua_skin = rhokap(zcell) - albedo(zcell)
                mu_tot  = albedo(zcell) + mua_skin

                do i = 1, size(f_array)
                    f_array(i)%mua = ln10 * f_array(i)%getFluroWave(wave) * conc(zcell, i)
                    mu_tot = mu_tot + f_array(i)%mua
                end do

                ran = ran2(jseed)

                if(ran < f_array(1)%mua / mu_tot)then
                    !do nadh
                    call sample(f_array(1)%emission, f_array(1)%cdf, wave, jseed)
                     call opt_set(wave, f_array)
                    f_array(1)%bool = .true.
                    e = 1
                    f_array(2:)%bool = .false.
                    hgg = 0.d0
                    g2 = 0.d0
                    call stokes(jseed)
                    hgg = gethgg(wave)
                    g2 = hgg**2
                elseif(ran < (f_array(1)%mua + f_array(2)%mua)/ mu_tot)then
                    !do fad
                    call sample(f_array(2)%emission, f_array(2)%cdf,wave, jseed)
                    call opt_set(wave, f_array)
                    f_array(2)%bool = .true.
                    e = 2
                    f_array(1)%bool = .false.
                    f_array(3)%bool = .false.
                    hgg = 0.d0
                    g2 = 0.d0
                    call stokes(jseed)
                    hgg = gethgg(wave)
                    g2 = hgg**2 
                elseif(ran < (f_array(1)%mua + f_array(2)%mua + f_array(3)%mua)/ mu_tot)then
                    !do ribofake
                    call sample(f_array(3)%emission, f_array(3)%cdf,wave, jseed)
                    call opt_set(wave, f_array)
                    f_array(3)%bool = .true.
                    e = 3
                    f_array(2)%bool = .false.
                    f_array(1)%bool = .false.
                    hgg = 0.d0
                    g2 = 0.d0
                    call stokes(jseed)
                    hgg = gethgg(wave)
                    g2 = hgg**2   
                elseif(ran < (f_array(1)%mua + f_array(2)%mua + f_array(3)%mua + albedo(zcell))/ mu_tot)then
                    !do scatter
                    call stokes(jseed)
                else
                    !absorb
                    tflag=.true.
                    exit
                end if

                !************ Find next scattering location
                call tauint1(xcell,ycell,zcell,tflag,jseed,delta, id, e)

            end do

            if(tflag .and. wave /= wave_in .and. zcell == -1 .and. zp > 0.)then
                    if(f_array(1)%bool)fluro_img(nint(wave),1) = fluro_img(nint(wave),1) + 1
                    if(f_array(2)%bool)fluro_img(nint(wave),2) = fluro_img(nint(wave),2) + 1
                    if(f_array(3)%bool)fluro_img(nint(wave),3) = fluro_img(nint(wave),3) + 1
            end if

        end do      ! end loop over nph photons

        ! allocate(jmeanGLOBAL(nxg, nyg, nzg,2))


        call mpi_allreduce(fluro_img, fluroglobal, size(fluroglobal), mpi_integer, mpi_sum, comm)
        ! call mpi_allreduce(jmean, jmeanglobal, size(jmeanglobal), mpi_double_precision, mpi_sum, comm)

        src = fluroglobal(:,1) + fluroglobal(:,2) + fluroglobal(:,3)
        if(id==0)print*,sum(src)

        ! src = src / max(real(maxval(src)), 1.)
        if(id == 0)then

            ! open(newunit=j,file=trim(fileplace)//"all_out1.dat")
            ! do i = 1, size(fluro_img,1)
            !     write(j,*)i,fluro_img(i,:)
            ! end do
            ! close(j)
            call writer(src)
        end if

        call cpu_time(finish)

        call dealloc_array()
        call mpi_barrier(comm)
end subroutine mcpolar
end module monte




                ! ! print*,ran,albedo, mu_nadh/kappa
                ! if(ran < albedo(zcell))then!interacts with tissue
                !     call stokes(jseed)
                ! else
                !     if(ran < (ln10 * NADH(wave) * conc(zcell,2))/rhokap(zcell) + albedo(zcell))then
                !         call sample(nadh_fluro, nadh_cdf, wave, jseed)
                !         call opt_set(wave)
                !         nadh_bool=.true.;fad_bool=.false.;collagen_bool=.false.

                !     elseif(ran < ((ln10 * NADH(wave) * conc(zcell,2)) + &
                !                   (ln10 * fadflavin(wave) * conc(zcell,3)))/rhokap(zcell) + albedo(zcell))then
                !         call sample(riboflavin_fluro, riboflavin_cdf, wave, jseed)
                !         call opt_set(wave)
                !         riboflavin_bool=.true.;nadh_bool=.false.;tyrosine_bool=.false.
                !     elseif(ran < (ln10 * NADH(wave) * conc(zcell, 2)/rhokap(zcell)) + &
                !                  (ln10 * riboflavin(wave) * conc(zcell,3)) + &
                !                  (ln10 * tyrosine(wave) * conc(zcell,1)) + albedo(zcell))then
                !         call sample(tyrosine_fluro, tyrosine_cdf, wave, jseed)
                !         call opt_set(wave)
                !         riboflavin_bool=.false.;nadh_bool=.false.;tyrosine_bool=.true.
                !     else
                !         tflag=.true.
                !         exit
                !     end if
                !     hgg = 0.
                !     call stokes(jseed)
                !     hgg = .8
                ! end if