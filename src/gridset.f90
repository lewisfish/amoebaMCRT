module gridset_mod

    implicit none

    contains

        subroutine gridset(f_array, id)
        !set up concentrion of fluorophores into layers
        !calculate layer ranges in voxels for speed up in other parts of the code

            use fluorophores, only : fluro
            use constants,    only : nxg, nyg, nzg, xmax, ymax, zmax
            use iarray,       only : xface, yface, zface, conc, refrac
            use opt_prop,     only : Strat_range, epi_range, pap_range, ret_range, hypo_range

            implicit none

            type(fluro), intent(IN) :: f_array(:)

            real    :: z
            integer :: i, strat_last, epi_last, pap_last, ret_last, hypo_last,id, j
            integer :: strat_first, epi_first, pap_first, ret_first, hypo_first
            logical :: strat_bool, epi_bool, pap_bool, ret_bool, hypo_bool

            strat_first = 0
            epi_first = 0
            pap_first = 0
            ret_first = 0
            hypo_first = 0
            strat_last = 0
            epi_last = 0
            pap_last = 0
            ret_last = 0
            hypo_last = 0

            strat_bool = .false.
            epi_bool = .false.
            pap_bool = .false.
            ret_bool = .false.
            hypo_bool = .false.

            !**********  Linear Cartesian grid. Set up grid faces ****************
            do i = 1 , nxg + 1
                xface(i) = (i - 1.d0) * 2.d0 *xmax/nxg
            end do
            do i = 1 , nyg + 1
                yface(i) = (i - 1.d0) * 2.d0 *ymax/nyg
            end do
            do i = 1 , nzg + 1
                zface(i) = (i - 1.d0) * 2.d0 *zmax/nzg
            end do

            refrac(nzg + 1) = 1.0
            conc = 0.d0

            do i = 1, nzg
                z = zface(i)+2.*zmax/nzg
                !distance here are converted from mm to cm.src is louise age paper
                if(z >= 2.*zmax - 0.02d-1)then
                    !Strat corenum
                    if(.not. strat_bool)then
                        strat_bool = .true.
                        strat_first = i
                    end if
                    strat_last = i
                    refrac(i) = 1.5d0 
                elseif(z >= 2.*zmax - (0.02d-1 + 0.08d-1))then
                    !Living Epidermis
                    if(.not. epi_bool)then
                        epi_bool = .true.
                        epi_first = i
                    end if
                    epi_last = i
                    refrac(i) = 1.34d0
                elseif(z >= 2.*zmax - (0.02d-1 + 0.08d-1 + 0.18d-1))then
                    !PaPillary Dermis
                    if(.not. pap_bool)then
                        pap_bool = .true.
                        pap_first = i
                    end if
                    pap_last = i
                    refrac(i) = 1.4d0
                elseif(z >= 2.*zmax - (0.02d-1 + 0.08d-1 + 0.18d-1 + 1.82d-1))then
                    !Reticular Dermis   
                    if(.not. ret_bool)then
                        ret_bool = .true.
                        ret_first = i
                    end if 
                    ret_last = i
                    refrac(i) = 1.395d0
                else
                    !Hypodermis
                    if(.not. hypo_bool)then
                        hypo_bool = .true.
                        hypo_first = i
                    end if
                    hypo_last = i
                    refrac(i) = 1.41d0
                end if
        end do

        Strat_range = [strat_first, strat_last]
        epi_range = [epi_first, epi_last]
        pap_range = [pap_first, pap_last]
        ret_range = [ret_first, ret_last]
        hypo_range = [hypo_first, hypo_last]


        do i = 1, size(f_array)
            do j = 1, 5
                if(f_array(i)%location(j:j) == "1" .and. j == 1)then
                    !Strat
                    conc(strat_first:strat_last, i) = f_array(i)%concs(1)
                elseif(f_array(i)%location(j:j) == "1" .and. j == 2)then
                    !epi
                    conc(epi_first:epi_last, i) = f_array(i)%concs(2)
                elseif(f_array(i)%location(j:j) == "1" .and. j == 3)then
                    !pap
                    conc(pap_first:pap_last, i) = f_array(i)%concs(3)
                elseif(f_array(i)%location(j:j) == "1" .and. j == 4)then
                    !ret
                    conc(ret_first:ret_last, i) = f_array(i)%concs(4)
                elseif(f_array(i)%location(j:j) == "1" .and. j == 5)then
                    !hypo
                    conc(hypo_first:hypo_last, i) = f_array(i)%concs(5)
                end if
            end do
        end do

        end subroutine gridset
end module gridset_mod