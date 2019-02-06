module opt_prop

    implicit none

    real    :: mua, mus, g2, hgg, kappa, wave, n1, n2, wave_incd
    real    :: opt_array(5, 2) !array of optical properties at wave_incd. mua, mus
    real    :: blood_ret, blood_pap, blood_hypo
    integer :: strat_range(2), epi_range(2), pap_range(2), ret_range(2), hypo_range(2) 
end module opt_prop
