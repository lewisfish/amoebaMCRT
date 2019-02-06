module iarray
!
!  Contains all array var names.
!
    implicit none

    real, allocatable :: xface(:), yface(:), zface(:)
    real, allocatable :: rhokap(:), albedo(:), conc(:, :)
    ! real, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)
    real, allocatable :: refrac(:)


    ! real, allocatable :: nadh_array(:, :), nadh_fluro(:,:), nadh_cdf(:)
    ! real, allocatable :: fad_array(:, :), fad_fluro(:,:), fad_cdf(:)
    ! real, allocatable :: tyrosine_array(:, :), tyro_fluro(:,:), tyro_cdf(:)
    real, allocatable :: water_array(:,:), carotene_array(:,:), bilirubin_array(:,:), fat_array(:,:)
    real, allocatable :: oxy_hb_array(:, :), deoxy_hb_array(:, :)
end module iarray
