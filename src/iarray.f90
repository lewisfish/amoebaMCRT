module iarray
!
!  Contains all array var names.
!
    implicit none

    real, allocatable :: xface(:), yface(:), zface(:)
    real, allocatable :: rhokap(:), albedo(:), conc(:, :)
    real, allocatable :: refrac(:)
    real, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)

    real, allocatable :: water_array(:,:), carotene_array(:,:), bilirubin_array(:,:), fat_array(:,:)
    real, allocatable :: oxy_hb_array(:, :), deoxy_hb_array(:, :)
end module iarray
