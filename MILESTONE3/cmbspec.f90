program cmbspec
    use healpix_types
    use params
    use time_mod
    use rec_mod
    use evolution_mod
    implicit none

    ! Initialize time grids
    call initialize_time_mod
    call initialize_rec_mod
    call initialize_perturbation_eqns
    call integrate_perturbation_eqns

    ! Output to file desired quantities here
    call write_mk3

end program cmbspec
