!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
!!$omp default(none)                              &
!$omp firstprivate(p_rhs)                         &
!$omp firstprivate(nrhs, n)                       &
!$omp firstprivate(sub_task_manager)              &
!$omp firstprivate(p_bc, blk_en)                  &
!$omp firstprivate(en, sa, ndep_lvl, chunk)       &
!$omp shared(fkeep)                               &
!$omp depend(inout: p_bc(blk_en))                 &
!$omp private(i, trace_id)
