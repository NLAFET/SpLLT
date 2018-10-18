!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
module spllt_tree_stat_mod
  implicit none

  type spllt_tree_stat_t
    integer           :: max_dep                ! max #dep of a task
    integer           :: nblk_kdep              ! #block with more than k dep 
    !  (k = 2 by default)
  end type spllt_tree_stat_t

contains

  subroutine spllt_init_tree_stat(tree_stat)
    type(spllt_tree_stat_t), intent(inout)  :: tree_stat

    tree_stat%nblk_kdep = 0
    tree_stat%max_dep   = 0

  end subroutine spllt_init_tree_stat



  subroutine spllt_update_tree_stat(tree_stat, ndep)
    use spllt_data_mod
    type(spllt_tree_stat_t),  intent(inout) :: tree_stat
    integer,                  intent(in)    :: ndep   ! #dep of the block

    tree_stat%nblk_kdep = tree_stat%nblk_kdep + &
      merge(1, 0, ndep .gt. k_dep)
    tree_stat%max_dep   = merge(ndep, tree_stat%max_dep, &
      ndep .gt. tree_stat%max_dep)

  end subroutine spllt_update_tree_stat



  subroutine print_tree_stat(tree_stat, msg)
    use spllt_data_mod
    type(spllt_tree_stat_t),    intent(in)  :: tree_stat
    character(len=*), optional, intent(in)  :: msg

    if(present(msg)) print *, msg
    print '(a, i9)', "max #dep of a blk   : ", tree_stat%max_dep
    print '(a, i1, a, i9)', "#blk with #dep>", k_dep, &
      "    : ", tree_stat%nblk_kdep

  end subroutine print_tree_stat


end module spllt_tree_stat_mod
