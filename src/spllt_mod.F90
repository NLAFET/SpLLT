module spllt_mod
  use spllt_data_mod
#if defined(SPLLT_USE_STARPU)
  use iso_c_binding
  use starpu_f_mod
#elif defined(SPLLT_USE_PARSEC)
  use parsec_f08_interfaces
#endif
  use spllt_analyse_mod
  implicit none
  
  ! Read matrix in Matrix Market foramt  
  interface mm_read
     module procedure mm_double_read 
  end interface mm_read

  ! convert matrix from COO to CSC format
  interface coo_to_csc
     module procedure coo_to_csc_double
  end interface coo_to_csc

contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Initialize SpLLT.
  !>
  !> @param options User-supplied options.
  subroutine spllt_init(options)
#if defined(SPLLT_USE_STARPU)
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE) 
    use trace_mod
#endif
#elif defined(SPLLT_USE_PARSEC)
    use parsec_f08_interfaces
    use spllt_parsec_mod
#endif

#if defined(SPLLT_USE_GPU)
    use magma
#endif
    implicit none

    type(spllt_options) :: options
#if defined(SPLLT_USE_STARPU)
    integer :: start_starpuinit_t, stop_starpuinit_t, rate_starpuinit_t
    integer(c_int) :: ret
#endif

#if defined(SPLLT_USE_STARPU)

    call system_clock(start_starpuinit_t, rate_starpuinit_t)
    ! initialize starpu
    ret = starpu_f_init(options%ncpu)
    call system_clock(stop_starpuinit_t)
    write(*,'("[>] [spllt_test_mat] StarPU init time: ", es10.3, " s")') &
         &(stop_starpuinit_t - start_starpuinit_t)/real(rate_starpuinit_t)
    call starpu_f_fxt_start_profiling()

#elif defined(SPLLT_USE_OMP)

#if defined(SPLLT_OMP_TRACE)

    call trace_init(omp_get_num_threads())

#endif

#elif defined(SPLLT_USE_PARSEC)

    ctx = spllt_parsec_init(options%ncpu, nds, rank)
    write(*,'("Parsec init    nodes: ", i6, ", rank: ", i6)') nds, rank

#endif

#if defined(SPLLT_USE_GPU)

    call magmaf_init()

#endif

    return
  end subroutine spllt_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Shutdowns SpLLT.
  subroutine spllt_finalize()
#if defined(SPLLT_USE_STARPU)
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE) 
    use trace_mod
#endif
#elif defined(SPLLT_USE_PARSEC)
    use parsec_f08_interfaces
    use spllt_parsec_mod
#endif
    ! use spllt_factorization_task_mod 
    implicit none

    integer :: start_starpushutdown_t, stop_starpushutdown_t, rate_starpushutdown_t

#if defined(SPLLT_USE_GPU)
    call magmaf_finalize()
#endif

#if defined(SPLLT_USE_STARPU)

    call system_clock(start_starpushutdown_t, rate_starpushutdown_t)
    call starpu_f_shutdown()
    call system_clock(stop_starpushutdown_t)
    write(*,'("[>] [spllt_test_mat] StarPU shutdown time: ", es10.3, " s")') &
         &(stop_starpushutdown_t - start_starpushutdown_t)/real(rate_starpushutdown_t)

#elif defined(SPLLT_USE_OMP)

#if defined(SPLLT_OMP_TRACE) 
    call trace_log_dump_paje('trace')
#endif

#elif defined(SPLLT_USE_PARSEC)

    ! call dague_fini(ctx)
    call spllt_parsec_fini(ctx)

#endif
    return
  end subroutine spllt_finalize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Factorizes the input matrix.
  subroutine spllt_factor(akeep, fkeep, options, val, info)
#if defined(SPLLT_USE_STF) || defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)
    use spllt_stf_mod
#elif defined(SPLLT_USE_PARSEC)
    use spllt_ptg_mod
#endif    
    implicit none

    type(spllt_akeep), intent(in)    :: akeep
    type(spllt_fkeep), intent(inout) :: fkeep
    type(spllt_options), intent(in)  :: options ! User-supplied options
    real(wp), intent(in)             :: val(:) ! Matrix values
    type(spllt_inform), intent(out)  :: info 


#if defined(SPLLT_USE_STF) || defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)

    ! Call the STF-based factorize routine.
    call spllt_stf_factorize(akeep, fkeep, options, val, info)

#elif defined(SPLLT_USE_PARSEC)

    ! Call the PTG-based factorize routine.
    call spllt_ptg_factorize(akeep, fkeep, options, val, info)

#endif

  end subroutine spllt_factor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Wait for completion of previously submitted tasks.
  subroutine spllt_wait()
    use spllt_data_mod
    implicit none
#if defined(SPLLT_USE_STARPU)
    call starpu_f_task_wait_for_all()
#elif defined(SPLLT_USE_OMP)
    !$omp taskwait
#elif defined(SPLLT_USE_PARSEC)
    call parsec_context_wait(ctx)
#endif
  end subroutine spllt_wait

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Prints the assemlby tree.  
  !>
  !> @param akeep Symbolic factorization data.
  !> @param fkeep Factorization data.
  !> @param options User-supplied options. 
  subroutine spllt_print_atree(akeep, fkeep, options)
    use spllt_data_mod
    implicit none    

    type(spllt_akeep), intent(in)  :: akeep
    type(spllt_fkeep), target, intent(in) :: fkeep
    type(spllt_options), intent(in) :: options

    integer :: snode, num_nodes
    type(spllt_node), pointer     :: node ! node in the atree
    integer :: m, n

    num_nodes = fkeep%info%num_nodes

    open(2, file="atree.dot")

    write(2, '("graph atree {")')
    write(2, '("node [")')
    write(2, '("style=filled")')
    write(2, '("]")')
    
    do snode=1,num_nodes

       node => fkeep%nodes(snode)
       m = size(node%index)
       n = node%en - node%sa + 1

       write(2, '(i10)', advance="no")snode
       write(2, '(" ")', advance="no")
       write(2, '("[")', advance="no")
       if (options%prune_tree) then
          if (akeep%small(snode) .eq. 1) then
             write(2, '("fillcolor=grey ")', advance="no")
          else
             write(2, '("fillcolor=white ")', advance="no")
          endif
       else
          write(2, '("fillcolor=white ")', advance="no")
       endif
       write(2, '("label=""")', advance="no")  
       write(2, '("node:", i5,"\n")', advance="no")snode
       write(2, '("m:", i5,"\n")', advance="no")m
       write(2, '("n:", i5,"\n")', advance="no")n
       if (options%prune_tree) then
          write(2, '("small:", i5,"\n")', advance="no")akeep%small(snode)
          write(2, '("weight:", f5.1,"\n")', advance="no") &
               & 100.0 * (real(akeep%weight(snode), kind(1d0)) / real(akeep%weight(num_nodes+1), kind(1d0)))
       end if
       write(2, '("""")', advance="no")         
       write(2, '("]")', advance="no")
       write(2, '(" ")')

       if(fkeep%nodes(snode)%parent .ne. -1) write(2, '(i10, "--", i10)')fkeep%nodes(snode)%parent, snode
    end do

    write(2, '("}")')

    close(2)
    
    return
  end subroutine spllt_print_atree



  function spllt_dep_in_add(dep_in_arr, id_jk, id_ik, csrc, rsrc, sync)
    implicit none

    type(spllt_dep_upd_in), dimension(:), pointer :: dep_in_arr
    integer(long) :: id_jk, id_ik
    integer :: csrc, rsrc
    integer :: spllt_dep_in_add
    logical, optional :: sync

    integer :: sz
    logical :: s
    type(spllt_dep_upd_in), dimension(:), pointer :: new_dep_in_arr => null()

    if (present(sync)) then
       s = sync
    else
       s = .true.
    end if
    
    if (.not. associated(dep_in_arr)) then
       allocate(dep_in_arr(1))
       dep_in_arr(1)%id_jk = id_jk
       dep_in_arr(1)%id_ik = id_ik
       dep_in_arr(1)%csrc = csrc
       dep_in_arr(1)%rsrc = rsrc
       dep_in_arr(1)%sync = s
    else
       sz = size(dep_in_arr)
       allocate(new_dep_in_arr(sz+1))
       new_dep_in_arr(1:sz) = dep_in_arr(1:sz)       
       new_dep_in_arr(sz+1)%id_jk = id_jk
       new_dep_in_arr(sz+1)%id_ik = id_ik
       new_dep_in_arr(sz+1)%csrc = csrc
       new_dep_in_arr(sz+1)%rsrc = rsrc
       new_dep_in_arr(sz+1)%sync = s
       deallocate(dep_in_arr)
       dep_in_arr => new_dep_in_arr
    end if

    spllt_dep_in_add = size(dep_in_arr)

  end function spllt_dep_in_add

  function spllt_dep_out_add(dep_out_arr, id_ij, flow)
    implicit none

    type(spllt_dep_upd_out), dimension(:), pointer :: dep_out_arr
    integer(long) :: id_ij
    integer :: flow
    integer :: spllt_dep_out_add

    integer :: sz
    type(spllt_dep_upd_out), dimension(:), pointer :: new_dep_out_arr => null()
    
    if (.not. associated(dep_out_arr)) then
       allocate(dep_out_arr(1))
       dep_out_arr(1)%id_ij = id_ij
       dep_out_arr(1)%flow = flow
    else
       sz = size(dep_out_arr)
       allocate(new_dep_out_arr(sz+1))
       new_dep_out_arr(1:sz) = dep_out_arr(1:sz)       
       new_dep_out_arr(sz+1)%id_ij = id_ij
       new_dep_out_arr(sz+1)%flow = flow
       deallocate(dep_out_arr)
       dep_out_arr => new_dep_out_arr
    end if

    spllt_dep_out_add = size(dep_out_arr)

  end function spllt_dep_out_add
  


  subroutine spllt_parse_args(options, matfile, nrhs)
    implicit none

    type(spllt_options), intent(inout) :: options
    character(len=200), intent(inout) :: matfile
    integer, intent(inout) :: nrhs ! Number of right-hand side
    
    integer :: argnum, narg
    character(len=200) :: argval

    narg = command_argument_count()
    argnum = 1
    do while(argnum <= narg)
       call get_command_argument(argnum, argval)
       argnum = argnum + 1
       select case(argval)
       case("--mat")
          ! input matrix in Rutherford Boeing format
          options%fmt = 'csc'
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) matfile
       case("--nb")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%nb
       case("--ncpu")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%ncpu
       case("--nemin")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%nemin
       case("--prune-tree")
          options%prune_tree = .true.
       case("--no-prune-tree")
          options%prune_tree = .false.
       case("--mm")
          ! input matrix in Matrix Market format
          options%fmt = 'coo'
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) matfile
       case("--nrhs")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) nrhs
          print *, 'solving for', nrhs, 'right-hand sides'
       case("--nb-range")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%nb_min
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%nb_max
       case("--nrhs-range")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%nrhs_min
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%nrhs_max
          print *, 'solving from ', options%nrhs_min, &
            ' to ', options%nrhs_max, ' RHS'
        case("--nb-linear")
          options%nb_linear_comp = .true.
        case("--no-nb-linear")
          options%nb_linear_comp = .false.
        case("--nrhs-linear")
          options%nrhs_linear_comp = .true.
        case("--no-nrhs-linear")
          options%nrhs_linear_comp = .false.

       case default
          write(*,'("Unrecognised command line argument: ", a20)') argval
       end select
    end do

  end subroutine spllt_parse_args



  ! Read matrix in Matrix Market foramt
  subroutine mm_double_read(matfile, m, n, nnz, indx, jndx, val, info)
    use spral_random, only : random_state, random_real
    implicit none
    
    character(len=*), intent(in)   :: matfile ! matrix file name
    integer, intent(out) :: m
    integer, intent(out) :: n
    integer, intent(out) :: nnz
    integer, dimension(:), allocatable, intent(out) :: indx ! row index array
    integer, dimension(:), allocatable, intent(out) :: jndx ! column index array
    real(wp), dimension(:), allocatable, intent(out) :: val ! entry array
    integer, intent(out) :: info

    integer :: i
    character(len=20) :: rep ! Representation e.g. Coordinate
    character(len=20) :: field, symm, typ, fmt
    logical :: values ! values provided
    type(random_state) :: state  
    integer :: err ! error code

    info = 0

    ! initialize matrix description
    rep   = ''
    field = ''
    symm  = ''
    typ   = ''
    fmt   = ''
  
    open(4,file=matfile, status='OLD', action='READ', iostat=err)
    if (err.gt.0) then
       goto 100
    end if

    read(4,*)rep,typ,fmt,field,symm

    read(4,*)rep
    do
       if(rep(1:1) .ne. '%') exit
       read(4,*)rep
    end do

    backspace(4)

    read(4,*)m,n,nnz

    values = field .ne. 'pattern'
    ! print *, values
    allocate(indx(nnz), jndx(nnz), val(nnz))

    if(values) then
       do i=1, nnz
          read(4,*)indx(i), jndx(i), val(i)
       end do
    else
       ! make up values
       do i=1, nnz
          read(4,*)indx(i), jndx(i)
          val(i) = random_real(state, .false.)
       end do
    end if

    close(4)

100 continue

    info = err
    return

  end subroutine mm_double_read



  ! convert matrix from CSC to CSR format
  subroutine csc_to_csr_double(m, n, nnz, colPtr, rowInd, val_in,       &
    rowPtr, colInd, val_out, info)
    use utils_mod
    implicit none

    integer,                intent(in)  :: m         ! Number of rows
    integer,                intent(in)  :: n         ! Number of columns
    integer,                intent(in)  :: nnz       ! Number of entries
    integer,                intent(in)  :: colPtr(:) ! Elements that index rowInd array
    integer,                intent(in)  :: rowInd(:) ! Stores the row index
    real(wp),               intent(in)  :: val_in(:) ! Entry array in CSC format
    integer,  allocatable,  intent(out) :: rowPtr(:) ! Elements that index colInd array
    integer,  allocatable,  intent(out) :: colInd(:) ! Stores the col index
    real(wp), allocatable,  intent(out) :: val_out(:)! Entry array in CSR format
    integer,                intent(out) :: info      ! status

    integer, allocatable :: offset(:) ! counter
    integer :: i, j, k
    integer :: st

   !info = 0 ! init to success

   !! allocate ptr 
   !allocate(rowPtr(m + 1))
   !allocate(offset(m))
   ! 
   !rowPtr = 0
   !offset = 0

   !allocate(colInd(nnz),   stat=st)
   !info = info + st
   !allocate(val_out(nnz),  stat=st)
   !info = info + st

    ! Count
   !do i = 1, n
   !  do j = colPtr(i), colPtr(i+1) - 1
   !    k = rowInd(j)
   !    rowPtr(k+1) = rowPtr(k+1) + 1
   !  end do
   !end do
   !
   !! Sum
   !rowPtr(1) = 1
   !do j = 1, m
   !  rowPtr(j + 1) = rowPtr(j + 1) + rowPtr(j)
   !end do

   !! Store
   !do i = 1, n
   !  do j = colPtr(i), colPtr(i+1) - 1
   !    k = rowInd(j)
   !    colInd(rowPtr(k) + offset(k))   = i
   !    val_out(rowPtr(k) + offset(k))  = val_in(j)
   !    offset(k)                       = offset(k) + 1
   !  end do
   !end do

   !deallocate(offset)

  end subroutine csc_to_csr_double



  ! convert matrix from COO to CSC format
  subroutine coo_to_csc_double(m, n, nnz, indx_in, jndx, val_in, & 
       indx_out, val_out, ptr, info)
    implicit none

    integer, intent(in) :: m ! number of rows
    integer, intent(in) :: n ! number of columns
    integer, intent(in) :: nnz ! number of entries
    integer, dimension(nnz), intent(in) :: indx_in ! row index in COO format
    integer, dimension(nnz), intent(in) :: jndx ! row index in COO format
    real(wp), dimension(nnz), intent(in) :: val_in ! entry array in COO format
    integer, dimension(:), allocatable, intent(out) :: indx_out ! row index in CSC format
    real(wp), dimension(:), allocatable, intent(out) :: val_out ! entry array in CSC format
    integer, dimension(:), allocatable, intent(out) :: ptr ! pointer array in CSC format
    integer, intent(out) :: info ! status

    integer, dimension(:), allocatable :: work ! counter
    integer :: i, j, k

    info = 0 ! init to success

    ! allocate ptr 
    allocate(ptr(n+1))
    ! allocate temporary array
    allocate(work(n))
     
    work = 0 ! init

    ! count number of entries per columns
    do k = 1, nnz
       i = indx_in(k)
       j = jndx(k)
       work(j) = work(j)+1 
    end do

    ! create ptr array, such that row indexe for column j are in
    ! indx(ptr(i):ptr(i+1)-1) and entries in val(ptr(i):ptr(i+1)-1)
    ptr(1) = 1
    do k = 2, n+1
       ptr(k) = ptr(k-1) + work(k-1)
    end do

    ! create new indx and val array
    allocate(indx_out(nnz), val_out(nnz))
    work = 0
    do k = 1, nnz
       i = indx_in(k)
       j = jndx(k)
       indx_out(ptr(j)+work(j)) = i
       val_out(ptr(j)+work(j)) = val_in(k)
       work(j) = work(j)+1       
    end do

    ! deallocate temporary array
    deallocate(work)

  end subroutine coo_to_csc_double

end module spllt_mod
