!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#define CP_L_LESS_Q(el1,el2) ( (.not.el1).and.el2 )


! *****************************************************************************
!> \brief various utilities that regard array of different kinds:
!>      output, allocation,...
!> \note
!>     ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _
!>    |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!>    | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!>    | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!>    |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>                                                        |___/
!>      ____ _                  ___                              _ _       _       _
!>     / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _
!>    | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!>    | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!>     \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!>                                                                                   |___/
!>     _____ _     _       _____ _ _      _
!>    |_   _| |__ (_)___  |  ___(_) | ___| |
!>      | | | '_ \| / __| | |_  | | |/ _ \ |
!>      | | | | | | \__ \ |  _| | | |  __/_|
!>      |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>
!>      This is a template
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify the .template and .instantition instead
!>
!>      maybe it is not a good idea mixing output and memeory utils...
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_array_utils_logical
  USE cp_log_handling,                 ONLY: cp_to_string
! dummy line needed for coverage broadcasting
  USE machine,                         ONLY: m_flush
#include "../base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_array_utils_logical'

  ! types
  PUBLIC :: cp_1d_logical_p_type, cp_2d_logical_p_type, cp_3d_logical_p_type

  ! generic interfaces
  PUBLIC :: cp_guarantee_size

  ! the underlying procedures
  PUBLIC :: cp_1d_logical_guarantee_size, &
       cp_1d_logical_write, cp_2d_logical_write,&
       cp_2d_logical_guarantee_size,&
       cp_1d_logical_bsearch

  INTERFACE cp_guarantee_size
     MODULE PROCEDURE cp_1d_logical_guarantee_size,&
          cp_2d_logical_guarantee_size
  END INTERFACE

!***

! *****************************************************************************
!> \brief represent a pointer to a 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE cp_1d_logical_p_type
     LOGICAL, DIMENSION(:), POINTER :: array
  END TYPE cp_1d_logical_p_type

! *****************************************************************************
!> \brief represent a pointer to a 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE cp_2d_logical_p_type
     LOGICAL, DIMENSION(:,:), POINTER :: array
  END TYPE cp_2d_logical_p_type

! *****************************************************************************
!> \brief represent a pointer to a 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE cp_3d_logical_p_type
     LOGICAL, DIMENSION(:,:,:), POINTER :: array
  END TYPE cp_3d_logical_p_type

CONTAINS

! *****************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated paretized list
! *****************************************************************************
SUBROUTINE cp_1d_logical_write(array, unit_nr, el_format)
    LOGICAL, INTENT(in)                      :: array(:)
    INTEGER, INTENT(in)                      :: unit_nr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: el_format

    CHARACTER(len=*), PARAMETER :: defaultFormat = "(l1)", &
      routineN = 'cp_1d_logical_write', routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iostat

  WRITE(unit=unit_nr,fmt="('( ')",advance="no", iostat=iostat)
  CPASSERT(iostat==0)
  IF (PRESENT(el_format)) THEN
     DO i=1,SIZE(array)-1
        WRITE(unit=unit_nr,fmt=el_format,advance="no") array(i)
        IF (MOD(i,5).EQ.0) THEN  ! only a few elements per line
           WRITE(unit=unit_nr,fmt="(',')")
        ELSE
           WRITE(unit=unit_nr,fmt="(',')",advance="no")
        ENDIF
     END DO
     IF (SIZE(array)>0) &
          WRITE(unit=unit_nr,fmt=el_format,advance="no") array(SIZE(array))
  ELSE
     DO i=1,SIZE(array)-1
        WRITE(unit=unit_nr,fmt=defaultFormat,advance="no") array(i)
        IF (MOD(i,5).EQ.0) THEN  ! only a few elements per line
           WRITE(unit=unit_nr,fmt="(',')")
        ELSE
           WRITE(unit=unit_nr,fmt="(',')",advance="no")
        ENDIF
     END DO
     IF (SIZE(array)>0) &
          WRITE(unit=unit_nr,fmt=defaultFormat,advance="no") array(SIZE(array))
  END IF
  WRITE(unit=unit_nr,fmt="(' )')")
  CALL m_flush(unit_nr)

END SUBROUTINE cp_1d_logical_write

! *****************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated parentized list
! *****************************************************************************
SUBROUTINE cp_2d_logical_write(array, unit_nr, el_format)
    LOGICAL, INTENT(in)                      :: array(:,:)
    INTEGER, INTENT(in)                      :: unit_nr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: el_format

    CHARACTER(len=*), PARAMETER :: defaultFormat = "(l1)", &
      routineN = 'cp_2d_logical_write', routineP = moduleN//':'//routineN

    CHARACTER(len=10)                        :: nRiga
    CHARACTER(len=200)                       :: fmtstr
    INTEGER                                  :: i, iostat

  nRiga=cp_to_string(SIZE(array,2))
  DO i=1,SIZE(array,1)
     IF (PRESENT(el_format)) THEN
        fmtstr='(" ",'//nRiga//el_format//')'
        WRITE(unit=unit_nr,fmt=fmtstr,iostat=iostat) array(i,:)
     ELSE
        fmtstr = '(" ",'//nRiga//defaultFormat//')'
        WRITE(unit=unit_nr,fmt=fmtstr,iostat=iostat) array(i,:)
     END IF
     CPASSERT(iostat==0)
  END DO
  CALL m_flush(unit_nr)
END SUBROUTINE cp_2d_logical_write

! *****************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n the wanted size
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! *****************************************************************************
  SUBROUTINE cp_1d_logical_guarantee_size(array, n)
    LOGICAL, POINTER                         :: array(:)
    INTEGER, INTENT(in)                      :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_1d_logical_guarantee_size', &
      routineP = moduleN//':'//routineN

    CPASSERT(n>=0)
    IF (ASSOCIATED(array)) THEN
       IF (SIZE(array) /= n) THEN
          CPWARN('size has changed')
          DEALLOCATE(array)
       END IF
    END IF
    IF (.NOT.ASSOCIATED(array)) THEN
       ALLOCATE(array(n))
    END IF
  END SUBROUTINE cp_1d_logical_guarantee_size

! *****************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n_rows the wanted number of rows
!> \param n_cols the wanted number of cols
!> \par History
!>      5.2001 first version [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! *****************************************************************************
  SUBROUTINE cp_2d_logical_guarantee_size(array, n_rows, n_cols)
    LOGICAL, POINTER                         :: array(:,:)
    INTEGER, INTENT(in)                      :: n_rows, n_cols

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_2d_logical_guarantee_size', &
      routineP = moduleN//':'//routineN

    CPASSERT(n_cols>=0)
    CPASSERT(n_rows>=0)
    IF (ASSOCIATED(array)) THEN
       IF (SIZE(array,1) /= n_rows .OR. SIZE(array,2) /= n_cols) THEN
          CPWARN('size has changed')
          DEALLOCATE(array)
       END IF
    END IF
    IF (.NOT.ASSOCIATED(array)) THEN
       ALLOCATE(array(n_rows,n_cols))
    END IF
  END SUBROUTINE cp_2d_logical_guarantee_size

! *****************************************************************************
!> \brief returns the index at which the element el should be inserted in the
!>      array to keep it ordered (array(i)>=el).
!>      If the element is bigger than all the elements in the array returns
!>      the last index+1.
!> \param array the array to search
!> \param el the element to look for
!> \param l_index the lower index for binary search (defaults to 1)
!> \param u_index the upper index for binary search (defaults to size(array))
!> \retval res ...
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the array should be ordered in growing order
! *****************************************************************************
FUNCTION cp_1d_logical_bsearch(array, el, l_index, u_index)&
     RESULT(res)
    LOGICAL, INTENT(in)                      :: array(:), el
    INTEGER, INTENT(in), OPTIONAL            :: l_index, u_index
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_1d_logical_bsearch', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aindex, lindex, uindex

  lindex=1
  uindex=SIZE(array)
  IF (PRESENT(l_index)) lindex=l_index
  IF (PRESENT(u_index)) uindex=u_index
  DO WHILE (lindex <= uindex)
     aindex=(lindex+uindex)/2
     IF (CP_L_LESS_Q(array(aindex),el)) THEN
        lindex=aindex+1
     ELSE
        uindex=aindex-1
     END IF
  END DO
  res=lindex
END FUNCTION cp_1d_logical_bsearch

  ! template def put here so that line numbers in template and derived
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,USE,defaultFormatType1,lessQ,defines)]
! ARGS:
!  USE = "! dummy line needed for coverage broadcasting"
!  defaultFormatType1 = ""(l1)""
!  defines = 
!    "#define CP_L_LESS_Q(el1,el2) ( (.not.el1).and.el2 )
!     "
!  lessQ = "CP_L_LESS_Q"
!  nametype1 = "logical"
!  type1 = "logical"


END MODULE cp_array_utils_logical
