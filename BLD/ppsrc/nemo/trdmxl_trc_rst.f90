# 1 "trdmxl_trc_rst.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "trdmxl_trc_rst.F90"
MODULE trdmxl_trc_rst
   !!======================================================================
   !!                       ***  MODULE  trdmxl_rst  ***
   !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
   !!======================================================================
   !! History :  9.0  ! 07-03 (C. Deltel) Original code
   !!----------------------------------------------------------------------
  
# 189 "trdmxl_trc_rst.F90"
  !!=================================================================================
  !!                       ***  MODULE  trdmxl_rst  ***
  !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
  !!=================================================================================
CONTAINS
  SUBROUTINE trd_mxl_trc_rst_opn( kt )
    IMPLICIT NONE
    INTEGER, INTENT( in ) :: kt
    WRITE(*,*) 'trd_mxl_trc_rst_opn: You should not have seen this print! error?', kt
  END SUBROUTINE trd_mxl_trc_rst_opn
  SUBROUTINE trd_mxl_trc_rst_write( kt )           !  No ML diags ==> empty routine
    IMPLICIT NONE
    INTEGER, INTENT( in ) :: kt
    WRITE(*,*) 'trd_mxl_trc_rst_wri: You should not have seen this print! error?', kt
  END SUBROUTINE trd_mxl_trc_rst_write
  SUBROUTINE trd_mxl_trc_rst_read                  !  No ML Diags ==> empty routine
    IMPLICIT NONE
  END SUBROUTINE trd_mxl_trc_rst_read


  !!=================================================================================
END MODULE trdmxl_trc_rst
