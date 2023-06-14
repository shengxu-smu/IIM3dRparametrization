c----------------------------------------------------------------------
c
      subroutine singular_call

      call singular_force

      call jc_rhs
      call jc_first
      call jc_second

      call euler_link
      call jc_velocity

      call correction_interpolate
      call correction_strain

      call u_interpolate
      call u_strain
      call udu_surface

      call v_interpolate
      call v_strain
      call vdv_surface

      call w_interpolate
      call w_strain
      call wdw_surface        

      call jc_pressure
      call correction_difference

      return
      end


c-----------------------------------------------------------------------
