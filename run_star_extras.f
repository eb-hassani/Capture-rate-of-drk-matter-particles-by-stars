! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
      
            subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 7
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      use math_lib, only: safe_log10
      use chem_def, only: ih1, ihe3, ihe4, ic12, in14, io16, ine20, img24

      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp), parameter :: frac = 0.90
      integer :: i
      real(dp) :: edot, edot_partial
      !!!!!!!!!!!!!!!   star parameters
      REAL(dp) ::  avogad = 6.02214179D23   ! avogadro number
      REAL(dp) ::  At_h1    = 1            ! Atomic number of h1
      REAL(dp) ::  At_he4   = 2            ! Atomic number of he4_cap_rate
      REAL(dp) ::  At_c12   = 6            ! Atomic number of c12
      REAL(dp) ::  At_n14   = 7            ! Atomic number of n14
      REAL(dp) ::  At_o16   = 8            ! Atomic number of o16
      REAL(dp) ::  At_ne20  = 10           ! Atomic number of ne20_cap_rate
      REAL(dp) ::  At_mg24  = 12           ! Atomic number of mg24

      !!!!!!!!!!!!!!!   Other parameters
      REAL(dp) :: part3_h1
      REAL(dp) :: part3_he4
      REAL(dp) :: part3_c12
      REAL(dp) :: part3_n14
      REAL(dp) :: part3_o16
      REAL(dp) :: part3_ne20
      REAL(dp) :: part3_mg24
     
      ! part_3_h1
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      edot = dot_product(s% dm(1:s% nz), s% eps_nuc(1:s% nz))

      edot_partial = 0
      do i = s% nz, 1, -1
      edot_partial = edot_partial + s% dm(i) * s% eps_nuc(i)
      if (edot_partial .ge. (frac * edot)) exit
      end do
      
      part3_h1 = 0
      do i = s% nz, 2, -1
      part3_h1 = part3_h1 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(ih1),i)) &
      / (At_h1))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_h1 = part3_h1 + part3_h1 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(ih1),1)) &
      / (At_h1))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))
   
      ! part_3_he4      
      part3_he4 = 0
      do i = s% nz, 2, -1
      part3_he4 = part3_he4 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(ihe4),i)) &
      / (At_he4))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_he4 = part3_he4 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(ihe4),1)) &
      / (At_he4))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))
     
      ! part_3_c12
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      edot = dot_product(s% dm(1:s% nz), s% eps_nuc(1:s% nz))

      edot_partial = 0
      do i = s% nz, 1, -1
      edot_partial = edot_partial + s% dm(i) * s% eps_nuc(i)
      if (edot_partial .ge. (frac * edot)) exit
      end do
      
      part3_c12 = 0
      do i = s% nz, 2, -1
      part3_c12 = part3_c12 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(ic12),i)) &
      / (At_c12))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_c12 = part3_c12 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(ic12),1)) &
      / (At_c12))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))

      ! part_3_n14      
      part3_n14 = 0
      do i = s% nz, 2, -1
      part3_n14 = part3_n14 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(in14),i)) &
      / (At_n14))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_n14 = part3_n14 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(in14),1)) &
      / (At_n14))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))

      ! part_3_o16      
      part3_o16 = 0
      do i = s% nz, 2, -1
      part3_o16 = part3_o16 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(io16),i)) &
      / (At_o16))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_o16 = part3_o16 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(io16),1)) &
      / (At_o16))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))
      
      ! part_3_ne20      
      part3_ne20 = 0
      do i = s% nz, 2, -1
      part3_ne20 = part3_ne20 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(ine20),i)) &
      / (At_ne20))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_ne20 = part3_ne20 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(ine20),1)) &
      / (At_ne20))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))
    
      ! part_3_mg24      
      part3_mg24 = 0
      do i = s% nz, 2, -1
      part3_mg24 = part3_mg24 + ((avogad * (s% rho(i)) &
      *(s% xa(s% net_iso(img24),i)) &
      / (At_mg24))) * ((s% r(i))**2) * ((s% r(i-1)) - (s% r(i)))
      end do
      part3_mg24 = part3_mg24 + ((avogad * (s% rho(1)) &
      *(s% xa(s% net_iso(img24),1)) &
      / (At_mg24))) * ((s% r(1))**2) * ((s% r(1)) - (s% r(2)))
     
      ! column n(h1)*r^2
      names(1) = "n(h1)*r^2"
      vals(1)  =  part3_h1
      PRINT*, 'n(h1)*r^2=', vals(1)
      
      ! column n(he4)*r^2
      names(2) = "n(he4)*r^2"
      vals(2)  = part3_he4
      PRINT*, 'n(he4)*r^2=', vals(2)
      
      ! column n(c12)*r^2
      names(3) = "n(c12)*r^2"
      vals(3)  = part3_c12
      !PRINT*, 'n(c12)*r^2=', vals(3)
      
      ! column n(n14)*r^2
      names(4) = "n(n14)*r^2"
      vals(4)  = part3_n14
      !PRINT*, 'n(n14)*r^2=', vals(4) 
      
      ! column n(o16)*r^2
      names(5) = "n(o16)*r^2"
      vals(5)  = part3_o16
      !PRINT*, 'n(o16)*r^2=', vals(5)  
      
      ! column n(ne20)*r^2
      names(6) = "n(ne20)*r^2"
      vals(6)  = part3_ne20
      !PRINT*, 'n(ne20)*r^2=', vals(6)  
      
      ! column n(mg24)*r^2
      names(7) = "n(mg24)*r^2"
      vals(7)  = part3_mg24
      !PRINT*, 'n(mg24)*r^2=', vals(7) 

      ierr = 0
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      end module run_star_extras
