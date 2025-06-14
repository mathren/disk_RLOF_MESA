! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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
      module run_binary_extras

      use star_lib
      use star_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib

      implicit none

      contains

      ! include custom functions for Paczynski 1991
      include 'binary_disk.inc'

      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         b% warn_binary_extra =.false.

         ! point to function defined in binary_disk.inc
         b% other_accreted_material_j => disk_accreted_material_j
         b% other_extra_jdot => L2_extra_jdot
         b% other_adjust_mdots => L2_adjust_mdots

      end subroutine extras_binary_controls

      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 5
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n), fL2, xL2, omega_orb, log_q
         real(dp), parameter  :: g_b = 1.0d0
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         names(1) = 'fL2'
         names(2) = 'mdot_L2'
         names(3) = 'extra_jdot'
         names(4) = 'jdot_acc'
         names(5) = 'jdot_L2'
         if (b% accretion_mode /= 2) then ! no disk
            vals(1) = 0.0d0
            vals(2) = 0.0d0
         else
            call eval_L2_mass_loss_fraction(b% s_donor% m(1)/Msun, &
                 b% s_accretor% m(1)/Msun, &
                 b% mtransfer_rate/(Msun/secyer), &
                 b% separation/Rsun, &
                 0.1d0, & ! disk alpha viscosity
                 1.3d0/2.4d0, &  ! full ionization
                 fL2, ierr)
            vals(1) = fL2
            vals(2) = (b% mtransfer_rate * fL2)/(Msun/secyer)
         end if
         vals(3) = b% extra_jdot
         vals(4) = b% s_accretor% accreted_material_j
         ! recalculate L2 loss, may differ from extra_jdot for AM not accreted when disk accretion
         log_q = log10(b% s_accretor% m(1)/ b% s_donor% m(1))
         ! position of L2 w.r.t. center of mass according to Lu et al. 23 fit in units of separation
         xL2 = 0.0756_dp * log_q**2 - 0.424_dp * abs(log_q) + 1.699_dp  ! xL2 = rL2/a
         omega_orb = 2*pi/b% period ! 1/sec
         vals(5) = ((b% mtransfer_rate * fL2)) * g_b * & ! amount of mass lost at L2, note that mtransfer_rate is negative
              omega_orb * pow2(b% separation) * pow2(xL2-b% s_accretor% m(1) / (b%s_accretor% m(1)+ b% s_donor% m(1)))
       end subroutine data_for_extra_binary_history_columns


      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         b% lxtra(1) = .false. ! flag for end of donor's main sequence
         b% lxtra(2) = .false. ! flag for beginning RLOF
         b% lxtra_old(2) = .false.

         if (b% x_ctrl(2) - b% x_ctrl(1) < 0) then
            print *, "To interpolate we need b% x_ctrl(2) > b% x_ctrl(1)"
            print *, "Please fix your inlist"
            STOP "interp values"
         end if

         print *, "accretion mode", b% accretion_mode

         extras_binary_startup = keep_going
      end function  extras_binary_startup

      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

      end function  extras_binary_start_step

      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         extras_binary_check_model = keep_going

      end function extras_binary_check_model


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         character (len=200) :: fname
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         extras_binary_finish_step = keep_going

         ! find beginning RLOF
         if ((b% lxtra(2) .eqv. .false.) .and. &
              (b% lxtra_old(2) .eqv. .false.)) then
            ! RLOF has not started before
            if (b% rl_relative_gap(b% d_i) > 0) then
               write(fname, fmt="(a20)") 'donor_onset_RLOF.mod'
               call star_write_model(b% star_ids(1), fname, ierr)
               write(fname, fmt="(a23)") 'accretor_onset_RLOF.mod'
               call star_write_model(b% star_ids(2), fname, ierr)
               b% lxtra(2) = .true.
               b% lxtra_old(2) = .true.
            end if
         end if

         ! find donor's TAMS
         if ((b% s_donor% xa(b% s_donor% net_iso(ih1), b% s_donor% nz) < 1d-4) .and. &
              (b% lxtra(1) .eqv. .false.)) then
            b% lxtra(1) =.true.
            b% xtra(1) = b% s_donor% r(1)
            print *, "saved donor radius at TAMS", b% xtra(1)
            write(fname, fmt="(a14)") 'donor_TAMS.mod'
            call star_write_model(b% star_ids(1), fname, ierr)
         end if

      end function extras_binary_finish_step

      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
      end subroutine extras_binary_after_evolve

      end module run_binary_extras
