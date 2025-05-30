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
         b% other_adjust_mdots => L2_mdot

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
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

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

         print *, "accretion mode", b% accretion_mode
         print *, "use other j", b% use_other_accreted_material_j
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


! other_accreted_material_j
subroutine disk_accreted_material_j(binary_id, ierr)
  use binary_def, only : binary_info, binary_ptr
  use const_def, only: dp, standard_cgrav
  integer, intent(in) :: binary_id
  integer, intent(out) :: ierr
  real(dp) :: mdot ! mass lost to RLOF by the donor
  real(dp) :: qratio, min_r ! needed for reimplementation of UB76 fit to LS75
  type (binary_info), pointer :: b
  ierr = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if

  mdot = b% mtransfer_rate ! g/s
  print *, "in disk accreted j"
  print *, "mdot=", mdot/(Msun/secyer)

  ! reimplement Ulrich & Burger 1976 's fit to Lubow & Shu 1975
  ! as in de Mink+13, copying code from $MESA_DIR/binary/private/binary_mdot.f90

  qratio = b% m(b% a_i) / b% m(b% d_i)
  qratio = min(max(qratio,0.0667d0),15d0)
  min_r = 0.0425d0*b% separation*pow(qratio+qratio*qratio, 0.25d0)
  print *, "min_r", min_r/Rsun

  if (b% r(b% a_i) < min_r) then
     b% accretion_mode = 2 ! means there is a disk
     print *, "There is a disk!"
     ! TODO: implement once we have mdot(jdot) relation from
     ! Paczinsky's 1991 paper
     print *, "No AM accretion!"
     b% acc_am_div_kep_am = 0.0d0
  else
     b% accretion_mode = 1
     b% s_accretor% accreted_material_j = &
          sqrt(standard_cgrav * b% m(b% a_i) * 1.7d0*min_r)
     b% acc_am_div_kep_am = b% s_accretor% accreted_material_j / &
          sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i))
  end if

end subroutine disk_accreted_material_j


      end module run_binary_extras
