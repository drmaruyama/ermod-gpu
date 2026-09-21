! -*- F90 -*-
! ERmod - Energy Representation Module
! Copyright (C) 2000- The ERmod authors
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

! =====================================================================
! AMD/ROCm port: cuFFT -> hipFFT (via hipfort). hipFFT's API mirrors
! cuFFT's, so the cufftXxx/CUFFT_XXX -> hipfftXxx/HIPFFT_XXX rename is
! mostly mechanical, with one notable difference: cuFFT's Fortran plan
! handle is a plain "integer", while hipfort's is "type(c_ptr)"
! (iso_c_binding) -- fft_handle's "plan" component is typed accordingly.
! hipfort's hipfftExec* interfaces take the array directly (using
! c_loc internally), which requires the "target" attribute on the
! actual array argument, so it has been added to the in/out dummy
! arguments here.
!
! The device-pointer handoff (!$omp target data ... use_device_addr)
! is unchanged from the cuFFT version; in/out (the actual arguments
! are recpcal.F90's cnvslt/rcpslt) are still expected to already be
! resident on the device via recpcal_init's target enter data.
!
! Building this requires hipfort's .mod files and library
! (hipfort-amdgcn etc.) plus the hipFFT/HIP runtime libraries; see
! configure.ac's --with-hipfort option.
! =====================================================================

module fft_iface
  use precision_kinds, only: wp
  use iso_c_binding, only: c_ptr, c_null_ptr
  use hipfort_hipfft
  implicit none

  integer :: fftsize(3)

  type fft_handle
     type(c_ptr) :: plan = c_null_ptr
     ! transforms per execution: 1 for hipfftPlan3d plans, nbatch for
     ! the batched hipfftPlanMany plans below
     integer :: nbatch = 1
  end type fft_handle

contains 

  subroutine fft_set_size(fftsize_in)
    integer, intent(in) :: fftsize_in(3)
    fftsize(:) = fftsize_in(:)
  end subroutine fft_set_size

  ! Check the return status of a hipFFT call and abort (on every rank, under
  ! MPI) with a diagnostic message if it did not succeed. Every hipfftPlan3d/
  ! hipfftExec*/hipfftDestroy call below used to discard its status entirely
  ! (in the original cuFFT version), so a GPU/library-level failure (e.g. an
  ! unsupported transform size, or running out of device memory) would
  ! silently continue with whatever garbage was left in the output array,
  ! instead of stopping right away with a message that says what went wrong.
  !
  ! [location] should identify the failing call (e.g. "hipfftPlan3d (fft_init_rtc)").
  ! hipFFT mirrors cuFFT's status codes; see:
  ! https://rocm.docs.amd.com/projects/hipFFT/en/latest/
  ! They are spelled out here as literals rather than named constants,
  ! since not every "hipfort_hipfft" module version exposes all of them by name.
  subroutine check_hipfft_status(stat, location)
    use engmain, only: stdout
    use mpiproc, only: mpi_abend
    implicit none
    integer, intent(in) :: stat
    character(len=*), intent(in) :: location
    character(len=64) :: msg

    if (stat == HIPFFT_SUCCESS) return

    select case(stat)
    case(1);  msg = "invalid plan handle (HIPFFT_INVALID_PLAN)"
    case(2);  msg = "failed to allocate GPU or CPU memory (HIPFFT_ALLOC_FAILED)"
    case(3);  msg = "invalid transform type (HIPFFT_INVALID_TYPE)"
    case(4);  msg = "invalid pointer or parameter (HIPFFT_INVALID_VALUE)"
    case(5);  msg = "internal driver error (HIPFFT_INTERNAL_ERROR)"
    case(6);  msg = "failed to execute the FFT on the GPU (HIPFFT_EXEC_FAILED)"
    case(7);  msg = "the hipFFT library failed to initialize (HIPFFT_SETUP_FAILED)"
    case(8);  msg = "invalid transform size (HIPFFT_INVALID_SIZE)"
    case(9);  msg = "unaligned data (HIPFFT_UNALIGNED_DATA)"
    case(10); msg = "missing parameters in call (HIPFFT_INCOMPLETE_PARAMETER_LIST)"
    case(11); msg = "plan executed on a different GPU than it was created on (HIPFFT_INVALID_DEVICE)"
    case(12); msg = "internal plan database error (HIPFFT_PARSE_ERROR)"
    case(13); msg = "no workspace provided prior to plan execution (HIPFFT_NO_WORKSPACE)"
    case(14); msg = "functionality not implemented for the given parameters (HIPFFT_NOT_IMPLEMENTED)"
    case(15); msg = "license error (HIPFFT_LICENSE_ERROR)"
    case(16); msg = "operation not supported for the given parameters (HIPFFT_NOT_SUPPORTED)"
    case default; msg = "unrecognized hipFFT status code"
    end select

    write(stdout, "(A,A,A,I0,A,A,A)") " hipFFT error in ", trim(location), ": status = ", stat, " (", trim(msg), ")"
    call mpi_abend()
    stop "hipFFT call failed"
  end subroutine check_hipfft_status

  ! 3D-FFT, hipFFT version

  subroutine fft_init_rtc(handle, in, out)
    type(fft_handle), intent(out) :: handle
    real(wp), intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex(wp), intent(out) :: out(fftsize(1)/2+1, fftsize(2), fftsize(3))
    integer :: stat
#ifdef DP
    stat = hipfftPlan3d(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         HIPFFT_D2Z)
#else
    stat = hipfftPlan3d(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         HIPFFT_R2C)
#endif
    call check_hipfft_status(stat, "hipfftPlan3d (fft_init_rtc)")
  end subroutine fft_init_rtc

  subroutine fft_init_ctr(handle, in, out)
    type(fft_handle), intent(out) :: handle
    complex(wp), intent(in) :: in(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real(wp), intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
    integer :: stat
#ifdef DP
    stat = hipfftPlan3d(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         HIPFFT_Z2D)
#else
    stat = hipfftPlan3d(handle%plan, fftsize(1), fftsize(2), fftsize(3), &
         HIPFFT_C2R)
#endif
    call check_hipfft_status(stat, "hipfftPlan3d (fft_init_ctr)")
  end subroutine fft_init_ctr

  subroutine fft_rtc(handle, in, out)
    use hipfort_hipfft
    type(fft_handle), intent(in) :: handle
    real(wp), intent(in), target :: in(fftsize(1), fftsize(2), fftsize(3))
    complex(wp), intent(out), target :: out(fftsize(1)/2+1, fftsize(2), fftsize(3))
    integer :: stat
    !$omp target data map(alloc: in, out) use_device_addr(in, out)
#ifdef DP
    stat = hipfftExecD2Z(handle%plan, in, out)
#else
    stat = hipfftExecR2C(handle%plan, in, out)
#endif
    !$omp end target data
    call check_hipfft_status(stat, "hipfftExecR2C/D2Z (fft_rtc)")
  end subroutine fft_rtc

  subroutine fft_ctr(handle, in, out)
    use hipfort_hipfft
    type(fft_handle), intent(in) :: handle
    complex(wp), intent(in), target :: in(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real(wp), intent(out), target :: out(fftsize(1), fftsize(2), fftsize(3))
    integer :: stat
    !$omp target data map(alloc: in, out) use_device_addr(in, out)
#ifdef DP
    stat = hipfftExecZ2D(handle%plan, in, out)
#else
    stat = hipfftExecC2R(handle%plan, in, out)
#endif
    !$omp end target data
    call check_hipfft_status(stat, "hipfftExecC2R/Z2D (fft_ctr)")
  end subroutine fft_ctr

  ! ---- batched transforms (used by recpcal.F90's recpcal_prepare_solute_refs) ----
  !
  ! One hipfftPlanMany plan transforms nbatch grids per execution. The
  ! batch stride is one whole grid: Fortran stores a(:,:,:,cnt) with cnt
  ! slowest, so each trial's grid is already contiguous.
  !
  ! Dimension order: hipFFT takes n() slowest-axis first and halves the
  ! LAST axis for R2C. Our arrays halve the FIRST Fortran axis (the
  ! fastest in memory), so n() = (/ fftsize(3), fftsize(2), fftsize(1) /).
  ! (The unbatched hipfftPlan3d calls above pass the sizes the other
  ! way round, which only coincides with this for a cubic grid.)
  !
  ! Two pitfalls found the hard way, both handled below:
  !  * the rank-4 buffers must be passed as c_ptr (hipfort's generics
  !    stop at rank 3), and c_loc() inside a use_device_addr region
  !    returns the HOST address with this compiler -- on an APU that
  !    silently transforms the host copy -- so the device address is
  !    fetched explicitly with omp_get_mapped_ptr;
  !  * hipFFT runs on HIP's stream, unordered with OpenMP's, so each
  !    execution is followed by a device-wide synchronize.

  subroutine fft_init_rtc_batch(handle, nbatch)
    type(fft_handle), intent(out) :: handle
    integer, intent(in) :: nbatch
    integer :: stat, n(3), inembed(3), onembed(3), idist, odist

    n(1) = fftsize(3) ; n(2) = fftsize(2) ; n(3) = fftsize(1)
    inembed(:) = n(:)
    onembed(1) = fftsize(3) ; onembed(2) = fftsize(2)
    onembed(3) = fftsize(1)/2 + 1
    idist = fftsize(1) * fftsize(2) * fftsize(3)
    odist = (fftsize(1)/2 + 1) * fftsize(2) * fftsize(3)
#ifdef DP
    stat = hipfftPlanMany(handle%plan, 3, n, inembed, 1, idist, &
         onembed, 1, odist, HIPFFT_D2Z, nbatch)
#else
    stat = hipfftPlanMany(handle%plan, 3, n, inembed, 1, idist, &
         onembed, 1, odist, HIPFFT_R2C, nbatch)
#endif
    call check_hipfft_status(stat, "hipfftPlanMany (fft_init_rtc_batch)")
    handle%nbatch = nbatch
  end subroutine fft_init_rtc_batch

  subroutine fft_init_ctr_batch(handle, nbatch)
    type(fft_handle), intent(out) :: handle
    integer, intent(in) :: nbatch
    integer :: stat, n(3), inembed(3), onembed(3), idist, odist

    n(1) = fftsize(3) ; n(2) = fftsize(2) ; n(3) = fftsize(1)
    inembed(1) = fftsize(3) ; inembed(2) = fftsize(2)
    inembed(3) = fftsize(1)/2 + 1
    onembed(:) = n(:)
    idist = (fftsize(1)/2 + 1) * fftsize(2) * fftsize(3)
    odist = fftsize(1) * fftsize(2) * fftsize(3)
#ifdef DP
    stat = hipfftPlanMany(handle%plan, 3, n, inembed, 1, idist, &
         onembed, 1, odist, HIPFFT_Z2D, nbatch)
#else
    stat = hipfftPlanMany(handle%plan, 3, n, inembed, 1, idist, &
         onembed, 1, odist, HIPFFT_C2R, nbatch)
#endif
    call check_hipfft_status(stat, "hipfftPlanMany (fft_init_ctr_batch)")
    handle%nbatch = nbatch
  end subroutine fft_init_ctr_batch

  subroutine fft_rtc_batch(handle, in, out)
    use hipfort_hipfft
    use iso_c_binding, only: c_ptr, c_loc, c_associated
    use omp_lib, only: omp_get_mapped_ptr, omp_get_default_device
    type(fft_handle), intent(in) :: handle
    real(wp), intent(in), target :: &
         in(fftsize(1), fftsize(2), fftsize(3), handle%nbatch)
    complex(wp), intent(out), target :: &
         out(fftsize(1)/2+1, fftsize(2), fftsize(3), handle%nbatch)
    integer :: stat, dev
    type(c_ptr) :: p_in, p_out

    dev = omp_get_default_device()
    p_in  = omp_get_mapped_ptr(c_loc(in(1,1,1,1)),  dev)
    p_out = omp_get_mapped_ptr(c_loc(out(1,1,1,1)), dev)
    if (.not. c_associated(p_in) .or. .not. c_associated(p_out)) then
       call check_hipfft_status(4, "fft_rtc_batch: buffer not mapped on device")
       return
    end if
#ifdef DP
    stat = hipfftExecD2Z(handle%plan, p_in, p_out)
#else
    stat = hipfftExecR2C(handle%plan, p_in, p_out)
#endif
    call check_hipfft_status(stat, "hipfftExecR2C/D2Z (fft_rtc_batch)")
    call sync_after_fft("fft_rtc_batch")
  end subroutine fft_rtc_batch

  subroutine fft_ctr_batch(handle, in, out)
    use hipfort_hipfft
    use iso_c_binding, only: c_ptr, c_loc, c_associated
    use omp_lib, only: omp_get_mapped_ptr, omp_get_default_device
    type(fft_handle), intent(in) :: handle
    complex(wp), intent(in), target :: &
         in(fftsize(1)/2+1, fftsize(2), fftsize(3), handle%nbatch)
    real(wp), intent(out), target :: &
         out(fftsize(1), fftsize(2), fftsize(3), handle%nbatch)
    integer :: stat, dev
    type(c_ptr) :: p_in, p_out

    dev = omp_get_default_device()
    p_in  = omp_get_mapped_ptr(c_loc(in(1,1,1,1)),  dev)
    p_out = omp_get_mapped_ptr(c_loc(out(1,1,1,1)), dev)
    if (.not. c_associated(p_in) .or. .not. c_associated(p_out)) then
       call check_hipfft_status(4, "fft_ctr_batch: buffer not mapped on device")
       return
    end if
#ifdef DP
    stat = hipfftExecZ2D(handle%plan, p_in, p_out)
#else
    stat = hipfftExecC2R(handle%plan, p_in, p_out)
#endif
    call check_hipfft_status(stat, "hipfftExecC2R/Z2D (fft_ctr_batch)")
    call sync_after_fft("fft_ctr_batch")
  end subroutine fft_ctr_batch

  subroutine sync_after_fft(location)
    use hipfort, only: hipDeviceSynchronize, hipSuccess
    use engmain, only: stdout
    use mpiproc, only: mpi_abend
    implicit none
    character(len=*), intent(in) :: location
    integer :: stat
    stat = hipDeviceSynchronize()
    if (stat /= hipSuccess) then
       write(stdout, "(A,A,A,I0)") " hipDeviceSynchronize failed after ", &
            trim(location), ": status = ", stat
       call mpi_abend()
       stop "hipDeviceSynchronize failed"
    end if
  end subroutine sync_after_fft

  subroutine fft_cleanup_rtc(handle)
    type(fft_handle), intent(in) :: handle
    integer :: stat
    stat = hipfftDestroy(handle%plan)
    call check_hipfft_status(stat, "hipfftDestroy (fft_cleanup_rtc)")
  end subroutine fft_cleanup_rtc

  subroutine fft_cleanup_ctr(handle)
    type(fft_handle), intent(in) :: handle
    integer :: stat
    stat = hipfftDestroy(handle%plan)
    call check_hipfft_status(stat, "hipfftDestroy (fft_cleanup_ctr)")
  end subroutine fft_cleanup_ctr

end module fft_iface
