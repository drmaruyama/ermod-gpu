#                                               -*- Autoconf -*-

# autoconf <= 2.64 has a bug in AC_OPENMP
AC_PREREQ([2.65])
AC_INIT(erdst-openacc, 1.0.0, Matubayasi Laboratory)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE

AC_ARG_ENABLE([debug-symbol],
              [AS_HELP_STRING([--enable-debug-symbol],
                        [Enable debug symbols (-g)])],
              [enable_debug=yes], [enable_debug=no])

AS_IF([test x$enable_debug = xno],
      [ac_cv_prog_cc_g=no
       ac_cv_prog_fc_g=no])

AC_ARG_WITH([gpu],
    [AS_HELP_STRING([--with-gpu=nvidia|amd|none],
        [Select GPU offload toolchain. Default: auto-detect from the
         Fortran compiler (nvfortran -> nvidia, amdflang/flang-new -> amd)])],
    [with_gpu=$withval], [with_gpu=auto])

AC_ARG_WITH([offload-arch],
    [AS_HELP_STRING([--with-offload-arch=ARCH],
        [AMD GPU target architecture passed to --offload-arch, e.g. gfx90a
         (MI200), gfx942 (MI300). Only used when --with-gpu=amd. Default:
         "native" (auto-detect the GPU on the build machine); use
         "rocminfo | grep gfx" to find the code by hand if "native" is not
         accepted by your compiler version.])],
    [with_offload_arch=$withval], [with_offload_arch=native])

AC_ARG_WITH([fc-optlevel],
    [AS_HELP_STRING([--with-fc-optlevel=LEVEL],
        [Optimization level (0/1/2/3) passed to the Fortran compiler.
         Default: 3 for NVIDIA, 2 for AMD. The AMD default avoids a known
         amdflang issue where device kernel symbols can be silently
         optimized away at -O3 (see llvm/llvm-project issue #61483),
         causing HSA_STATUS_ERROR_INVALID_SYMBOL_NAME / "requested object
         was not found in the binary image" at run time; -O2 was not
         observed to trigger it. Once your ROCm/amdflang version has the
         -O3 issue fixed, raise this and re-verify correctness at
         production scale before trusting any performance gain.])],
    [with_fc_optlevel=$withval], [with_fc_optlevel=auto])

# Checks for programs.
AC_PROG_CPP
# nvc (NVIDIA HPC SDK's C compiler) is only relevant/available on NVIDIA
# systems, and on other systems (or if a stray nvc happens to be left on
# PATH, e.g. from an unrelated module) it can be broken or reference
# paths that don't exist there. Skip it whenever the environment already
# points at an AMD toolchain, either via --with-gpu=amd or via FC already
# naming an AMD compiler (as when a "module load rocm/..." sets FC).
cc_candidates="nvc gcc cc"
AS_IF([test "x$with_gpu" = "xamd"], [cc_candidates="gcc clang cc"])
AS_CASE([$FC],
        [*amdflang*|*flang-new*|*flang*], [cc_candidates="gcc clang cc"])
AC_PROG_CC([$cc_candidates])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
# nvfortran: NVIDIA HPC SDK. amdflang / flang-new: AMD's ROCm/LLVM-based
# Fortran compiler (loaded via e.g. "module load rocm/7.2.0" or one of the
# "..._amdflang" combo modules). If FC is already set in the environment
# (as when a module sets it, or the person exports FC=amdflang by hand),
# that value is used as-is and this candidate list is ignored.
AC_PROG_FC([nvfortran amdflang flang-new])
AC_PROG_LN_S

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran)
ACX_MPI_MY
AC_LANG(Fortran)
ACX_MPI_MY([use_mpi=yes], [use_mpi=no])

AS_IF([test "x$FC" = "x"],
      [AC_MSG_ERROR("Failed to find Fortran Compiler (Fortran 90). Set environment variable FC to work around")]
      [])

AC_ARG_ENABLE([vmdplugin],
	AS_HELP_STRING([--disable-vmdplugin],
		[disable VMD plugin]))

################################

AC_ARG_ENABLE([double],
    [AS_HELP_STRING([--enable-double],
        [Enable conversion of single precision real to double precision (real(8))])],
    [enable_double=yes],
    [enable_double=no])

# Precision (single vs. double) is selected entirely by the -DDP
# preprocessor definition, consumed by precision.F90's "wp" kind
# parameter. Unlike before, no compiler flag that redefines what a
# bare "real" means (e.g. nvfortran's -Mr8) is needed anymore, since
# every real/complex declaration in the source spells out "real(wp)"
# explicitly.
AS_IF([test "x$enable_double" = "xyes"],
      [FCFLAGS="$FCFLAGS -DDP"])

################################

AC_LANG_PUSH(Fortran)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [dnl
      implicit none
      integer(8) :: a = 1
      real(8) :: b = 1.0d0
      
      a = transfer(b, a)
      print *, a
])], [my_check_have_transfer=yes], [my_check_have_transfer=no])

echo -n $my_check_have_transfer
AC_LANG_POP(Fortran)

AM_CONDITIONAL([USE_TRANSFER], [test "x$my_check_have_transfer" = "xyes"])

################################

AC_LANG_PUSH(Fortran)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [dnl
      implicit none

      open(10,file='/tmp/configure.ermod.flush.test')
      write(10,*) 342
      flush(10)
      close(10)
])], [my_check_have_flush=yes], [my_check_have_flush=no])

echo -n $my_check_have_flush
AC_LANG_POP(Fortran)

AM_CONDITIONAL([USE_FLUSH], [test "x$my_check_have_flush" = "xyes"])

################################

# set to f90 as default
AC_FC_SRCEXT(f90)

AC_ARG_ENABLE([mpi],
	AS_HELP_STRING([--disable-mpi],
		[disable MPI even if exists]),
	AS_IF([test x$enable_mpi = xno],
		     [use_mpi=no]),
	[])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# FIXME: check cpp flags

AC_ARG_ENABLE([defaults],
	AS_HELP_STRING([--disable-defaults],
		[(deprcated)]))

AS_IF([test x$enable_defaults != xno],
CFLAGS="$CFLAGS"
FCFLAGS="$FCFLAGS")

AM_CONDITIONAL(USE_MPI, [test "x$use_mpi" = "xyes"])

AS_IF([test "x$use_mpi" = "xyes"],
	    [mpi_result="enabled"
	     FC_CANONICAL="$FC"
	     FC="$MPIFC"
             FCFLAGS="-DMPI $FCFLAGS"],
	    [mpi_result="disabled"
	     FC_CANONICAL="$FC"])
AC_SUBST([FC_CANONICAL])

# optimization / GPU offload flags
#
# with_gpu=auto detects the vendor from the (pre-MPI-wrapping) Fortran
# compiler name recorded in FC_CANONICAL just above. Override with
# --with-gpu=nvidia / --with-gpu=amd / --with-gpu=none if the name-based
# guess is wrong (e.g. a custom-named wrapper script).
AS_IF([test "x$with_gpu" = "xauto"],
      [AS_CASE([$FC_CANONICAL],
               [*nvfortran*], [with_gpu=nvidia],
               [*amdflang*|*flang-new*|*flang*], [with_gpu=amd],
               [with_gpu=none])])

AC_ARG_WITH([hipfort],
    [AS_HELP_STRING([--with-hipfort=DIR],
        [Root directory containing hipfort's AMD ("amdgcn") build, i.e.
         the directory such that DIR/include/hipfort/amdgcn/hipfort_hipfft.mod
         and DIR/lib/libhipfort-amdgcn.a both exist. Default: $ROCM_PATH,
         since ROCm 7.x ships hipfort bundled inside the ROCm tree itself.])],
    [with_hipfort=$withval], [with_hipfort=auto])

GPU_LIBS=
AS_CASE([$with_gpu],
  [nvidia], [dnl
    AS_IF([test "x$with_fc_optlevel" = "xauto"], [with_fc_optlevel=3])
    # -acc is kept only for backward compatibility while any file still
    # has "!$acc" directives; as of this branch every file has been
    # converted to OpenMP target offload, so -acc is a no-op and can be
    # dropped once you are sure no OpenACC directives remain anywhere.
    FCFLAGS="$FCFLAGS -fast -Minfo=accel -mp=gpu -cudalib=cufft"
  ],
  [amd], [dnl
    AS_IF([test "x$with_fc_optlevel" = "xauto"], [with_fc_optlevel=2])
    # ROCm's ROCM_PATH is normally exported by "module load rocm/x.y.z";
    # fall back to the conventional /opt/rocm install location otherwise.
    AS_IF([test "x$ROCM_PATH" = "x"], [ROCM_PATH=/opt/rocm])
    FCFLAGS="$FCFLAGS -O$with_fc_optlevel -cpp -fopenmp --offload-arch=$with_offload_arch"

    # As of ROCm 7.x, hipfort ships bundled inside the ROCm tree itself
    # (not as a separate package), with prebuilt .mod files for both the
    # "amdgcn" (AMD GPU) and "nvptx" (NVIDIA GPU, via hipfort's CUDA
    # backend) targets side by side:
    #   $ROCM_PATH/include/hipfort/amdgcn/hipfort_hipfft.mod
    #   $ROCM_PATH/lib/libhipfort-amdgcn.a
    # --with-hipfort=DIR overrides the root DIR to search under (still
    # expecting the same include/hipfort/amdgcn and lib layout below it)
    # for installations that keep hipfort somewhere else.
    AS_IF([test "x$with_hipfort" = "xauto"], [with_hipfort=$ROCM_PATH])

    AS_IF([test ! -f "$with_hipfort/include/hipfort/amdgcn/hipfort_hipfft.mod"],
      [AC_MSG_ERROR([dnl
Could not find $with_hipfort/include/hipfort/amdgcn/hipfort_hipfft.mod.
Locate hipfort by hand, e.g.:
  find $ROCM_PATH -iname "hipfort_hipfft.mod" 2>/dev/null
  find / -iname "hipfort_hipfft.mod" 2>/dev/null
then re-run configure with --with-hipfort=DIR, where DIR is the
directory that contains both "include/hipfort/amdgcn" and "lib".])])

    FCFLAGS="$FCFLAGS -I$with_hipfort/include/hipfort/amdgcn"
    GPU_LIBS="-L$with_hipfort/lib -lhipfort-amdgcn -lhipfft -lamdhip64"

    # History note: builds for a while forced the final link step (where
    # amdflang re-optimizes every file's OpenMP device code together)
    # down to -O1, because linking at -O2 gave rare, grossly wrong
    # energies at production scale. That turned out to be triggered by
    # per-thread private *arrays* in realcal.F90's kernels (xst(3)),
    # which the device compiler places in scratch memory; with those
    # rewritten as scalars, the whole program links correctly at the
    # same level it compiles at, verified at production scale with
    # several random seeds. Linking at -O1 had been costing far more
    # than any single kernel: refs went from 14m11s to 4m26s and soln
    # from 34.0s to 19.7s once it was removed. If you add GPU kernels,
    # avoid private arrays in them for the same reason.
  ],
  [none], [dnl
    AC_MSG_WARN([Could not auto-detect a GPU vendor from FC="$FC_CANONICAL". Building without any GPU offload flags. Pass --with-gpu=nvidia or --with-gpu=amd to force a toolchain.])
  ],
  [AC_MSG_ERROR([Unknown value for --with-gpu: $with_gpu (expected nvidia, amd, or none)])]
)
AC_SUBST([GPU_LIBS])
AC_DEFUN([MY_FC_OPT_ACCEPT_NATIVE],[
AC_CACHE_CHECK([fortran compiler accepts -march=native],
my_fc_opt_accept_native_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
flag="-march=native"
FCFLAGS="$save_FCFLAGS $flag"
AC_COMPILE_IFELSE([      program main
      end program], [my_fc_opt_accept_native=yes], [my_fc_opt_accept_native=no])
echo -n $my_fc_opt_accept_native
FCFLAGS=$save_FCFLAGS
AC_LANG_POP(Fortran)
])])

# AM_CONDITIONAL([MPI], [test x$use_mpi = xyes])
AM_CONDITIONAL([USE_VMDPLUGIN], [test x$enable_vmdplugin = xyes || test x$enable_vmdplugin = x])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AC_MSG_RESULT([dnl
MPI: $mpi_result
GPU vendor: $with_gpu
FC(for ermod)=$FC
FCFLAGS=$FCFLAGS
GPU_LIBS=$GPU_LIBS
LDFLAGS=$LDFLAGS
LIBS=$LIBS

Configuration complete; type \"make\" to start compiling])
