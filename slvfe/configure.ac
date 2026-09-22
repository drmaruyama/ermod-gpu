#                                               -*- Autoconf -*-

# autoconf <= 2.64 has a bug in AC_OPENMP
AC_PREREQ([2.65])
AC_INIT(slvfe, 1.0.0, Matubayasi Laboratory)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR([slvfe.F90])
AM_INIT_AUTOMAKE

AC_ARG_ENABLE([debug-symbol],
              [AS_HELP_STRING([--enable-debug-symbol],
                        [Enable debug symbols (-g)])],
              [enable_debug=yes], [enable_debug=no])

AS_IF([test x$enable_debug = xno],
       ac_cv_prog_fc_g=no])

AC_ARG_WITH([gpu],
    [AS_HELP_STRING([--with-gpu=nvidia|amd],
        [GPU toolchain for the linear-algebra solver. Default: detected
         from the Fortran compiler (nvfortran -> nvidia, amdflang or
         flang-new -> amd). nvidia uses OpenACC + cuSOLVER, amd uses
         OpenMP target offload + hipSOLVER.])],
    [with_gpu=$withval], [with_gpu=auto])

AC_ARG_WITH([offload-arch],
    [AS_HELP_STRING([--with-offload-arch=ARCH],
        [AMD GPU architecture passed to --offload-arch, e.g. gfx942
         (MI300A/MI300X), gfx90a (MI210/MI250). Only used for
         --with-gpu=amd. Default: native, which needs a GPU on the build
         machine; set it explicitly when building on a login node.])],
    [with_offload_arch=$withval], [with_offload_arch=native])

# Checks for programs.
AC_PROG_FC([nvfortran amdflang flang-new])
AC_PROG_LN_S

# Check for buggy variant of ifort.
AS_IF([test x$FC = xifort],
      [version=`$FC --version`
       case "x$version" in # (
       *12.1.*):
         AC_MSG_ERROR([This program does not run correctly with Intel Compiler 12.1.x])
       ;;
       esac
       ])

# modified MPI detector (searches mpiifort first)
AC_LANG(Fortran)

: ${enable_openmp=no}
AC_ARG_ENABLE([openmp],
              AS_HELP_STRING([--enable-openmp],
                             [use OpenMP (default = no)]))
AC_OPENMP

AS_IF([test "x$FC" = "x"],
      [AC_MSG_ERROR("Failed to find Fortran Compiler (Fortran 90). Set environment variable FC to work around")]
      [])

################################

AC_DEFUN([MY_FC_REAL_DOUBLE],[
AC_CACHE_CHECK([how to read real as real(8) as default],
my_fc_real_double_cv_flag,
[AC_LANG_PUSH(Fortran)
save_FCFLAGS=$FCFLAGS
# for default, gfortran, ifort, xlf, Fujitsu f90, Fujitsu f90 (old), Portland, Cray
for flag in "" "-fdefault-real-8 -fdefault-double-8" "-fdefault-real-8" "-autodouble" "-qrealsize=8" "-CcdRR8" "-Ad" "-Mr8" "-s default64"; do
FCFLAGS="$save_FCFLAGS $flag"
# This is very tricky indeed...
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [dnl
      implicit none
      real :: a
      integer, parameter :: len = kind(a)
      integer, parameter :: x(len) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
      print *, "Success"
])], [my_check_realasdouble=yes], [my_check_realasdouble=no])
  if test $my_check_realasdouble = yes; then
    break
  fi
done
if test $my_check_realasdouble = no; then
  AC_MSG_ERROR([unable to find compiler flag to read real as real(8)])
fi
FCFLAGS="$save_FCFLAGS $flag"
if test "x$flag" = x; then
  echo -n "(none)"
else
  echo -n $flag
fi
AC_LANG_POP(Fortran)
])])

MY_FC_REAL_DOUBLE
FCFLAGS="$FCFLAGS"

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

# GPU toolchain
AS_IF([test "x$with_gpu" = "xauto"],
      [AS_CASE([$FC],
               [*nvfortran*], [with_gpu=nvidia],
               [*amdflang*|*flang-new*|*flang*], [with_gpu=amd],
               [AC_MSG_ERROR([dnl
Could not tell which GPU toolchain to use from FC="$FC". Pass
--with-gpu=nvidia (nvfortran) or --with-gpu=amd (amdflang).])])])

GPU_LIBS=
AS_CASE([$with_gpu],
  [nvidia], [dnl
    FCFLAGS="$FCFLAGS -fast -Minfo=accel -acc -cudalib=cublas,cusolver"
  ],
  [amd], [dnl
    AS_IF([test "x$ROCM_PATH" = "x"], [ROCM_PATH=/opt/rocm])
    # hipSOLVER is called through slvfe.F90's own bind(C) interfaces;
    # only the enumerator values come from hipfort, hence the include
    # path and no hipfort library.
    AS_IF([test ! -f "$ROCM_PATH/include/hipfort/amdgcn/hipfort_hipsolver_enums.mod"],
      [AC_MSG_ERROR([dnl
Could not find $ROCM_PATH/include/hipfort/amdgcn/hipfort_hipsolver_enums.mod.
Load your ROCm environment (which normally sets ROCM_PATH) and re-run.])])
    # -O2 rather than -O3: amdflang can drop GPU kernels from the
    # executable at -O3 (llvm/llvm-project issue #61483).
    FCFLAGS="$FCFLAGS -O2 -fopenmp --offload-arch=$with_offload_arch -DUSE_HIP"
    FCFLAGS="$FCFLAGS -I$ROCM_PATH/include/hipfort/amdgcn"
    GPU_LIBS="-L$ROCM_PATH/lib -lhipsolver -lamdhip64"
  ],
  [AC_MSG_ERROR([Unknown value for --with-gpu: $with_gpu (expected nvidia or amd)])]
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

AC_ARG_ENABLE([opt],
	AS_HELP_STRING([--disable-opt],
		[disable automatic edit to optimization flags]))

# The AMD branch above sets its own optimisation level; the generic
# editor below would append a bare -O for a non-GNU, non-Intel compiler,
# which would override that -O2, so it is skipped there.
AS_IF([test "x$with_gpu" != "xamd" && (test x$enable_opt = xyes || test x$enable_opt = x)],
  [FCFLAGS="$FCFLAGS"
   AS_IF([test x$ac_compiler_gnu = xyes],
        [optflag="-O3 -ffast-math -fno-finite-math-only"
         MY_FC_OPT_ACCEPT_NATIVE
         AS_IF([test x$my_fc_opt_accept_native = xyes],
                optflag="$optflag -march=native")],
        [case x$FC in #(
           xmpiifort|xifort)
             optflag="-O3 -no-prec-div -xHost"
             ;; #(
           *)
             optflag="-O"
             ;;
         esac])
   FCFLAGS="$FCFLAGS $optflag"])

# Check OpenMP
AS_IF([test "$enable_openmp" != "no"],
  [FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
   # Check OpenMP version
   AC_MSG_CHECKING([whether OpenMP version >= 3.0])
   AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM([], [dnl
!\$omp task
  print *, "test"
!\$omp end task
])], 
     [AC_MSG_RESULT([yes])], 
     [AC_MSG_RESULT([no])
      AC_MSG_ERROR([OpenMP must be >= 3.0])])])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AC_MSG_RESULT([dnl
OpenMP: $enable_openmp
GPU: $with_gpu
FC(for ermod)=$FC
FCFLAGS=$FCFLAGS
GPU_LIBS=$GPU_LIBS

Configuration complete; type \"make\" to start compiling])

