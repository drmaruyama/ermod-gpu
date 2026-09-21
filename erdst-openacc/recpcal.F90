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
! [OpenMP target offload 検証版]
! realcal.F90と同じ方針でOpenACCからOpenMP target offloadへ変換した
! 検証用ブランチのバージョンです。数値計算ロジックは変更していません。
!
! realcal.F90との違い:
! この module reciprocal が所有するデータ(slvtag, splslv, grdslv,
! cnvslt, engfac, rcpslt)については、確保(enter data)・転送(update)を
! 含めて全てOpenMP targetに統一しています。したがってこれらのデータに
! 関してはOpenACC/OpenMPの presence table 共有(realcal.F90のコメント
! 参照)には依存していません。
! 一方、engmain(mol_begin_index, charge, numsite, sluvid)や
! engproc.F90/呼び出し元(uvengy, tagpt, sltlist)側のデータは
! 引き続きsetconf.F90/engproc.F90でOpenACC管理のままなので、そちらは
! 従来どおりpresence table共有に依存した map(alloc: ...) を使っています。
!
! 変換にあたっての注意点:
! 1. map(alloc: ...) は当初 map(present, alloc: ...) を使う想定でしたが、
!    realcal.F90の変換時にnvfortranがOpenMP 5.1のpresent修飾子付き構文で
!    コンパイルエラーになることが分かったため、ここでも "alloc" のみに
!    しています。対象データが本当にデバイス上に存在するかは実行時の
!    暗黙のチェックに委ねられる点は realcal.F90 と同じです。
! 2. "!$acc loop seq"(cg1ループ)には対応するOpenMP指示文がありません。
!    指示文を付けないループは、外側の "target teams distribute parallel do"
!    が割り当てたスレッド内で逐次実行されるため、単純に指示文を削除して
!    います(realcal_energy_soln/refs内の該当箇所にコメントを残しています)。
! 3. reduction(+:self_energy) を伴う "target teams distribute parallel do"
!    は、reduction句がACCと同様にサポートされている前提です。realcal.F90
!    では reduction を使うカーネルが無かったため、ここが今回新たに検証が
!    必要な箇所です。数値一致の確認を特に注意深く行ってください。
! 4. [-Minfo=mpログで判明した点] splval/grdval(ホスト側で
!    calc_spline_molecule(_refs)が書き込む局所配列)をmap句に明示していな
!    かったため、nvfortranが "implicit map(tofrom: splval, grdval)" を
!    生成し、カーネル呼び出しのたびに読み取り専用のはずのこれらの配列を
!    コピーイン+コピーアウトの両方向で転送していた。OpenACC版の対応箇所は
!    copyinのみだったため、この分だけOpenMP版が余分にコストを払っていた
!    とみられる。読み取り専用であることを明示する "map(to: splval, grdval)"
!    を追加して修正済み。
!    ただし -Minfo=mp は(-Minfo=accelと異なり)teams/threadsへの割り当て
!    ("grid"/"block"に相当する情報)を出力しないため、recpcal_prepare_solute
!    /_refs で観測された全体の実行時間差(約3倍)がこの転送コストだけで
!    説明しきれるかは未確認。実際のカーネル構成の比較にはnsys profileでの
!    実測が別途必要。
! =====================================================================

module reciprocal
  use precision_kinds, only: wp
  use fft_iface, only: fft_handle
  implicit none
  integer :: rc1min, rc1max, rc2min, rc2max, rc3min, rc3max
  integer :: ccesize, ccemax
  integer, allocatable :: slvtag(:)
  real(wp),    allocatable :: engfac(:,:,:)
  real(wp),    allocatable :: gf_b(:)
  complex(wp), allocatable :: rcpslt(:,:,:)
  real(wp),    allocatable :: splslv(:,:,:)
  integer, allocatable :: grdslv(:,:)
  real(wp),    allocatable :: cnvslt(:,:,:,:)
  real(wp),    allocatable :: splfc1(:), splfc2(:), splfc3(:)
  complex(wp), allocatable :: fft_buf(:, :, :)

  ! One reciprocal-space grid, and one self-energy value, per "slot":
  ! either per solute molecule (SLT_SOLN) or per insertion trial
  ! (SLT_REFS_*). slttype is fixed for the whole run, so exactly one of
  ! these two uses is ever active -- there is no need for separate
  ! cnvslt/cnvslt_r or solute_self_energy/_refs arrays.
  real(wp),    allocatable :: solute_self_energy(:)

  type(fft_handle) :: handle_c2r, handle_r2c

contains
  subroutine recpcal_init(slvmax, tagpt)
    use engmain, only:  nummol, numsite, splodr, ms1max, ms2max, ms3max, &
         maxins, slttype, SLT_SOLN, numslt
    use spline, only: spline_init
    use fft_iface, only: fft_init_ctr, fft_init_rtc, fft_set_size
    implicit none
    integer, intent(in) :: slvmax, tagpt(:)
    integer :: m, k
    integer :: gridsize(3), ptrnk

    allocate( slvtag(nummol) )
    !$omp target enter data map(alloc: slvtag)
    slvtag(:) = -1
    ptrnk = 0
    do k = 1, slvmax
       m = tagpt(k)
       slvtag(m) = ptrnk + 1
       ptrnk = ptrnk + numsite(m)
    enddo
    !$omp target update to(slvtag)

    allocate(gf_b(splodr)) ! for PPPM Green's function
    call calc_gfb_pppm()   ! calc gf_b (function of splodr)
    rc1min = 0 ; rc1max = ms1max - 1
    rc2min = 0 ; rc2max = ms2max - 1
    rc3min = 0 ; rc3max = ms3max - 1
    ccesize = ms1max / 2 + 1; ccemax = ccesize - 1
    call spline_init(splodr)
    allocate(splslv(0:splodr-1, 3, ptrnk), grdslv(3, ptrnk))
    !$omp target enter data map(alloc: splslv, grdslv)
    ! One reciprocal-space grid per "slot" (solute molecule for
    ! SLT_SOLN, insertion trial for SLT_REFS_*), so each can be
    ! spread/FFT'd/evaluated without clobbering the others (see
    ! recpcal_prepare_solute).
    if (slttype == SLT_SOLN) then
       allocate(solute_self_energy(numslt))
       allocate(cnvslt(rc1min:rc1max, rc2min:rc2max, rc3min:rc3max, numslt))
    else
       allocate(solute_self_energy(maxins))
       allocate(cnvslt(rc1min:rc1max, rc2min:rc2max, rc3min:rc3max, maxins))
    end if
    !$omp target enter data map(alloc: cnvslt)
    ! initialize spline table for all axes
    allocate( splfc1(rc1min: rc1max) )
    allocate( splfc2(rc2min: rc2max) )
    allocate( splfc3(rc3min: rc3max) )
    call init_spline_axis(rc1min, rc1max, splfc1(rc1min:rc1max))
    call init_spline_axis(rc2min, rc2max, splfc2(rc2min:rc2max))
    call init_spline_axis(rc3min, rc3max, splfc3(rc3min:rc3max))
    gridsize(1) = ms1max
    gridsize(2) = ms2max
    gridsize(3) = ms3max
    call fft_set_size(gridsize)
    allocate( engfac(rc1min:ccemax, rc2min:rc2max, rc3min:rc3max) )
    allocate( rcpslt(rc1min:ccemax, rc2min:rc2max, rc3min:rc3max) )
    !$omp target enter data map(alloc: engfac, rcpslt)
    ! init fft (each grid slice has the same shape; slice 1 always exists)
    call fft_init_rtc(handle_r2c, cnvslt(:,:,:,1), rcpslt)
    call fft_init_ctr(handle_c2r, rcpslt, cnvslt(:,:,:,1))
  end subroutine recpcal_init

  subroutine init_spline_axis(imin, imax, splfc)
    use engmain, only: splodr, PI
    use spline, only: spline_value
    implicit none
    integer, intent(in) :: imin, imax
    real(wp), intent(out) :: splfc(imin:imax)
    real(wp) :: chr, factor, rtp2
    real(wp) :: cosk, sink
    complex(wp) :: rcpi
    integer :: rci, spi
    do rci = imin, imax
       rcpi = (0.0_wp, 0.0_wp)
       do spi = 0, splodr - 2
          chr = spline_value(real(spi + 1, wp))
          rtp2 = 2.0_wp * PI * real(spi * rci, wp) / real(imax + 1, wp)
          cosk = chr * cos(rtp2)
          sink = chr * sin(rtp2)
          rcpi = rcpi + cmplx(cosk, sink, wp)
       end do
       factor = real(rcpi * conjg(rcpi), wp)
       splfc(rci) = factor
    end do
  end subroutine init_spline_axis

  subroutine recpcal_spline_greenfunc()
    use engmain, only: invcl, ms1max, ms2max, ms3max, splodr, volume, screen, PI
    implicit none
    integer :: rc1, rc2, rc3, rci, m, rcimax
    real(wp) :: factor, rtp2, chr
    real(wp) :: inm(3), xst(3)
    logical :: at_nyquist
    do rc3 = rc3min, rc3max
       do rc2 = rc2min, rc2max
          do rc1 = rc1min, ccemax
             factor = 0.0_wp
             if (rc1 == 0 .and. rc2 == 0 .and. rc3 == 0) cycle
             at_nyquist = .false.
             do m = 1, 3
                if (m == 1) rci = rc1
                if (m == 2) rci = rc2
                if (m == 3) rci = rc3

                if (m == 1) rcimax = ms1max
                if (m == 2) rcimax = ms2max
                if (m == 3) rcimax = ms3max

                if ((mod(splodr, 2) == 1) .and. (2*abs(rci) == rcimax)) then
                   ! at the Nyquist frequency along this axis: leave factor at 0
                   at_nyquist = .true.
                   exit
                end if
                if (rci <= rcimax / 2) then
                   inm(m) = real(rci, wp)
                else
                   inm(m) = real(rci - rcimax, wp)
                endif
             end do
             if (.not. at_nyquist) then
                do m = 1, 3
                   xst(m) = dot_product(invcl(:, m), inm(:))
                end do
                rtp2 = sum(xst(1:3) ** 2)
                chr = (PI ** 2) * rtp2 / (screen ** 2 )
                factor = exp(-chr) / rtp2 / PI / volume
                rtp2 = splfc1(rc1) * splfc2(rc2) * splfc3(rc3)
                factor = factor / rtp2
             end if
             engfac(rc1, rc2, rc3) = factor
          end do
       end do
    end do
    engfac(0, 0, 0) = 0.0
    !$omp target update to(engfac)
  end subroutine recpcal_spline_greenfunc

  ! implemented only for PPPM
  function factorial(n) result(fact)
    implicit none
    integer, intent(in) :: n
    integer(kind=8) :: fact, i
    if (n < 0) stop "Error"
    fact = 1
    if (n > 0) then
       do i = 1, n
          fact = fact * i
       end do
    end if
  end function factorial

  ! calc gf_b needed by denominator in Green's func for PPPM
  subroutine calc_gfb_pppm()
    use engmain, only: splodr
    implicit none
    integer :: l, m
    real(wp) :: gaminv

    gf_b(2:splodr) = 0.
    gf_b(1) = 1.
    do m = 1, splodr-1
       do l = m, 1, -1
          gf_b(l+1)=4.*(gf_b(l+1)*(l-m)*(l-m-0.5)-gf_b(l)*(l-m-1)*(l-m-1))
       end do
       gf_b(1)=4.*(gf_b(1)*(l-m)*(l-m-0.5))
    end do
    gaminv = 1./ factorial(2*splodr-1)
    do m = 1, splodr
       gf_b(m) = gf_b(m) * gaminv
    end do
  end subroutine calc_gfb_pppm

  ! calc denominator in Green's func for PPPM
  function gf_denom(x,y,z) result(denom)
    use engmain, only: splodr
    implicit none
    real(wp), intent(in) :: x, y, z
    real(wp) :: s(3), denom
    integer :: l
    s(1:3) = 0._wp
    do l = splodr, 1, -1
       s(1) = gf_b(l) + s(1) * x
       s(2) = gf_b(l) + s(2) * y
       s(3) = gf_b(l) + s(3) * z
    end do
    denom = product(s)**2
  end function gf_denom
    
  ! calc green's function (optimal influence function) for PPPM
  subroutine recpcal_pppm_greenfunc()
    use engmain, only: &
         invcl, ms1max, ms2max, ms3max, splodr, volume, screen, PI, cell
    implicit none
    integer :: rc1, rc2, rc3, i
    integer :: mx, my, mz, mx_max, my_max, mz_max, mN(3)
    real(wp) :: factor
    real(wp) :: inm(3) ! folded rc1, rc2 or rc3
    real(wp) :: k(3), sin2kh2(3), k2, sum_dru2, tmp, kmk(3), km(3)
    real(wp) :: gam(3), m_max(3), ukm(3)
    real(wp),parameter :: EPS_HOC = 0.0000001_wp ! this value is adopted in LAMMPS 10Feb15
    m_max(1) = screen/(PI*ms1max)*((-log(EPS_HOC))**0.25) ! LAMMPS' form
    m_max(2) = screen/(PI*ms2max)*((-log(EPS_HOC))**0.25)
    m_max(3) = screen/(PI*ms3max)*((-log(EPS_HOC))**0.25)
    mx_max = int(dot_product(cell(:, 1), m_max(:)))
    my_max = int(dot_product(cell(:, 2), m_max(:)))
    mz_max = int(dot_product(cell(:, 3), m_max(:)))

!    write (6,*) "recpcal_pppm_spline_greenfunc"
!    write (6,*) m_max
!    write (6,*) cell(:,1)
!    write (6,*) cell(:,2)
!    write (6,*) cell(:,3)
!    write (6,*) dot_product(cell(:, 1), m_max(:))
!    write (6,*) 'Calculating PPPM Green function...'
!    write (6,*) ' Range for m vector:'
!    write (6,*) '  [-', mx_max,':',mx_max,']'
!    write (6,*) '  [-', my_max,':',my_max,']'
!    write (6,*) '  [-', mz_max,':',mz_max,']'

    if (mx_max < 0 .or. my_max < 0 .or. mz_max < 0) &
         stop "cannot set range of m vectors for Green\'s function correctly"
    if (mod(splodr,2) == 1) &
         stop "interpolation order for PPPM must be even"

    do rc3 = rc3min, rc3max
       if (rc3 <= ms3max/2) then
          inm(3) = rc3
       else
          inm(3) = rc3 - ms3max
       end if
       sin2kh2(3) = (sin(PI*inm(3)/ms3max))**2 !! sin^2(kz*hz/2)

       do rc2 = rc2min, rc2max
          if (rc2 <= ms2max/2) then
             inm(2) = rc2
          else
             inm(2) = rc2 - ms2max
          end if
          sin2kh2(2) = (sin(PI*inm(2)/ms2max))**2 !! sin^2(ky*hy/2)

          do rc1 = rc1min, ccemax
             if (rc1 <= ms1max/2) then
                inm(1) = rc1
             else
                inm(1) = rc1 - ms1max
             end if
             sin2kh2(1) = (sin(PI*inm(1)/ms1max))**2 !! sin^2(kx*hx/2)

             factor = 0.0
             if (rc1 == 0 .and. rc2 == 0 .and. rc3 == 0) cycle
             do i = 1, 3
                k(i) = 2 * PI * dot_product(invcl(:, i), inm(:)) ! k = 2*pi*m
             end do
             k2 = dot_product(k,k) ! k^2

             if (k2 > 0.) then
                sum_dru2 = 0.

                do mx = -mx_max, mx_max
                   tmp = PI*(inm(1)/ms1max+mx)
                   if (tmp == 0.) then
                      ukm(1) = 1.
                   else
                      ukm(1) = (sin(tmp)/tmp)**(2*splodr) ! product(ukm) = U(km)^2
                   end if
                   mN(1) = mx*ms1max

                   do my = -my_max, my_max
                      tmp = PI*(inm(2)/ms2max+my)
                      if (tmp == 0.) then
                         ukm(2) = 1.
                      else
                         ukm(2) = (sin(tmp)/tmp)**(2*splodr)
                      end if
                      mN(2) = my*ms2max

                      do mz = -mz_max, mz_max
                         tmp = PI*(inm(3)/ms3max+mz)
                         if (tmp == 0.) then
                            ukm(3) = 1.
                         else
                            ukm(3) = (sin(tmp)/tmp)**(2*splodr)
                         end if
                         mN(3) = mz*ms2max

                         do i = 1, 3
                            kmk(i) = 2 * PI * dot_product(invcl(:, i), mN(:))
                            ! km - k
                         end do
                         km(1:3) = k(1:3) + kmk(1:3) ! km
                         gam(1:3) = exp(-0.25*(km(1:3)/screen)**2)
                         ! product(gam) = gamma(km)

                         sum_dru2 = sum_dru2 &
                              + dot_product(k,km)*product(gam)*product(ukm) &
                              /dot_product(km,km)
                      end do
                   end do
                end do
                factor = sum_dru2*4.*PI &
                     / (k2*gf_denom(sin2kh2(1),sin2kh2(2),sin2kh2(3))) &
                     / volume
             end if
             engfac(rc1, rc2, rc3) = factor
          end do
       end do
    end do
    engfac(0, 0, 0) = 0.0
    !$omp target update to(engfac)
  end subroutine recpcal_pppm_greenfunc

  ! note: this routine is named as "solvent", but may include solute molecule, when mutiple solute is used.
  subroutine recpcal_prepare_solvent(tagpt, slvmax)
    use engmain, only: numsite
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: tagpt(:), slvmax
    integer :: i, k, svi, stmax

    do k = 1, slvmax
       i = tagpt(k)
       svi = slvtag(i)
       if (svi <= 0) call halt_with_error('rcp_cns')

       stmax = numsite(i)
       call calc_spline_molecule(i, stmax, splslv(:,:,svi:svi+stmax-1), &
            grdslv(:,svi:svi+stmax-1))
    end do
    !$omp target update to(splslv, grdslv)
  end subroutine recpcal_prepare_solvent

  subroutine recpcal_prepare_solute(sltlist, maxdst)
    use engmain, only: ms1max, ms2max, ms3max, sitepos, invcl, numsite, splodr, charge, mol_begin_index
    use fft_iface, only: fft_ctr, fft_rtc
    implicit none
    integer, intent(in) :: sltlist(:), maxdst
    integer :: i, j, k, cnt, tagslt
    integer :: rc1, rc2, rc3, sid, ati, cg1, cg2, cg3, stmax
    real(wp) :: factor, chr
    real(wp) :: self_energy
    logical :: ms1max_even
    real(wp), allocatable, save :: splval(:,:,:)
    integer, allocatable, save :: grdval(:,:)
    logical, save :: initialized = .false.

    ! SLT_SOLN only ever tracks one solute species, so every molecule
    ! in sltlist(1:maxdst) shares the same topology (stmax does not
    ! vary with cnt).
    stmax = numsite(sltlist(1))
    if(.not. initialized) then
       allocate( splval(0:splodr-1, 3, stmax), grdval(3, stmax) )
       initialized = .true.
    end if
    ms1max_even = (mod(ms1max, 2) == 0)
    !$omp target teams distribute parallel do collapse(4) map(alloc: cnvslt)
    do cnt = 1, maxdst
       do k = rc3min, rc3max
          do j = rc2min, rc2max
             do i = rc1min, rc1max
                cnvslt(i, j, k, cnt) = 0.0_wp
             end do
          end do
       end do
    end do
    !$omp end target teams distribute parallel do
    do cnt = 1, maxdst
       tagslt = sltlist(cnt)

       call calc_spline_molecule(tagslt, stmax, splval(:,:,1:stmax), grdval(:,1:stmax))
       !$omp target teams distribute parallel do &
       !$omp&   map(to: splval, grdval) &
       !$omp&   map(alloc: mol_begin_index, charge, cnvslt, rcpslt)
       do sid = 1, stmax
!         ati = specatm(sid, tagslt)
          ati = mol_begin_index(tagslt) + (sid - 1)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             do cg2 = 0, splodr - 1
                do cg1 = 0, splodr - 1
                   rc1 = modulo(grdval(1, sid) - cg1, ms1max)
                   rc2 = modulo(grdval(2, sid) - cg2, ms2max)
                   rc3 = modulo(grdval(3, sid) - cg3, ms3max)
                   factor = chr * splval(cg1, 1, sid) * splval(cg2, 2, sid) &
                        * splval(cg3, 3, sid)
                   !$omp atomic update
                   cnvslt(rc1, rc2, rc3, cnt) = cnvslt(rc1, rc2, rc3, cnt) + factor
                end do
             end do
          end do
       end do
       !$omp end target teams distribute parallel do

       call fft_rtc(handle_r2c, cnvslt(:,:,:,cnt), rcpslt)                 ! 3D-FFT

       ! see the comment at the equivalent computation in
       ! recpcal_prepare_solute_refs: computed directly on the device
       ! via a reduction so rcpslt/engfac never need to leave the GPU.
       !
       ! original form is:
       ! 0.5 * sum(engfac(:, :, :) * real(rcpslt_c(:, :, :)) * conjg(rcpslt_c(:, :, :)))
       ! where rcpslt_c(rc1, rc2, rc3) = conjg(rcpslt_buf(ms1max - rc1, ms2max - rc2, ms3max - rc3))
       ! Here we use symmetry of engfac to calculate efficiently: every
       ! plane gets full weight, except rc1==0 (and rc1==ccemax when
       ! ms1max is even), which get half weight.
       self_energy = 0.0_wp
       !$omp target teams distribute parallel do collapse(3) reduction(+:self_energy) &
       !$omp&   map(alloc: engfac, rcpslt)
       do k = rc3min, rc3max
          do j = rc2min, rc2max
             do i = rc1min, ccemax
                if (i == 0 .or. (ms1max_even .and. i == ccemax)) then
                   self_energy = self_energy + &
                        0.5_wp * engfac(i, j, k) * real(rcpslt(i, j, k) * conjg(rcpslt(i, j, k)), wp)
                else
                   self_energy = self_energy + &
                        engfac(i, j, k) * real(rcpslt(i, j, k) * conjg(rcpslt(i, j, k)), wp)
                end if
             end do
          end do
       end do
       !$omp end target teams distribute parallel do
       solute_self_energy(cnt) = self_energy

       ! see the comment at the equivalent loop in recpcal_prepare_solute_refs
       !$omp target teams distribute parallel do collapse(3) map(alloc: engfac, rcpslt)
       do k = rc3min, rc3max
          do j = rc2min, rc2max
             do i = rc1min, ccemax
                rcpslt(i, j, k) = engfac(i, j, k) * rcpslt(i, j, k)
             end do
          end do
       end do
       !$omp end target teams distribute parallel do

       call fft_ctr(handle_c2r, rcpslt, cnvslt(:,:,:,cnt))                    ! 3D-FFT

    end do
  end subroutine recpcal_prepare_solute

  subroutine recpcal_prepare_solute_refs(tagslt, maxdst)
    use engmain, only: ms1max, ms2max, ms3max, sitepos, invcl, numsite, splodr, charge, mol_begin_index
    use fft_iface, only: fft_ctr, fft_rtc
    implicit none
    integer, intent(in) :: tagslt, maxdst
    integer :: i, j, k, cnt
    integer :: rc1, rc2, rc3, sid, ati, cg1, cg2, cg3, stmax
    real(wp) :: factor, chr
    real(wp) :: self_energy
    logical :: ms1max_even
    real(wp), allocatable, save :: splval(:,:,:)
    integer, allocatable, save :: grdval(:,:)
    logical, save :: initialized = .false.

    stmax = numsite(tagslt)
    if(.not. initialized) then
       allocate( splval(0:splodr-1, 3, stmax), grdval(3, stmax) )
       initialized = .true.
    end if
    ms1max_even = (mod(ms1max, 2) == 0)
    !$omp target teams distribute parallel do collapse(4) map(alloc: cnvslt)
    do cnt = 1, maxdst
       do k = rc3min, rc3max
          do j = rc2min, rc2max
             do i = rc1min, rc1max
                cnvslt(i, j, k, cnt) = 0.0_wp
             end do
          end do
       end do
    end do
    !$omp end target teams distribute parallel do
    do cnt = 1, maxdst

       call calc_spline_molecule_refs(tagslt, cnt, stmax, &
            splval(:,:,1:stmax), grdval(:,1:stmax))
       !$omp target teams distribute parallel do &
       !$omp&   map(to: splval, grdval) &
       !$omp&   map(alloc: mol_begin_index, charge, cnvslt, rcpslt)
       do sid = 1, stmax
          ! ati = specatm(sid, tagslt)
          ati = mol_begin_index(tagslt) + (sid - 1)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             do cg2 = 0, splodr - 1
                do cg1 = 0, splodr - 1
                   rc1 = modulo(grdval(1, sid) - cg1, ms1max)
                   rc2 = modulo(grdval(2, sid) - cg2, ms2max)
                   rc3 = modulo(grdval(3, sid) - cg3, ms3max)
                   factor = chr * splval(cg1, 1, sid) * splval(cg2, 2, sid) &
                        * splval(cg3, 3, sid)
                   !$omp atomic update
                   cnvslt(rc1, rc2, rc3, cnt) = &
                        cnvslt(rc1, rc2, rc3, cnt) + factor
                end do
             end do
          end do
       end do
       !$omp end target teams distribute parallel do

       call fft_rtc(handle_r2c, cnvslt(:,:,:,cnt), rcpslt)

       ! see the comment at the equivalent computation in
       ! recpcal_prepare_solute: computed directly on the device via a
       ! reduction so rcpslt/engfac never need to leave the GPU, instead
       ! of syncing rcpslt back to the host on every one of the
       ! (potentially thousands of) insertion trials in this loop.
       self_energy = 0.0_wp
       !$omp target teams distribute parallel do collapse(3) reduction(+:self_energy) &
       !$omp&   map(alloc: engfac, rcpslt)
       do k = rc3min, rc3max
          do j = rc2min, rc2max
             do i = rc1min, ccemax
                if (i == 0 .or. (ms1max_even .and. i == ccemax)) then
                   self_energy = self_energy + &
                        0.5_wp * engfac(i, j, k) * real(rcpslt(i, j, k) * conjg(rcpslt(i, j, k)), wp)
                else
                   self_energy = self_energy + &
                        engfac(i, j, k) * real(rcpslt(i, j, k) * conjg(rcpslt(i, j, k)), wp)
                end if
             end do
          end do
       end do
       !$omp end target teams distribute parallel do
       solute_self_energy(cnt) = self_energy

       ! see the comment at the equivalent loop in recpcal_prepare_solute
       !$omp target teams distribute parallel do collapse(3) map(alloc: engfac, rcpslt)
       do k = rc3min, rc3max
          do j = rc2min, rc2max
             do i = rc1min, ccemax
                rcpslt(i, j, k) = engfac(i, j, k) * rcpslt(i, j, k)
             end do
          end do
       end do
       !$omp end target teams distribute parallel do

       call fft_ctr(handle_c2r, rcpslt, cnvslt(:,:,:,cnt))

    end do
  end subroutine recpcal_prepare_solute_refs

  subroutine calc_spline_molecule(imol, stmax, store_spline, store_grid)
    use engmain, only: ms1max, ms2max, ms3max, splodr, specatm, sitepos, invcl
    use spline, only: spline_values_all
    implicit none

    integer, intent(in) :: imol, stmax
    real(wp), intent(out) :: store_spline(0:splodr-1, 3, 1:stmax)
    integer, intent(out) :: store_grid(3, 1:stmax)
    
    integer :: sid, ati, rcimax, m, k, rci
    real(wp) :: xst(3), inm(3)
    real(wp) :: factor, u

    do sid = 1, stmax
       ati = specatm(sid, imol)
       xst(:) = sitepos(:, ati)
       do k = 1, 3
          factor = dot_product(invcl(k,:), xst(:))
          factor = factor - floor(factor)
          inm(k) = factor
       end do
       do m = 1, 3
          if (m == 1) rcimax = ms1max
          if (m == 2) rcimax = ms2max
          if (m == 3) rcimax = ms3max
          factor = inm(m) * real(rcimax, wp)
          rci = int(factor)
          u = factor - real(rci, wp)
          call spline_values_all(u, store_spline(:, m, sid))
          store_grid(m, sid) = rci
       end do
    end do
  end subroutine calc_spline_molecule

  subroutine calc_spline_molecule_refs(imol, cntdst, stmax, store_spline, store_grid)
    use engmain, only: ms1max, ms2max, ms3max, splodr, specatm, sitepos, invcl
    use spline, only: spline_values_all
    implicit none

    integer, intent(in) :: imol, cntdst, stmax
    real(wp), intent(out) :: store_spline(0:splodr-1, 3, 1:stmax)
    integer, intent(out) :: store_grid(3, 1:stmax)

    integer :: sid, ati, ati_ext, rcimax, m, k, rci
    real(wp) :: xst(3), inm(3)
    real(wp) :: factor, u

    do sid = 1, stmax
       ati = specatm(sid, imol)
       ati_ext = ati + (cntdst - 1) * stmax
       xst(:) = sitepos(:, ati_ext)
       do k = 1, 3
          factor = dot_product(invcl(k,:), xst(:))
          factor = factor - floor(factor)
          inm(k) = factor
       end do
       do m = 1, 3
          if (m == 1) rcimax = ms1max
          if (m == 2) rcimax = ms2max
          if (m == 3) rcimax = ms3max
          factor = inm(m) * real(rcimax, wp)
          rci = int(factor)
          u = factor - real(rci, wp)
          call spline_values_all(u, store_spline(:, m, sid))
          store_grid(m, sid) = rci
       end do
    end do
  end subroutine calc_spline_molecule_refs

  function recpcal_self_energy(cnt) result(pairep)
    implicit none
    integer, intent(in) :: cnt
    real(wp) :: pairep

    pairep = solute_self_energy(cnt)
  end function recpcal_self_energy

  function recpcal_self_energy_refs(cnt) result(pairep)
    implicit none
    integer, intent(in) :: cnt
    real(wp) :: pairep

    pairep = solute_self_energy(cnt)
  end function recpcal_self_energy_refs

  subroutine recpcal_energy_soln(sltlist, maxdst, tagpt, slvmax, uvengy)
    use engmain, only: ms1max, ms2max, ms3max, splodr, numsite, sluvid, charge, mol_begin_index
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: sltlist(:), maxdst, tagpt(:), slvmax
    real(wp), intent(inout) :: uvengy(:, :)

    real(wp) :: pairep
    integer :: cg1, cg2, cg3, i, k, cnt, tagslt
    integer :: rc1, rc2, rc3, ptrnk, sid, ati, svi, stmax
    real(wp) :: fac1, fac2, fac3, chr
    integer :: grid1

    !$omp target teams distribute parallel do collapse(2) &
    !$omp&   map(alloc: uvengy, mol_begin_index, tagpt, sltlist, charge, numsite, sluvid, slvtag, splslv, grdslv, cnvslt)
    do cnt = 1, maxdst
    do k = 1, slvmax
       tagslt = sltlist(cnt)
       i = tagpt(k)
       if (i == tagslt) cycle

       pairep = 0.0_wp
       svi = slvtag(i)
       if (svi <= 0) stop  ! call halt_with_error('rcp_cns')
       stmax = numsite(i)
       do sid = 1, stmax
          ptrnk = svi + sid - 1
          ati = mol_begin_index(i) + (sid - 1) ! = specatm(sid, i)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             fac1 = chr * splslv(cg3, 3, ptrnk)
             rc3 = modulo(grdslv(3, ptrnk) - cg3, ms3max)
             do cg2 = 0, splodr - 1
                fac2 = fac1 * splslv(cg2, 2, ptrnk)
                rc2 = modulo(grdslv(2, ptrnk) - cg2, ms2max)
                grid1 = grdslv(1, ptrnk)
                if (grid1 >= splodr-1 .and. grid1 < ms1max) then
                   ! "!$acc loop seq" had no OpenMP equivalent; a loop with
                   ! no directive is already executed sequentially by the
                   ! thread that owns this (cnt, k) iteration.
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = grid1 - cg1
                      pairep = pairep + fac3 * cnvslt(rc1, rc2, rc3, cnt)
                   enddo
                else
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = mod(grid1 + ms1max - cg1, ms1max) ! speedhack
                      pairep = pairep + fac3 * cnvslt(rc1, rc2, rc3, cnt)
                   end do
                endif
             end do
          end do
       end do
       uvengy(k, cnt) = uvengy(k, cnt) + pairep
    end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine recpcal_energy_soln

  subroutine recpcal_energy_refs(tagslt, maxdst, slvmax, uvengy)
    use engmain, only: ms1max, ms2max, ms3max, splodr, numsite, sluvid, charge, mol_begin_index
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: tagslt, maxdst, slvmax
    real(wp), intent(inout) :: uvengy(:, :)

    real(wp) :: pairep
    integer :: cg1, cg2, cg3, i, k, cnt
    integer :: rc1, rc2, rc3, ptrnk, sid, ati, svi, stmax
    real(wp) :: fac1, fac2, fac3, chr
    integer :: grid1
    complex(wp) :: rcpt

    if (sluvid(tagslt) == 0) stop  ! call halt_with_error('rcp_fst')

    !$omp target teams distribute parallel do collapse(2) &
    !$omp&   map(alloc: uvengy, mol_begin_index, charge, numsite, sluvid, slvtag, splslv, grdslv, cnvslt)
    do cnt = 1, maxdst
    do i = 1, slvmax

       pairep = 0.0_wp
       svi = slvtag(i)
       if (svi <= 0) stop  ! call halt_with_error('rcp_cns')
       stmax = numsite(i)
       do sid = 1, stmax
          ptrnk = svi + sid - 1
          ati = mol_begin_index(i) + (sid - 1) ! = specatm(sid, i)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             fac1 = chr * splslv(cg3, 3, ptrnk)
             rc3 = modulo(grdslv(3, ptrnk) - cg3, ms3max)
             do cg2 = 0, splodr - 1
                fac2 = fac1 * splslv(cg2, 2, ptrnk)
                rc2 = modulo(grdslv(2, ptrnk) - cg2, ms2max)
                grid1 = grdslv(1, ptrnk)
                if (grid1 >= splodr-1 .and. grid1 < ms1max) then
                   ! see recpcal_energy_soln: no directive needed here,
                   ! this loop already runs sequentially within its thread.
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = grid1 - cg1
                      pairep = pairep &
                           + fac3 * cnvslt(rc1, rc2, rc3, cnt)
                   enddo
                else
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = mod(grid1 + ms1max - cg1, ms1max) ! speedhack
                      pairep = pairep &
                           + fac3 * cnvslt(rc1, rc2, rc3, cnt)
                   end do
                endif
             end do
          end do
       end do
       uvengy(i, cnt) = uvengy(i, cnt) + pairep
    end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine recpcal_energy_refs

end module reciprocal
