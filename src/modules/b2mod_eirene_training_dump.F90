! SOLPS-ITER--SOLSTICE neutral-source coupling project
! Project author and maintainer: Abdou Diaw
!
! Observation-only NetCDF writer for the B2.5--EIRENE coupling seam.
! This module does not modify EIRENE tallies or B2.5 source terms.

module b2mod_eirene_training_dump
  use b2mod_types, only : R8
#ifdef B25_EIRENE
  use eirmod_braeir
  use eirmod_eirbra
  use eirmod_mpi, only : mpi_comm_rank, mpi_comm_world
#endif
  implicit none
  private

  public :: write_eirene_training_event

contains

  subroutine write_eirene_training_event(b2_call, repeat_index, &
      repeat_count, event_kind, tflux_b2, flux_scale_b2, crcstra_b2)
    implicit none
    integer, intent(in) :: b2_call, repeat_index, repeat_count
    integer, intent(in) :: event_kind
    real(kind=R8), intent(in) :: tflux_b2(:), flux_scale_b2(:)
    character(len=1), intent(in) :: crcstra_b2(:)

#if defined(B25_EIRENE) && !defined(NO_CDF)
#include <netcdf.inc>
    integer :: status, ncid, mpi_rank, mpi_error
    integer :: eir_x_id, eir_y_id, fluid_id, momentum_slot_id
    integer :: bra_x_id, bra_y_id, bra_fluid_id
    integer :: stratum_id, active_stratum_id, three_id
    integer :: atom_id, molecule_id, ion_id
    integer :: d2(2), d3(3), d4(4)
    character(len=256) :: filename
    character(len=32) :: event_label
    character(len=32) :: b25_hash, solps_hash
    character(len=40) :: eirene_hash
    character(len=32) :: timestamp
    character(len=8) :: date_string
    character(len=10) :: time_string
    character(len=5) :: zone_string
    integer :: date_values(8)
    character(len=32) :: get_B25_hash, get_SOLPS_hash
    character(len=40) :: get_Eir_hash
    external :: check_cdf_status, get_B25_hash, get_SOLPS_hash
    external :: get_Eir_hash

    call mpi_comm_rank(mpi_comm_world, mpi_rank, mpi_error)
    if (mpi_rank /= 0) return

    if (.not.allocated(SNI)) then
      write(*,*) 'EIRENE training dump requested before SNI allocation'
      return
    end if

    if (event_kind == 0) then
      event_label = 'single_call'
    else
      event_label = 'average_used_by_b2'
    end if
    write(filename,'(a,i8.8,a,a,a,i4.4,a)') &
        'eirene_training_v1_b2call_', b2_call, '_', &
        trim(event_label), '_', repeat_index, '.nc'

    status = nf_create(trim(filename), NF_CLOBBER, ncid)
    call check_cdf_status(status)

    status = nf_def_dim(ncid, 'eir_x', size(SNI,1), eir_x_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'eir_y', size(SNI,2), eir_y_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'fluid_species', size(SNI,3), fluid_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'momentum_slot', size(SMO,3), &
        momentum_slot_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'stratum_slot', size(SNI,4), stratum_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'active_stratum', size(tflux_b2), &
        active_stratum_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'three', 3, three_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'bra_x', size(DNIB,1), bra_x_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'bra_y', size(DNIB,2), bra_y_id)
    call check_cdf_status(status)
    status = nf_def_dim(ncid, 'bra_fluid_species', size(DNIB,3), &
        bra_fluid_id)
    call check_cdf_status(status)

    atom_id = -1
    molecule_id = -1
    ion_id = -1
    if (size(SRCCRFC,1) > 0) then
      status = nf_def_dim(ncid, 'atom_species', size(SRCCRFC,1), atom_id)
      call check_cdf_status(status)
    end if
    if (size(SNA_PAML,3) > 0) then
      status = nf_def_dim(ncid, 'molecule_species', size(SNA_PAML,3), &
          molecule_id)
      call check_cdf_status(status)
    end if
    if (size(SNA_PAIO,3) > 0) then
      status = nf_def_dim(ncid, 'test_ion_species', size(SNA_PAIO,3), &
          ion_id)
      call check_cdf_status(status)
    end if

    b25_hash = get_B25_hash()
    solps_hash = get_SOLPS_hash()
    eirene_hash = get_Eir_hash()
    call date_and_time(date_string, time_string, zone_string, date_values)
    timestamp = date_string//'T'//time_string//zone_string

    call put_global_text(ncid, 'schema_name', &
        'solps_eirene_training_event')
    call put_global_text(ncid, 'schema_version', '1.0.0')
    call put_global_text(ncid, 'seam', &
        'after eirene_eirsrt; before B2 mapping, scaling, and linearization')
    call put_global_text(ncid, 'event_kind', trim(event_label))
    call put_global_text(ncid, 'created_local', trim(timestamp))
    call put_global_text(ncid, 'b2_5_git', trim(b25_hash))
    call put_global_text(ncid, 'solps_iter_git', trim(solps_hash))
    call put_global_text(ncid, 'eirene_git', trim(eirene_hash))
    call put_global_text(ncid, 'raw_value_policy', &
        'EIRBRA values are unscaled and unlinearized EIRENE estimators')
    call put_global_text(ncid, 'contract_scope', &
        'BRAEIR and EIRBRA; WNEUTRALS is not included in schema v1')
    call put_global_text(ncid, 'smo_slot_layout', &
        'parallel species, radial species, diamagnetic species')
    call put_global_text(ncid, 'stratum_slot_policy', &
        'EIRBRA allocation may include one reserved time-stratum slot')
    call put_global_int(ncid, 'b2_call_index', b2_call)
    call put_global_int(ncid, 'eirene_repeat_index', repeat_index)
    call put_global_int(ncid, 'eirene_repeat_count', repeat_count)

    status = nf_enddef(ncid)
    call check_cdf_status(status)

    d4 = (/eir_x_id, eir_y_id, fluid_id, stratum_id/)
    call put_real(ncid, 'eirbra_sni', SNI, d4, 4, &
        'raw EIRENE particle source estimator')
    d4 = (/eir_x_id, eir_y_id, momentum_slot_id, stratum_id/)
    call put_real(ncid, 'eirbra_smo', SMO, d4, 4, &
        'raw EIRENE vector momentum source estimator')
    d3 = (/eir_x_id, eir_y_id, stratum_id/)
    call put_real(ncid, 'eirbra_see', SEE, d3, 3, &
        'raw EIRENE electron-energy source estimator')
    call put_real(ncid, 'eirbra_sei', SEI, d3, 3, &
        'raw EIRENE ion-energy source estimator')

    call put_real(ncid, 'eirbra_sne_pael', SNE_PAEL, d3, 3, &
        'electron-particle channel: atom-electron')
    call put_real(ncid, 'eirbra_sne_pmel', SNE_PMEL, d3, 3, &
        'electron-particle channel: molecule-electron')
    call put_real(ncid, 'eirbra_sne_piel', SNE_PIEL, d3, 3, &
        'electron-particle channel: test-ion-electron')

    if (atom_id >= 0) then
      d4 = (/eir_x_id, eir_y_id, atom_id, stratum_id/)
      call put_real(ncid, 'eirbra_sna_paat', SNA_PAAT, d4, 4, &
          'atom particle channel PAAT')
      call put_real(ncid, 'eirbra_sna_pmat', SNA_PMAT, d4, 4, &
          'atom particle channel PMAT')
      call put_real(ncid, 'eirbra_sna_piat', SNA_PIAT, d4, 4, &
          'atom particle channel PIAT')
    end if
    if (molecule_id >= 0) then
      d4 = (/eir_x_id, eir_y_id, molecule_id, stratum_id/)
      call put_real(ncid, 'eirbra_sna_paml', SNA_PAML, d4, 4, &
          'molecule particle channel PAML')
      call put_real(ncid, 'eirbra_sna_pmml', SNA_PMML, d4, 4, &
          'molecule particle channel PMML')
      call put_real(ncid, 'eirbra_sna_piml', SNA_PIML, d4, 4, &
          'molecule particle channel PIML')
    end if
    if (ion_id >= 0) then
      d4 = (/eir_x_id, eir_y_id, ion_id, stratum_id/)
      call put_real(ncid, 'eirbra_sna_paio', SNA_PAIO, d4, 4, &
          'test-ion particle channel PAIO')
      call put_real(ncid, 'eirbra_sna_pmio', SNA_PMIO, d4, 4, &
          'test-ion particle channel PMIO')
      call put_real(ncid, 'eirbra_sna_piio', SNA_PIIO, d4, 4, &
          'test-ion particle channel PIIO')
    end if

    d4 = (/eir_x_id, eir_y_id, fluid_id, stratum_id/)
    call put_real(ncid, 'eirbra_sni_papl', SNI_PAPL, d4, 4, &
        'plasma particle channel PAPL')
    call put_real(ncid, 'eirbra_sni_pmpl', SNI_PMPL, d4, 4, &
        'plasma particle channel PMPL')
    call put_real(ncid, 'eirbra_sni_pipl', SNI_PIPL, d4, 4, &
        'plasma particle channel PIPL')
    call put_real(ncid, 'eirbra_sni_pppl', SNI_PPPL, d4, 4, &
        'plasma particle channel PPPL')
    call put_real(ncid, 'eirbra_smo_mapl', SMO_MAPL, d4, 4, &
        'plasma momentum channel MAPL')
    call put_real(ncid, 'eirbra_smo_mmpl', SMO_MMPL, d4, 4, &
        'plasma momentum channel MMPL')
    call put_real(ncid, 'eirbra_smo_mipl', SMO_MIPL, d4, 4, &
        'plasma momentum channel MIPL')
    call put_real(ncid, 'eirbra_smo_mppl', SMO_MPPL, d4, 4, &
        'plasma momentum channel MPPL')

    d3 = (/eir_x_id, eir_y_id, stratum_id/)
    call put_real(ncid, 'eirbra_see_eael', SEE_EAEL, d3, 3, &
        'electron-energy channel EAEL')
    call put_real(ncid, 'eirbra_see_emel', SEE_EMEL, d3, 3, &
        'electron-energy channel EMEL')
    call put_real(ncid, 'eirbra_see_eiel', SEE_EIEL, d3, 3, &
        'electron-energy channel EIEL')
    call put_real(ncid, 'eirbra_see_epel', SEE_EPEL, d3, 3, &
        'electron-energy channel EPEL')
    call put_real(ncid, 'eirbra_sei_eapl', SEI_EAPL, d3, 3, &
        'ion-energy channel EAPL')
    call put_real(ncid, 'eirbra_sei_empl', SEI_EMPL, d3, 3, &
        'ion-energy channel EMPL')
    call put_real(ncid, 'eirbra_sei_eipl', SEI_EIPL, d3, 3, &
        'ion-energy channel EIPL')
    call put_real(ncid, 'eirbra_sei_eppl', SEI_EPPL, d3, 3, &
        'ion-energy channel EPPL')

    call put_real(ncid, 'eirbra_volsumn', VOLSUMN, &
        (/stratum_id/), 1, 'particle-source volume sum')
    call put_real(ncid, 'eirbra_volsumm', VOLSUMM, &
        (/three_id, stratum_id/), 2, 'momentum-source volume sum')
    call put_real(ncid, 'eirbra_volsumei', VOLSUMEI, &
        (/stratum_id/), 1, 'ion-energy-source volume sum')
    call put_real(ncid, 'eirbra_volsumee', VOLSUMEE, &
        (/stratum_id/), 1, 'electron-energy-source volume sum')
    call put_real(ncid, 'eirbra_srcstrn', SRCSTRN, &
        (/stratum_id/), 1, 'neutral source intensity by stratum')
    call put_real(ncid, 'eirbra_flxspci', FLXSPCI, &
        (/fluid_id, stratum_id/), 2, &
        'plasma-ion flux to recycling surfaces')
    if (atom_id >= 0) then
      call put_real(ncid, 'eirbra_srccrfc', SRCCRFC, &
          (/atom_id, stratum_id/), 2, 'source correction factors')
    end if

    call put_real(ncid, 'b2_tflux', tflux_b2, &
        (/active_stratum_id/), 1, 'tflux argument passed to EIRENE')
    call put_real(ncid, 'b2_flux_scale', flux_scale_b2, &
        (/active_stratum_id/), 1, 'B2 post-EIRENE flux scale')
    call put_text(ncid, 'b2_crcstra', crcstra_b2, active_stratum_id, &
        'B2 stratum type code')

    d3 = (/bra_x_id, bra_y_id, bra_fluid_id/)
    call put_real(ncid, 'braeir_dni', DNIB, d3, 3, &
        'plasma species density supplied to EIRENE')
    call put_real(ncid, 'braeir_vv', VVB, d3, 3, &
        'plasma velocity quantity VVB supplied to EIRENE')
    call put_real(ncid, 'braeir_uu', UUB, d3, 3, &
        'plasma velocity quantity UUB supplied to EIRENE')
    call put_real(ncid, 'braeir_ww', WWB, d3, 3, &
        'plasma velocity quantity WWB supplied to EIRENE')
    call put_real(ncid, 'braeir_up', UPB, d3, 3, &
        'parallel velocity supplied to EIRENE')
    call put_real(ncid, 'braeir_fnix', FNIXB, d3, 3, &
        'x plasma particle flux supplied to EIRENE')
    call put_real(ncid, 'braeir_fniy', FNIYB, d3, 3, &
        'y plasma particle flux supplied to EIRENE')
    call put_real(ncid, 'braeir_vparx', VPARXB, d3, 3, &
        'x parallel-velocity projection supplied to EIRENE')
    call put_real(ncid, 'braeir_vpary', VPARYB, d3, 3, &
        'y parallel-velocity projection supplied to EIRENE')
    call put_real(ncid, 'braeir_vradx', VRADXB, d3, 3, &
        'x radial-velocity projection supplied to EIRENE')
    call put_real(ncid, 'braeir_vrady', VRADYB, d3, 3, &
        'y radial-velocity projection supplied to EIRENE')
    call put_real(ncid, 'braeir_uudia', UUDIAB, d3, 3, &
        'diamagnetic velocity quantity UUDIAB supplied to EIRENE')
    call put_real(ncid, 'braeir_vvdia', VVDIAB, d3, 3, &
        'diamagnetic velocity quantity VVDIAB supplied to EIRENE')
    call put_real(ncid, 'braeir_zi', ZIB, d3, 3, &
        'plasma charge quantity supplied to EIRENE')

    d2 = (/bra_x_id, bra_y_id/)
    call put_real(ncid, 'braeir_te', TEB, d2, 2, &
        'electron temperature supplied to EIRENE')
    call put_real(ncid, 'braeir_ti', TIB, d2, 2, &
        'ion temperature supplied to EIRENE')
    call put_real(ncid, 'braeir_pr', PRB, d2, 2, &
        'plasma pressure quantity supplied to EIRENE')
    call put_real(ncid, 'braeir_rr', RRB, d2, 2, &
        'major radius supplied to EIRENE')
    call put_real(ncid, 'braeir_feix', FEIXB, d2, 2, &
        'x ion-energy flux supplied to EIRENE')
    call put_real(ncid, 'braeir_feiy', FEIYB, d2, 2, &
        'y ion-energy flux supplied to EIRENE')
    call put_real(ncid, 'braeir_feex', FEEXB, d2, 2, &
        'x electron-energy flux supplied to EIRENE')
    call put_real(ncid, 'braeir_feey', FEEYB, d2, 2, &
        'y electron-energy flux supplied to EIRENE')
    call put_real(ncid, 'braeir_vol', VOLB, d2, 2, &
        'cell volume supplied to EIRENE')
    call put_real(ncid, 'braeir_bfield', BFELDB, d2, 2, &
        'magnetic-field magnitude supplied to EIRENE')
    call put_real(ncid, 'braeir_bpol', BPOLB, d2, 2, &
        'poloidal magnetic field supplied to EIRENE')
    call put_real(ncid, 'braeir_brad', BRADB, d2, 2, &
        'radial magnetic field supplied to EIRENE')
    call put_real(ncid, 'braeir_btor', BTORB, d2, 2, &
        'toroidal magnetic field supplied to EIRENE')
    call put_real(ncid, 'braeir_deltae_parx', DELTAE_PARXB, d2, 2, &
        'electron parallel-energy correction x')
    call put_real(ncid, 'braeir_deltae_pary', DELTAE_PARYB, d2, 2, &
        'electron parallel-energy correction y')
    call put_real(ncid, 'braeir_deltae_radx', DELTAE_RADXB, d2, 2, &
        'electron radial-energy correction x')
    call put_real(ncid, 'braeir_deltae_rady', DELTAE_RADYB, d2, 2, &
        'electron radial-energy correction y')
    call put_real(ncid, 'braeir_deltai_parx', DELTAI_PARXB, d2, 2, &
        'ion parallel-energy correction x')
    call put_real(ncid, 'braeir_deltai_pary', DELTAI_PARYB, d2, 2, &
        'ion parallel-energy correction y')
    call put_real(ncid, 'braeir_deltai_radx', DELTAI_RADXB, d2, 2, &
        'ion radial-energy correction x')
    call put_real(ncid, 'braeir_deltai_rady', DELTAI_RADYB, d2, 2, &
        'ion radial-energy correction y')
    call put_real(ncid, 'braeir_delta_sheathx', DELTA_SHEATHXB, d2, 2, &
        'sheath-energy correction x')
    call put_real(ncid, 'braeir_delta_sheathy', DELTA_SHEATHYB, d2, 2, &
        'sheath-energy correction y')
    call put_real(ncid, 'braeir_aiso', AISOB, d2, 2, &
        'isotropy quantity supplied to EIRENE')
    call put_real(ncid, 'braeir_po', POB, d2, 2, &
        'electrostatic potential supplied to EIRENE')

    status = nf_close(ncid)
    call check_cdf_status(status)
    write(*,*) 'Wrote EIRENE training event: ', trim(filename)
#else
    if (b2_call == 0 .and. repeat_index == 1) then
      write(*,*) 'eirene_training_dump ignored: NetCDF/EIRENE unavailable'
    end if
#endif
  end subroutine write_eirene_training_event

#if defined(B25_EIRENE) && !defined(NO_CDF)
  subroutine put_global_text(ncid, name, value)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name, value
    integer :: status
    external :: check_cdf_status

    status = nf_put_att_text(ncid, NF_GLOBAL, trim(name), &
        len_trim(value), trim(value))
    call check_cdf_status(status)
  end subroutine put_global_text

  subroutine put_global_int(ncid, name, value)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid, value
    character(len=*), intent(in) :: name
    integer :: status
    external :: check_cdf_status

    status = nf_put_att_int(ncid, NF_GLOBAL, trim(name), NF_INT, 1, value)
    call check_cdf_status(status)
  end subroutine put_global_int

  subroutine put_real(ncid, name, values, dimids, ndim, long_name)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid, dimids(*), ndim
    character(len=*), intent(in) :: name, long_name
    real(kind=R8), intent(in) :: values(*)
    integer :: status, varid
    external :: check_cdf_status

    status = nf_redef(ncid)
    call check_cdf_status(status)
    status = nf_def_var(ncid, trim(name), NF_DOUBLE, ndim, dimids, varid)
    call check_cdf_status(status)
    status = nf_put_att_text(ncid, varid, 'long_name', &
        len_trim(long_name), trim(long_name))
    call check_cdf_status(status)
    status = nf_put_att_text(ncid, varid, 'units', 22, &
        'EIRENE coupling native')
    call check_cdf_status(status)
    status = nf_enddef(ncid)
    call check_cdf_status(status)
    status = nf_put_var_double(ncid, varid, values)
    call check_cdf_status(status)
  end subroutine put_real

  subroutine put_text(ncid, name, values, dimid, long_name)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid, dimid
    character(len=*), intent(in) :: name, long_name
    character(len=1), intent(in) :: values(*)
    integer :: status, varid, dims(1)
    external :: check_cdf_status

    dims(1) = dimid
    status = nf_redef(ncid)
    call check_cdf_status(status)
    status = nf_def_var(ncid, trim(name), NF_CHAR, 1, dims, varid)
    call check_cdf_status(status)
    status = nf_put_att_text(ncid, varid, 'long_name', &
        len_trim(long_name), trim(long_name))
    call check_cdf_status(status)
    status = nf_enddef(ncid)
    call check_cdf_status(status)
    status = nf_put_var_text(ncid, varid, values)
    call check_cdf_status(status)
  end subroutine put_text
#endif

end module b2mod_eirene_training_dump
