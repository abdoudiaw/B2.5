#!/usr/bin/env python3
"""Validate SOLPS-ITER--SOLSTICE EIRENE training-event NetCDF files.

Project author and maintainer: Abdou Diaw
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import netCDF4
import numpy as np


CORE_OUTPUTS = (
    "eirbra_sni",
    "eirbra_smo",
    "eirbra_see",
    "eirbra_sei",
)
REQUIRED_INPUTS = ("braeir_dni", "braeir_te", "braeir_ti", "braeir_vol")
REQUIRED_METADATA = ("b2_tflux", "b2_flux_scale", "b2_crcstra")
REQUIRED_V3_METADATA = ("eirene_index_npoint",)
REQUIRED_WNEUTRALS = (
    "wneutrals_dab2",
    "wneutrals_dmb2",
    "wneutrals_tab2",
    "wneutrals_tmb2",
    "wneutrals_pfluxa",
    "wneutrals_rfluxa",
    "wneutrals_wldna",
    "wneutrals_wldpp",
    "wneutrals_wlpump",
    "cestim_pdena",
)
AVERAGED_WNEUTRALS = (
    "wneutrals_dab2",
    "wneutrals_tab2",
    "wneutrals_pfluxa",
    "wneutrals_rfluxa",
    "wneutrals_pefluxa",
    "wneutrals_refluxa",
    "wneutrals_dmb2",
    "wneutrals_tmb2",
    "wneutrals_pfluxm",
    "wneutrals_rfluxm",
    "wneutrals_pefluxm",
    "wneutrals_refluxm",
    "wneutrals_dib2",
    "wneutrals_tib2",
    "wneutrals_emiss",
    "wneutrals_emissmol",
    "wneutrals_srcml",
    "wneutrals_edissml",
    "wneutrals_wldnek",
    "wneutrals_wldnep",
    "wneutrals_wldna",
    "wneutrals_ewlda",
    "wneutrals_wldnm",
    "wneutrals_ewldm",
    "wneutrals_wldra",
    "wneutrals_wldrm",
    "wneutrals_wldpp",
    "wneutrals_wldpa",
    "wneutrals_wldpm",
    "wneutrals_wldpeb",
    "wneutrals_wldspt",
    "wneutrals_wldspta",
    "wneutrals_wldsptm",
    "wneutrals_eneutrad",
    "wneutrals_emolrad",
    "wneutrals_eionrad",
    "wneutrals_wlabsrp",
    "wneutrals_wldna_res",
    "wneutrals_wldnm_res",
    "wneutrals_ewlda_res",
    "wneutrals_ewldm_res",
    "wneutrals_ewldt_res",
    "wneutrals_ewldea_res",
    "wneutrals_ewldem_res",
    "wneutrals_ewldrp_res",
    "wneutrals_ewldmr_res",
    "wneutrals_wldspt_res",
    "wneutrals_wldspta_res",
    "wneutrals_wldsptm_res",
    "wneutrals_wlpump",
    "wneutrals_wlpump_res",
    "cestim_pdena",
    "cestim_pdenm",
    "cestim_pdeni",
    "cestim_edena",
    "cestim_edenm",
    "cestim_edeni",
    "cestim_vxdena",
    "cestim_vxdenm",
    "cestim_vxdeni",
    "cestim_vydena",
    "cestim_vydenm",
    "cestim_vydeni",
    "cestim_vzdena",
    "cestim_vzdenm",
    "cestim_vzdeni",
)


def event_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    return sorted(path.glob("eirene_training_v*_b2call_*.nc"))


def validate_file(
    path: Path,
) -> tuple[int, str, int, int, int, dict[str, np.ndarray]]:
    arrays: dict[str, np.ndarray] = {}
    with netCDF4.Dataset(path) as dataset:
        if dataset.getncattr("schema_name") != "solps_eirene_training_event":
            raise ValueError(f"{path}: unexpected schema_name")
        schema_version = str(dataset.getncattr("schema_version"))
        if schema_version not in ("1.0.0", "2.0.0", "3.0.0"):
            raise ValueError(f"{path}: unexpected schema_version")

        required = CORE_OUTPUTS + REQUIRED_INPUTS + REQUIRED_METADATA
        if schema_version in ("2.0.0", "3.0.0"):
            required += REQUIRED_WNEUTRALS
        if schema_version == "3.0.0":
            required += REQUIRED_V3_METADATA
            if dataset.getncattr("braeir_capture_phase") != (
                "immediately before eirene_eirsrt"
            ):
                raise ValueError(f"{path}: unexpected BRAEIR capture phase")
            if not str(dataset.getncattr("braeir_input_indexing")).startswith(
                "B2 indexing"
            ):
                raise ValueError(f"{path}: unexpected BRAEIR input indexing")
        missing = [name for name in required if name not in dataset.variables]
        if missing:
            raise ValueError(f"{path}: missing variables: {', '.join(missing)}")

        for name, variable in dataset.variables.items():
            if variable.dtype.kind not in "iufc":
                continue
            values = np.asarray(variable[:])
            if not np.all(np.isfinite(values)):
                count = int(values.size - np.count_nonzero(np.isfinite(values)))
                raise ValueError(f"{path}: {name} has {count} non-finite values")
            if (
                name in CORE_OUTPUTS
                or name.startswith("braeir_")
                or name in AVERAGED_WNEUTRALS
            ):
                arrays[name] = values

        b2_call = int(dataset.getncattr("b2_call_index"))
        kind = str(dataset.getncattr("event_kind"))
        repeat_index = int(dataset.getncattr("eirene_repeat_index"))
        repeat_count = int(dataset.getncattr("eirene_repeat_count"))
        result_used = int(
            dataset.getncattr("eirene_result_used_by_b2")
            if schema_version == "3.0.0"
            else -1
        )
        strata = "".join(dataset.variables["b2_crcstra"][:].astype(str))
        print(
            f"PASS {path.name}: kind={kind} b2_call={b2_call} "
            f"repeat={repeat_index}/{repeat_count} used_by_b2={result_used} "
            f"strata={strata!r} variables={len(dataset.variables)}"
        )
        for name in CORE_OUTPUTS:
            values = arrays[name]
            nonzero = np.count_nonzero(values)
            print(
                f"  {name:12s} shape={str(values.shape):18s} "
                f"nonzero={nonzero / values.size:7.3%} "
                f"min={values.min(): .5e} max={values.max(): .5e}"
            )

    return b2_call, kind, repeat_index, repeat_count, result_used, arrays


def validate_repeat_metadata(
    events: list[
        tuple[Path, int, str, int, int, int, dict[str, np.ndarray]]
    ],
) -> None:
    by_call: dict[
        int, list[tuple[Path, str, int, int, int, dict[str, np.ndarray]]]
    ] = defaultdict(list)
    for path, b2_call, kind, index, count, used, arrays in events:
        by_call[b2_call].append((path, kind, index, count, used, arrays))

    for b2_call, group in sorted(by_call.items()):
        singles = [item for item in group if item[1] == "single_call"]
        if not singles or any(item[4] < 0 for item in singles):
            continue
        counts = {item[3] for item in singles}
        indices = sorted(item[2] for item in singles)
        if len(counts) != 1 or indices != list(range(1, len(singles) + 1)):
            raise ValueError(f"B2 call {b2_call}: inconsistent repeat metadata")
        count = counts.pop()
        if count != len(singles):
            raise ValueError(
                f"B2 call {b2_call}: repeat_count={count}, "
                f"but found {len(singles)} single-call events"
            )

        averages = [item for item in group if item[1] == "average_used_by_b2"]
        used = [item for item in group if item[4] == 1]
        if len(used) != 1:
            raise ValueError(
                f"B2 call {b2_call}: expected exactly one result used by B2"
            )
        if averages and used[0][1] != "average_used_by_b2":
            raise ValueError(
                f"B2 call {b2_call}: averaged result is not marked used by B2"
            )


def validate_averages(
    events: list[
        tuple[Path, int, str, int, int, int, dict[str, np.ndarray]]
    ],
    rtol: float,
    atol: float,
) -> None:
    by_call: dict[int, list[tuple[Path, str, dict[str, np.ndarray]]]] = (
        defaultdict(list)
    )
    for path, b2_call, kind, _, _, _, arrays in events:
        by_call[b2_call].append((path, kind, arrays))

    for b2_call, group in sorted(by_call.items()):
        singles = [item for item in group if item[1] == "single_call"]
        averages = [item for item in group if item[1] == "average_used_by_b2"]
        if not averages:
            continue
        if len(averages) != 1 or not singles:
            raise ValueError(f"B2 call {b2_call}: incomplete averaging event set")

        average_arrays = averages[0][2]
        averaged_names = CORE_OUTPUTS + tuple(
            name for name in AVERAGED_WNEUTRALS if name in average_arrays
        )
        for name in averaged_names:
            expected = np.mean([item[2][name] for item in singles], axis=0)
            actual = average_arrays[name]
            if not np.allclose(actual, expected, rtol=rtol, atol=atol):
                error = float(np.max(np.abs(actual - expected)))
                raise ValueError(
                    f"B2 call {b2_call}: {name} average mismatch; "
                    f"maximum absolute error={error:.6e}"
                )
        print(
            f"PASS B2 call {b2_call}: reproduced {len(averaged_names)} "
            "averaged EIRBRA/WNEUTRALS arrays"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path, help="event file or SOLPS run directory")
    parser.add_argument("--rtol", type=float, default=1.0e-12)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args()

    files = event_files(args.path)
    if not files:
        parser.error(f"no EIRENE training-event files found below {args.path}")

    events = []
    for path in files:
        event = validate_file(path)
        events.append((path, *event))
    validate_repeat_metadata(events)
    validate_averages(events, args.rtol, args.atol)
    print(f"PASS validated {len(files)} EIRENE training event(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
