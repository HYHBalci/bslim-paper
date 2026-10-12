#!/usr/bin/env python3
"""Prepare all 100 public IHDP inputs before upload; no fitting or installation.

The historical row order is preserved: pinned CEVAE CSVs for replicates 1..10,
then the checksummed Johansson NPZ train/test concatenation for 11..100.
Use --offline to verify a completed input directory without network access.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import urllib.request

CEVAE_REVISION = "9081f863e24ce21bd34c8d6a41bf0edc7d1b65dd"
NPZ_FILES = {
    "train": ("ihdp_npci_1-100.train.npz", "750697c71b4f8d7a3aafff771b56a4ac4cd83ec649bf69afb04f8a5aee41a240"),
    "test": ("ihdp_npci_1-100.test.npz", "a70a8acbcc4e8deb677cc9bf9e9dabeb17caaa37cdbb1d7ba06be7ffb929c41c"),
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, path, expected=None):
    path = Path(path)
    if not path.exists():
        temporary = path.with_suffix(path.suffix + ".download")
        request = urllib.request.Request(url, headers={"User-Agent": "IHDP-reproducibility-bundle/8.3.5"})
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
            for block in iter(lambda: response.read(1024 * 1024), b""):
                output.write(block)
        if expected is not None and sha256(temporary) != expected:
            raise ValueError(f"Downloaded checksum mismatch: {url}; temporary file retained")
        temporary.replace(path)
    if expected is not None and sha256(path) != expected:
        raise ValueError(f"Existing input checksum mismatch; not overwritten: {path}")


def validate_matrix(matrix, np):
    if matrix.shape != (747, 30) or not np.isfinite(matrix).all():
        raise ValueError(f"Expected a finite 747 x 30 IHDP input; found {matrix.shape}")
    if not np.isin(matrix[:, 0], [0, 1]).all() or matrix[:, 0].sum() != 139:
        raise ValueError("Unexpected IHDP treatment coding or count")


def assemble(train, test, index, np):
    columns = [np.concatenate((train[key][:, index], test[key][:, index]))
               for key in ("t", "yf", "ycf", "mu0", "mu1")]
    covariates = np.concatenate((train["x"][:, :, index], test["x"][:, :, index]), axis=0)
    return np.column_stack(columns + [covariates])


def sorted_rows(matrix, np):
    return matrix[np.lexsort(tuple(matrix[:, j] for j in range(matrix.shape[1] - 1, -1, -1)))]


def verify_offline(directory):
    manifest = Path(directory) / "input_manifest.csv"
    if not manifest.is_file():
        raise ValueError(f"Missing input manifest: {manifest}")
    with manifest.open(newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    if len(records) != 100 or {int(row["replicate"]) for row in records} != set(range(1, 101)):
        raise ValueError("IHDP manifest must identify exactly the 100 distinct replicates")
    for row in records:
        expected_name = f'ihdp_npci_{int(row["replicate"])}.csv'
        if row["file"] != expected_name:
            raise ValueError("Unexpected filename in IHDP manifest")
        path = Path(directory) / expected_name
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise ValueError(f"Missing or modified prepared input: {path}")
    print(f"Verified all 100 prepared IHDP inputs in {directory}")


def prepare(directory, download_cache):
    import numpy as np

    directory = Path(directory)
    download_cache = Path(download_cache)
    directory.mkdir(parents=True, exist_ok=True)
    download_cache.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, (name, checksum) in NPZ_FILES.items():
        path = download_cache / name
        fetch(f"https://www.fredjo.com/files/{name}", path, checksum)
        paths[key] = path
    records = []
    with np.load(paths["train"]) as train, np.load(paths["test"]) as test:
        for replicate in range(1, 101):
            expected = assemble(train, test, replicate - 1, np)
            validate_matrix(expected, np)
            destination = directory / f"ihdp_npci_{replicate}.csv"
            if replicate <= 10:
                source = (f"https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/{CEVAE_REVISION}"
                          f"/datasets/IHDP/csv/ihdp_npci_{replicate}.csv")
                fetch(source, destination)
                observed = np.loadtxt(destination, delimiter=",")
                validate_matrix(observed, np)
                if not np.allclose(sorted_rows(observed, np), sorted_rows(expected, np), rtol=1e-12, atol=1e-10):
                    raise ValueError(f"CEVAE replicate {replicate} differs from checksummed NPZ contents")
                order = "Pinned CEVAE CSV, preserved"
            else:
                source = "https://www.fredjo.com/files/ihdp_npci_1-100.{train,test}.npz"
                if not destination.exists():
                    np.savetxt(destination, expected, delimiter=",", fmt="%.15g")
                observed = np.loadtxt(destination, delimiter=",")
                validate_matrix(observed, np)
                if not np.allclose(observed, expected, rtol=1e-12, atol=1e-10):
                    raise ValueError(f"Existing replicate {replicate} differs in values or row order; not overwritten")
                order = "NPZ train then test"
            records.append({"replicate": replicate, "file": destination.name,
                            "sha256": sha256(destination), "source": source, "row_order": order})
    with (directory / "input_manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (directory / "preparation.json").write_text(json.dumps({
        "cevae_revision": CEVAE_REVISION, "numpy_version": np.__version__,
        "npz_sha256": {key: sha256(path) for key, path in paths.items()},
        "replicates": 100, "shape": [747, 30],
        "row_order": "CEVAE 1..10, NPZ train then test 11..100; affects seeded folds",
    }, indent=2) + "\n", encoding="utf-8")
    verify_offline(directory)


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=root / "data" / "ihdp")
    parser.add_argument("--download-cache", type=Path, default=root / "data" / "_downloads")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    if args.offline:
        verify_offline(args.out)
    else:
        prepare(args.out, args.download_cache)


if __name__ == "__main__":
    main()
