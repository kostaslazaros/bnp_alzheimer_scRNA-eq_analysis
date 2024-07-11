import scanpy as sc
import hdf5plugin
from pathlib import Path
import os


def get_subdirs_in_file(root_dir="./data/1.GSE_RAW", filename="meta_data.txt"):
    entries = os.listdir(root_dir)
    sub_dirs = [
        entry for entry in entries if os.path.isdir(os.path.join(root_dir, entry))
    ]
    with open(filename, "w", encoding="utf-8") as fil:
        fil.write("\n".join(sub_dirs))


def parse_meta_data(metadata_file):
    batch_data = []
    with open(metadata_file, encoding="utf-8") as fil:
        for line in fil.readlines():
            line = line.strip()
            if line.startswith("#"):
                continue
            if len(line) < 1:
                continue
            fdir, condition = line.split()
            batch_name = get_batch_name(fdir)
            batch_data.append((fdir, condition, batch_name))
    return batch_data


def create_h5ad(raw_path, h5ad_path, fdir, condition, batch_name):
    print(f"Reading {fdir} ... ", end="")
    adata = sc.read_10x_h5(raw_path / f"{fdir}.h5")
    adata.var_names_make_unique()
    print("ok. ", end="")
    adata.obs["batch"] = batch_name
    adata.obs["condition"] = condition
    filename = h5ad_path / f"{batch_name}_{condition}.h5ad"
    print(f"Saving to: {filename} ... ", end="")
    adata.write_h5ad(
        filename,
        compression=hdf5plugin.FILTERS["zstd"],
        compression_opts=hdf5plugin.Zstd(clevel=5).filter_options,
    )
    print("Done!!!")


def get_batch_name(fdir: str) -> str:
    batch_name, *_ = fdir.split("_")
    return batch_name


def main(metadata_file, data_root_path):
    path = Path(data_root_path)
    print("Running parse_meta_data")
    batch_data = parse_meta_data(metadata_file)
    for batch in batch_data:
        fdir, condition, batch_name = batch
        create_h5ad(path / "1.GSE_RAW", path / "2.H5AD_standalone", fdir, condition, batch_name)


if __name__ == "__main__":
    main("./data/meta_data.txt", "./data")
