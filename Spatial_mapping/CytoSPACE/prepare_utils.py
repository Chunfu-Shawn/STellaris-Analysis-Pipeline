import os
import sys
import pandas as pd
import numpy as np
import anndata
import gzip
import fsspec
import scipy
import logging
import pandas._libs.json as ujson
from urllib.parse import urlparse

logger = logging.getLogger("prepare_utils")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)

def get_scheme(path):
    pr = urlparse(path)
    if len(pr.scheme) <= 1:  # for file paths: /foo/bar/test.h5ad or C:/foo/bar/test.h5ad
        return "file"
    return pr.scheme

fsspec_kwargs = dict(s3_additional_kwargs=dict(ACL="bucket-owner-full-control"))

def get_fs(path):
    return fsspec.filesystem(get_scheme(path), **fsspec_kwargs)

def open_file(urlpath, mode="rb", compression=None):
    return fsspec.open(urlpath, mode=mode, compression=compression, **fsspec_kwargs)

def to_json(data, orient="values"):
    return ujson.dumps(data, double_precision=2, orient=orient)

SPATIAL_HELP = (
    "Directory containing 10x visium spatial data (tissue_hires_image.png, scalefactors_json.json, "
    + "and tissue_positions_list.csv) "
    + "or a directory containing `image.png`, `positions.image.csv` "
    + "with headers barcode, x, and y, and optionally `diameter.image.txt` containing spot diameter"
)

def __add_visium(adata, spatial_directory):
    scale_factors_path = os.path.join(spatial_directory, "scalefactors_json.json")
    tissue_hires_image_path = os.path.join(spatial_directory, "tissue_hires_image.png")
    tissue_positions_list_path = os.path.join(spatial_directory, "tissue_positions_list.csv")
    is_visium = True
    for path in [scale_factors_path, tissue_hires_image_path, tissue_positions_list_path]:
        if not os.path.exists(path):
            is_visium = False
            break
    if is_visium:
        import json

        with open(os.path.join(spatial_directory, "scalefactors_json.json"), "rt") as f:
            scalefactors = json.load(f)
            # {"spot_diameter_fullres": 89.49502418224989, "tissue_hires_scalef": 0.17011142,
            # "fiducial_diameter_fullres": 144.56888521748058, "tissue_lowres_scalef": 0.051033426}
        # barcode, in_tissue, array_row, array_col, pxl_col_in_fullres, pxl_row_in_fullres
        positions = pd.read_csv(tissue_positions_list_path, header=None)
        positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]
        positions.index = positions["barcode"]
        positions = positions.reindex(adata.obs.index) # 起到取子集的效果
        spatial_coords = positions[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
        adata.obsm["tissue_hires"] = spatial_coords * scalefactors["tissue_hires_scalef"]
        adata.uns["images"] = [
            dict(
                type="image",
                name="tissue_hires",
                image=tissue_hires_image_path,
                spot_diameter=scalefactors["spot_diameter_fullres"]
                * scalefactors["tissue_hires_scalef"],
            )
        ]
        return True
    else:
        return False

def __add_generic_spatial(adata, spatial_directory):
    # positions.image.csv with barcode, x, y
    # image.png
    # diameter.image.txt (optional)
    image_extensions = {".png", ".jpeg", ".jpg"}
    meta_image_extensions = {".svg"}
    found = False
    for f in os.listdir(spatial_directory):
        name, ext = os.path.splitext(f)
        ext = ext.lower()
        if ext in image_extensions:
            positions_path = os.path.join(spatial_directory, "positions." + name + ".csv")
            if os.path.exists(positions_path):
                diameter_path = os.path.join(spatial_directory, "diameter." + name + ".txt")
                spot_diameter = None
                if os.path.exists(diameter_path):
                    with open(diameter_path, "rt") as diameter_in:
                        spot_diameter = float(diameter_in.readline().strip())

                positions = pd.read_csv(positions_path, index_col="barcode")
                if not pd.api.types.is_string_dtype(positions.index):
                    positions.index = positions.index.astype(str)
                positions = positions.reindex(adata.obs.index)
                adata.obsm[name] = positions[["x", "y"]].values
                images = adata.uns.get("images", [])
                images.append(
                    dict(
                        type="image",
                        name=name,
                        image=os.path.join(spatial_directory, f),
                        spot_diameter=spot_diameter,
                    )
                )
                adata.uns["images"] = images
                found = True
        elif ext in meta_image_extensions:
            svg_path = os.path.join(spatial_directory, f)
            if os.path.exists(svg_path):
                svg_path = os.path.abspath(svg_path)
            import json
            import xml.etree.ElementTree as ET

            tree = ET.parse(svg_path)
            attrs = tree.getroot().attrib
            if "data-group" in attrs and "data-selection" in attrs:
                found = True
                images = adata.uns.get("meta_images", [])
                selection = attrs["data-selection"].replace("'", '"')
                selection = json.loads(selection)

                images.append(
                    dict(
                        type="meta_image",
                        name=name,
                        image=svg_path,
                        attrs=dict(group=attrs["data-group"], selection=selection),
                    )
                )
                adata.uns["meta_images"] = images
    return found

def add_spatial(adata, spatial_directory):
    if not __add_visium(adata, spatial_directory):
        return __add_generic_spatial(adata, spatial_directory)
    else:
        return True

def add_registerd_spatial(adata, scale_factors_path,tissue_hires_image_path):

    import json
    try:
        with open(scale_factors_path, "rt") as f:
            scalefactors = json.load(f)
            # {"spot_diameter_fullres": 89.49502418224989, "tissue_hires_scalef": 0.17011142,
            # "fiducial_diameter_fullres": 144.56888521748058, "tissue_lowres_scalef": 0.051033426}
        adata.obsm['tissue_hires'] =  adata.obsm['spatial'] * scalefactors["tissue_hires_scalef"]
        adata.uns["images"] = [
            dict(
                type="image",
                name="tissue_hires",
                image=tissue_hires_image_path,
                spot_diameter=scalefactors["spot_diameter_fullres"] * scalefactors["tissue_hires_scalef"],
            )
        ]
    except:
        return False
    else:
        return True

def get_scanpy_marker_keys(dataset):
    marker_keys = []
    try:
        dataset.uns  # dataset can be AnnData or zarr group
    except AttributeError:
        return marker_keys
    from collections.abc import Mapping

    for key in dataset.uns.keys():
        rank_genes_groups = dataset.uns[key]
        if (
            isinstance(rank_genes_groups, Mapping)
            and "names" in rank_genes_groups
            and (
                "pvals" in rank_genes_groups
                or "pvals_adj" in rank_genes_groups
                or "scores" in rank_genes_groups
            )
            and len(rank_genes_groups["names"][0]) > 0
            and not isinstance(rank_genes_groups["names"][0][0], bytes)
        ):
            marker_keys.append(key)
    return marker_keys

def dataset_schema(dataset,n_features=10):
    """Gets dataset schema.
    Returns
        schema dict. Example:
        {"version":"1.0.0",
        "categoryOrder":{
            "louvain":["0","1","2","3","4","5","6","7"],
            "leiden":["0","1","2","3","4","5","6","7"]},
        "var":["TNFRSF4","CPSF3L","ATAD3C"],
        "obs":["percent_mito","n_counts"],
        "obsCat":["louvain","leiden"],
        "shape":[2638,1838],
        "embeddings":[{"name":"X_pca","dimensions":3},{"name":"X_pca","dimensions":2},{"name":"X_umap","dimensions":2}]
    }
    """
    obs_cat = []
    obs = []
    marker_results = []

    de_results = []  # array of dicts containing params logfoldchanges, pvals_adj, scores, names
    category_to_order = {}
    embeddings = []
    field_to_value_to_color = dict()  # field -> value -> color

    # 区分obs/obsm
    for key in dataset.obs.keys():

        val = dataset.obs[key]

        if (
            pd.api.types.is_categorical_dtype(val)
            or pd.api.types.is_bool_dtype(val)
            or pd.api.types.is_object_dtype(val)
        ):
            obs_cat.append(key)
        else:
            obs.append(key)
        if pd.api.types.is_categorical_dtype(val):
            categories = val.cat.categories
            if len(categories) < 100:  # preserve order
                category_to_order[key] = dataset.obs[key].cat.categories

    # image信息
    images_node = dataset.uns.get("images", [])  # list of {type:image or meta_image, name:image name, image:path to image, spot_diameter:Number}
    image_names = list(map(lambda x: x["name"], images_node))
    layers = []
    try:
        dataset.layers  # dataset can be AnnData or zarr group
        layers = list(dataset.layers.keys())
    except AttributeError:
        pass

    for key in dataset.obsm.keys():
        dim = dataset.obsm[key].shape[1]
        if 1 < dim <= 3:
            embedding = dict(name=key, dimensions=dim)
            if dim == 2:
                try:
                    image_index = image_names.index(key)
                    embedding["spatial"] = images_node[image_index]
                except ValueError:
                    pass

            embeddings.append(embedding)

    # Marker genes
    for scanpy_marker_key in get_scanpy_marker_keys(dataset):
        rank_genes_groups = dataset.uns[scanpy_marker_key]
        has_fc = "logfoldchanges" in rank_genes_groups
        min_fold_change = 1
        params = rank_genes_groups["params"]
        if not isinstance(params, dict):
            from anndata._io.zarr import read_attribute

            params = {k: read_attribute(params[k]) for k in params.keys()}

        # pts and pts_rest in later scanpy versions
        rank_genes_groups_keys = list(rank_genes_groups.keys())
        for k in ["params", "names"]:
            if k in rank_genes_groups_keys:
                rank_genes_groups_keys.remove(k)
        if "pvals" in rank_genes_groups_keys and "pvals_adj" in rank_genes_groups_keys:
            rank_genes_groups_keys.remove("pvals")
        category = "{} ({})".format(params["groupby"], scanpy_marker_key)
        de_result_name = category
        group_names = rank_genes_groups["names"].dtype.names

        for group_name in group_names:
            group_df = pd.DataFrame(index=rank_genes_groups["names"][group_name][...])
            group_df = group_df[group_df.index != "nan"]
            for rank_genes_groups_key in rank_genes_groups_keys:
                values = rank_genes_groups[rank_genes_groups_key][group_name][...]
                column_name = "{}:{}".format(group_name, rank_genes_groups_key)
                group_df[column_name] = values

            if n_features > 0:
                markers_df = group_df
                if has_fc:
                    markers_df = group_df[
                        group_df["{}:logfoldchanges".format(group_name)] > min_fold_change
                    ]
                if len(markers_df) > 0:
                    marker_results.append(
                        dict(
                            category=de_result_name,
                            name=str(group_name),
                            features=markers_df.index[:n_features],
                        )
                    )

    # Summary
    schema_dict = {"version": "1.0.0"}
    schema_dict["results"] = de_results
    schema_dict["colors"] = field_to_value_to_color
    schema_dict["markers"] = marker_results
    schema_dict["embeddings"] = embeddings
    schema_dict["categoryOrder"] = category_to_order
    schema_dict["layers"] = layers

    var_df = dataset.var
    if not isinstance(var_df, pd.DataFrame):
        from anndata._io.zarr import read_attribute

        var_df = read_attribute(dataset.var)
    var_df.index.name = "id"
    schema_dict["var"] = var_df.reset_index().to_dict(orient="records")

    schema_dict["obs"] = obs
    schema_dict["obsCat"] = obs_cat
    shape = dataset.shape if isinstance(dataset, anndata.AnnData) else dataset.X.attrs.shape
    schema_dict["shape"] = shape
    return schema_dict

LINE_END = "\n".encode("UTF-8")

def write_jsonl(d, f, name, index, compress=False):
    output = {}
    output[name] = d
    c = ujson.dumps(output, double_precision=2, orient="values").encode("UTF-8")
    if compress:
        c = gzip.compress(c)
    start = f.tell()
    end = start + len(c)
    index[name] = [start, end - 1]
    f.write(c)
    f.write(LINE_END)

def save_adata_X(adata, f, index, compress):
    adata_X = adata.X
    names = adata.var.index
    is_sparse = scipy.sparse.issparse(adata_X)
    if is_sparse and scipy.sparse.isspmatrix_csr(adata_X):
        adata_X = adata_X.tocsc()

    for j in range(adata_X.shape[1]):
        X = adata_X[:, j]
        if is_sparse:
            X = X.toarray().flatten()
            indices = np.where(X != 0)[0]
            values = X[indices]
            write_jsonl(dict(index=indices, value=values), f, names[j], index, compress)
        else:
            write_jsonl(dict(value=X), f, names[j], index, compress)
        if j > 0 and (j + 1) % 1000 == 0:
            logger.info("Wrote adata X {}/{}".format(j + 1, adata_X.shape[1]))


def save_data_obsm(adata, f, index, compress):
    logger.info("writing adata obsm")

    for name in adata.obsm.keys():
        value = adata.obsm[name]
        dim = value.shape[1]
        d = {}
        for i in range(dim):
            d[name + "_" + str(i + 1)] = value[:, i]
        write_jsonl(d, f, name, index, compress)


def save_data_obs(adata, f, index, compress):
    logger.info("writing adata obs")
    for name in adata.obs:
        series = adata.obs[name]
        value = series
        if pd.api.types.is_categorical_dtype(series):
            value = dict(values=series.values.codes, categories=series.cat.categories.values)
        write_jsonl(value, f, name, index, compress)
    write_jsonl(adata.obs.index.values, f, "index", index, compress)

def save_dataset_jsonl(dataset, schema, output_dir, base_name, filesystem):
    compress = False
    index = {}  # key to byte start-end
    filesystem.makedirs(output_dir, exist_ok=True)
    jsonl_path = os.path.join(output_dir, base_name)
    with filesystem.open(jsonl_path, "wb") as f:
        save_adata_X(dataset, f, index, compress)
        save_data_obs(dataset, f, index, compress)
        save_data_obsm(dataset, f, index, compress)
        write_jsonl(schema, f, "schema", index)

    with filesystem.open(
        os.path.join(output_dir, base_name + ".idx.json"), "wt"
    ) as f:  # save index
        # json.dump(result, f)
        result = dict(index=index, file=os.path.basename(jsonl_path))
        f.write(ujson.dumps(result, double_precision=2, orient="values"))
