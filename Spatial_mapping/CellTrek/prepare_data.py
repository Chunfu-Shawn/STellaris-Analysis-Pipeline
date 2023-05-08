import os
import sys
import logging
import argparse
import anndata
import scipy.sparse

from prepare_utils import SPATIAL_HELP, add_spatial
from prepare_utils import dataset_schema
from prepare_utils import get_fs
from prepare_utils import save_dataset_jsonl

logger = logging.getLogger("prepare_data")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s:%(step)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)

def read_adata(path, spatial_directory=None, use_raw=False):
    #### Read h5ad
    logger.info("Read h5ad",extra={'step':'read_adata()'})
    adata = anndata.read(path)

    if use_raw and adata.raw is not None and adata.shape[0] == adata.raw.shape[0]:
        adata = anndata.AnnData(
            X=adata.raw.X, var=adata.raw.var, obs=adata.obs, obsm=adata.obsm, uns=adata.uns
        )
    #### Read image
    logger.info("Read Image",extra={'step':'read_adata()'})
    if spatial_directory is not None:
        if not add_spatial(adata, spatial_directory):
            print("No spatial data found in {}".format(spatial_directory))
        else:
            print("Spatial data found in {}".format(spatial_directory) + " and impoted successfully!")
    else:
        print("No spatial data was provided, skipping reading image")

    return adata

class PrepareData:
    def __init__(self,dataset,output,name,group=None,output_format="jsonl"):
        self.dataset = dataset
        self.output = output
        self.name = name
        self.group = group
        self.output_format = output_format

        logger.info("Initialize anndata",extra={'step':'PrepareData init'})
        #### Process obsm
        print("Process obsm")
        # Remove key in obsm with > 3 columns
        for key in list(dataset.obsm.keys()):
            m = dataset.obsm[key]
            dim = m.shape[1]
            if not (1 < dim <= 3):
                del dataset.obsm[key]

        #### Process var and obs
        print("Process var and obs")
        if "group" not in dataset.var:
            dataset.var["group"] = dataset.uns.get("name")

        dataset.var_names_make_unique() # make var index name unique
        dataset.var.index.name = "id" # set var index nmae as id
        dataset.obs.index.name = "id" # set obs index name as id

        #### Process .X
        print("Process X")
        if scipy.sparse.issparse(dataset.X) and not scipy.sparse.isspmatrix_csc(dataset.X):
            dataset.X = dataset.X.tocsc()
        for layer_name in dataset.layers.keys():
            X = dataset.layers[layer_name]
            if scipy.sparse.issparse(X) and not scipy.sparse.isspmatrix_csc(X):
                dataset.layers[layer_name] = X.tocsc()

    def execute(self):
        output_format = self.output_format
        dataset = self.dataset
        group = self.group

        logger.info("Convert to jsonl",extra={'step':'PrepareData execute'})

        #### Perform DEX
        if group != None:
            logger.info("Compute marker genes", extra={'step': 'PrepareData execute'})
            import scanpy as sc

            key_added = "rank_genes_" + str(group)
            value_counts = dataset.obs[group].value_counts()
            filtered_value_counts = value_counts[value_counts >= 3]
            if len(filtered_value_counts) >= 2:
                print("Computing markers for {}".format(group))
                sc.tl.rank_genes_groups(
                    dataset,
                    group,
                    key_added=key_added,
                    method="t-test",
                    groups=filtered_value_counts.index.to_list(),
                    use_raw=False
                )
        else:
            logger.info("Skip computing marker genes", extra={'step': 'PrepareData execute'})

        #### Get schema
        print("Get schema")
        schema = self.get_schema()

        #### Create output directory
        print("Create output directory")
        output_dir = os.path.dirname(self.output)
        base_output = os.path.basename(self.output)
        filesystem = get_fs(output_dir) # Access system platform, default is linux
        filesystem.makedirs(output_dir, exist_ok=True)

        #### Images
        print("Save image")
        images = dataset.uns.pop("images", None)
        if images is not None:
            image_dir = os.path.join(output_dir, "images")
            filesystem.makedirs(image_dir, exist_ok=True)
            for image in images:
                src = image["image"]
                dest = os.path.join(image_dir, os.path.basename(src))
                filesystem.copy(src, dest)
                image["image"] = "images/" + os.path.basename(src)

        #### Save jsonl
        print("Save jsonl")
        if output_format == "jsonl":
            save_dataset_jsonl(dataset, schema, output_dir, base_output, filesystem)
        else:
            raise ValueError("Unknown format")

    def get_schema(self):
        result = dataset_schema(self.dataset,n_features=10)
        result["format"] = self.output_format
        return result

def main(argsv):
    parser = argparse.ArgumentParser(description="Prepare a dataset for cirrocumulus server.")
    parser.add_argument("--dataset", metavar="<H5AD>", help="Path to a h5ad file (Only support one h5ad for now)")
    parser.add_argument("--outDir", help="Path to output directory")
    parser.add_argument("--name", help="Output name or ID (e.g., duplicate ID)")
    parser.add_argument("--group", help="Group to compute markers if specified",default=None)
    parser.add_argument("--format", help="Output format", choices=["jsonl"], default="jsonl")
    parser.add_argument("--spatial", help=SPATIAL_HELP)

    args = parser.parse_args()

    #### Check input
    logger.info("Check input",extra={'step':'Main'})

    output_format2extension = dict(parquet=".cpq", jsonl=".jsonl", zarr=".zarr", h5ad=".h5ad")
    output = os.path.join(args.outDir,args.name + output_format2extension[args.format])

    print("Dataset to process: ", args.dataset)
    print("Sample ID (e.g., section ID) to save: ", args.name)
    print("Output directory: ", args.outDir)
    print("Jsonl to save: ", output)
    print("Images provided: ", args.spatial != None)
    print("Compute marker genes: ", False if args.group == None else 'True, Using {} field'.format(args.group))

    #### Read dataset
    logger.info("Read dataset",extra={'step':'Main'})
    adata = read_adata(args.dataset, spatial_directory=args.spatial, use_raw=False)
    adata.uns["name"] = args.name

    #### Prepare data
    logger.info("Prepare data",extra={'step':'Main'})
    prepare_data = PrepareData(
        dataset=adata,
        output=output,
        name=args.name,
        group=args.group,
        output_format=args.format
    )

    prepare_data.execute()

    logger.info("All finished!",extra={'step':'Main'})

if __name__ == "__main__":
    import sys

    main(sys.argv)