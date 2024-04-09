import json
import os
import torch
import sys
import pandas as pd
import numpy as np
import utils
import tangram as tg
import assign_coordinates as ac


class Tangram:
    def __init__(self, ad_sc, ad_sp, output_dir='./', cell_type_key='cell_type', use_raw_sc=None,
                 use_raw_sp=None):
        # Read in ad_sc & ad_sp
        utils.info("Read in ad_sc and ad_sp")
        # unique cell_id
        ad_sc.obs['id'] = ac.make_elements_unique(ad_sc.obs.index)
        self.adata_sc = ad_sc
        self.adata_sp = ad_sp
        markers_df = pd.DataFrame(self.adata_sc.uns["rank_genes_groups_cell_type"]["names"]).iloc[0:100, :]
        self.markers = list(np.unique(markers_df.melt().value.values))
        # copy raw gene name
        self.adata_sc.var["name"] = self.adata_sc.var.index

        if 'mt' in ad_sp.var.columns.tolist():
            self.adata_sp = ad_sp[:, ~ad_sp.var['mt']].copy()
        if use_raw_sc:
            self.adata_sc = ad_sc.raw.to_adata().copy()
        if use_raw_sp:
            self.adata_sc = ad_sp.raw.to_adata().copy()

        self.cell_type_key = cell_type_key
        self.output_dir = output_dir
        self.mode = None
        self.adata_map = None
        self.adata_sc_ann = None
        self.adata_sp_ann = None
        self.adata_ge = None

    def run_tangram(self, marker=None, key_deg="rank_genes_groups", top_n_markers=100, gpu_index=None, n_threads=30,
                    mode="cells", num_epochs=1000):
        # Get marker
        utils.info("Get marker")
        if marker is not None:
            marker = pd.read_csv(marker, header=None)
            marker = list(marker[0])
        else:
            marker = self.markers
        # else:
        #    marker = pd.DataFrame(self.adata_sc.uns[key_deg]['names']).head(top_n_markers)
        #    marker = np.array(marker).flatten().tolist()

        # Run pp_adatas
        utils.info("Run pp_adatas")
        tg.pp_adatas(adata_sc=self.adata_sc, adata_sp=self.adata_sp, genes=marker)
        if not self.adata_sc.uns['training_genes'] == self.adata_sp.uns['training_genes']:
            print("Training genes in ad_sc and ad_sp are not identical!")
            sys.exit(1)

        utils.info("Run tangram mapping")
        # Mapping
        self.mode = mode
        if gpu_index is not None and torch.cuda.is_available():
            with torch.cuda.device(gpu_index):
                if self.mode == 'cells':
                    adata_map = tg.map_cells_to_space(
                        adata_sc=self.adata_sc,
                        adata_sp=self.adata_sp,
                        device="cuda",
                        n_threads=n_threads,
                        num_epochs=num_epochs,
                        mode=self.mode
                    )
                else:
                    adata_map = tg.map_cells_to_space(
                        adata_sc=self.adata_sc,
                        adata_sp=self.adata_sp,
                        device="cuda",
                        n_threads=n_threads,
                        mode=self.mode,
                        num_epochs=num_epochs,
                        cluster_label=self.cell_type_key
                    )
            torch.cuda.empty_cache()
        else:
            if self.mode == 'cells':
                adata_map = tg.map_cells_to_space(
                    adata_sc=self.adata_sc,
                    adata_sp=self.adata_sp,
                    device="cpu",
                    # n_threads=n_threads,
                    mode=self.mode,
                    num_epochs=num_epochs
                )
            else:
                adata_map = tg.map_cells_to_space(
                    adata_sc=self.adata_sc,
                    adata_sp=self.adata_sp,
                    device="cpu",
                    n_threads=n_threads,
                    mode=self.mode,
                    num_epochs=num_epochs,
                    cluster_label=self.cell_type_key
                )
        self.adata_map = adata_map
        # transfer to raw gene name
        self.adata_sc.var.index = self.adata_sc.var["name"]

    def project_cell_positions(self):
        utils.info("Project single cells to the positions of ST data")
        pred_spots = np.argmax(self.adata_map.X, axis=1)
        self.adata_sc_ann = self.adata_sc
        self.adata_sc_ann.obs["id_st"] = self.adata_sp.obs.iloc[pred_spots].index
        self.adata_sc_ann.obs['id_new'] = self.adata_sc_ann.obs['id']
        # raw x
        self.adata_sc_ann.obs["x"] = self.adata_sp.obsm["spatial"][pred_spots, 0]
        # raw y
        self.adata_sc_ann.obs["y"] = self.adata_sp.obsm["spatial"][pred_spots, 1]
        # noise
        noise_x_Y = ac.add_coord_noise(np.array(self.adata_sc_ann.obs[["x","y"]]), 10)
        # print(np.array(self.adata_sc_ann.obs[["x","y"]]))
        # print(noise_x_Y)
        self.adata_sc_ann.obsm["spatial"] = noise_x_Y
        # x with noise
        self.adata_sc_ann.obs["x_noise"] = self.adata_sc_ann.obsm["spatial"][:, 0]
        # y with noise
        self.adata_sc_ann.obs["y_noise"] = self.adata_sc_ann.obsm["spatial"][:, 1]

        # Generate mapping_summary
        map_factor = self.adata_sc.obs.loc[:, 'cell_type'].value_counts().index
        map_summary_mapped = self.adata_sc.obs.loc[
            self.adata_sc.obs['id'].isin(self.adata_sc_ann.obs['id']), 'cell_type'].value_counts()
        map_summary_mapped = map_summary_mapped.reindex(map_factor, fill_value=0)
        map_summary_failed = self.adata_sc.obs.loc[
            ~self.adata_sc.obs['id'].isin(self.adata_sc_ann.obs['id']), 'cell_type'].value_counts()
        map_summary_failed = map_summary_failed.reindex(map_factor, fill_value=0)
        map_summary = {'Cell_type': list(map_factor),
                       'Failed': list(map_summary_failed),
                       'Mapped': list(map_summary_mapped)}
        # save json
        with open(os.path.join(self.output_dir, 'out', 'json', 'mapping_summary.json'), 'w') as f:
            json.dump(map_summary, f)
        self.adata_sc_ann.obs.to_csv(os.path.join(self.output_dir, 'out', 'table', 'sc_coordinate.csv'), index=False)
        self.adata_sc_ann.obs[['id_new', 'cell_type']].to_csv(os.path.join(self.output_dir, 'meta.txt'),
                                                              index=False, sep="\t")

    def project_cell_annotations(self):
        utils.info("Project prediction cell types to ST data")
        self.adata_sp_ann = tg.project_cell_annotations(self.adata_map, self.adata_sp, annotation=self.cell_type_key)
        self.adata_sp.obs[self.cell_type_key] = self.adata_sp.obsm["tangram_ct_pred"].idxmax(axis=1)

    def project_genes(self):
        self.adata_ge = tg.project_genes(self.adata_map, self.adata_sc, cluster_label=self.cell_type_key)

    def return_adata_sp_ann(self):
        return self.adata_sp_ann

    def save_map_X(self):
        return self.adata_map.X

    def save(self):
        utils.info("Save h5ad and table results")
        # map_X = pd.DataFrame(self.adata_map.X, index=self.adata_sc.obs.index, columns=self.adata_sp.obs.index)
        # map_X.to_csv(self.output_dir + 'sc_st_map_pro.' + self.mode + '.csv', index=True)
        if self.adata_sc_ann:
            self.adata_sc_ann.write(os.path.join(self.output_dir, 'sc_registered.h5ad'))
        if self.adata_sp_ann:
            self.adata_sp_ann.write(os.path.join(self.output_dir, 'adata_sp_ann.' + self.mode + '.h5ad'))
        if self.adata_ge:
            self.adata_ge.write(os.path.join(self.output_dir, 'adata_ge.' + self.mode + '.h5ad'))
