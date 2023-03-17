Dr. Denis Tsygankov (2023)

Integrative Systems Biology Lab

Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology and Emory University SOM

If you use any part of the scripts in this package, please cite:

Kho M, Hladyshau S, Tsygankov D, Nie S. Coordinated regulation of Cdc42ep1, actin, and septin filaments during neural crest cell migration. Front Cell Dev Biol. 2023 Feb 27;11:1106595. doi: 10.3389/fcell.2023.1106595. PMID: 36923257; PMCID: PMC10009165.

# Filament_Structure_Classification

compute_features.m - image segmentation and features extraction for Septin and CEP1 data.

fourier_shape.m - supplemental function for computation of Fourier modes

## Septin filaments analysis

Septin_PCA_DBSCAN_classify.m - dimensionality reduction based on PCA analysis and structures classification with DBSCAN

Septin_PCA_DBSCAN_inspect.m - inspect results of structures classification, sort clusters based on area

Septin_PCA_DBSCAN_plot.m - plot results of classification for all clusters

Septin_PCA_DBSCAN_plot_cropped.m - plot selected regions of images for all clusters

Septin_PCA_DBSCAN_inspect_area_density.m - image-level analysis based on relative density of structures area

Septin_PCA_DBSCAN_inspect_join_area_density.m - join clusters (to simplify interpretation) and perform image-level analysis based on relative density of structures area

Septin_PCA_DBSCAN_plot_join_cropped.m  - plot selected regions of images for joined clusters

## CEP1 filaments analysis

CEP1_PCA_DBSCAN_classify.m - dimensionality reduction based on PCA analysis and structures classification with DBSCAN

CEP1_PCA_DBSCAN_inspect.m - inspect results of structures classification, sort clusters based on area

CEP1_PCA_DBSCAN_plot.m - plot results of classification for all clusters

CEP1_PCA_DBSCAN_plot_cropped.m - plot selected regions of images for all clusters

CEP1_PCA_DBSCAN_refine_GMM.m - perform refinement of selected clusters with GMM or with the threshold for selected principal component

CEP1_PCA_DBSCAN_refine_GMM_inspect.m - inspect results of cluster refinement for different phenotypes

CEP1_PCA_DBSCAN_refine_GMM_comb_plot_cropped.m - plot results of classification for combined clusters after refinement

CEP1_PCA_DBSCAN_refine_GMM_inspect_combined.m - combine results of cluster refinement with the original clustering and perform image-level analysis based on the relative density of structures

CEP1_PCA_DBSCAN_refine_GMM_inspect_combined_join_area_density.m - join clusters (to simplify interpretation) and perform image-level analysis based on relative density of structures area

CEP1_PCA_DBSCAN_refine_GMM_comb_join_plot_cropped.m  - plot selected regions of images for joined clusters
