clear; clc;

save_fold = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures';

dbscan_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\DBSCAN.mat';
load(dbscan_res_file);

cl_id_sel = 4;
cl_id_sel = 3;
cl_id_sel = 8;
cl_id_sel = 6;
cl_id_sel = 7;

load(fullfile(save_fold, strcat('PCA_DBSCAN_refine_GMM_cluster_', num2str(cl_id_sel), '.mat')), ...
    'phen_data', 'cl_id_sel', 'pca_comp_sel', 'sel_vect', 'f_vals', 'pcaval', ...
    'score', 'latent', 'explained', 'explained_cum', 'f_N', 'p_N', 'pca_sel_n', ...
    'score_sel', 'gmfit', 'ref_idx', 'b_IDX');

bLab{1} = {'001', '002', '005', '016'};
bLab{2} = {'001', '003', '004', '009', '010', '012'};
bLab{3} = {'007', '012', '013', '021', '022', '024'};

bLoc{1} = 'CEP1 Structures/Control CEP1 Structures/Ctrl_CEP1-EGFP_';
bLoc{2} = 'CEP1 Structures/Septin7-Knockdown CEP1 Structures/Sept7-KD_CEP1-EGFP_';
bLoc{3} = 'CEP1 Structures/100uM FCF CEP1 Structures/100uM FCF_CEP1-EGFP_';

exp_data = zeros(size(cl_idx));

cnt = 0;
cnt2 = 0;
cl_num = max(unique(cl_idx_sel));
for uu = 1:3
    ulab = bLab{uu}; 
    for w = 1:length(ulab)
        disp([uu,w]);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        ndat = size(features1,1);
        sel_indx = filt_vect(cnt+1:cnt+ndat);
        ndat2 = sum(sel_indx);
        exp_data(cnt2+1:cnt2+ndat2) = w;
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end

cl_vect = (cl_idx == cl_id_sel);
cl_exp_data = exp_data(cl_vect);
cl_phen_data = phen_data(cl_vect);

exp_n = length(unique(cl_exp_data));

f1_mean = zeros(3,1);
f1_CI_1 = zeros(3,1);
f1_CI_2 = zeros(3,1);
f2_mean = zeros(3,1);
f2_CI_1 = zeros(3,1);
f2_CI_2 = zeros(3,1);

f1_frac = zeros(3,1);
f2_frac = zeros(3,1);

for uu = 1:3
    ref_idx_i = ref_idx(cl_phen_data == uu);
    f1_frac(uu) = sum(ref_idx_i == 1)/length(ref_idx_i);
    f2_frac(uu) = sum(ref_idx_i == 2)/length(ref_idx_i);
end

x_tick_labels = {'CEP1 (Control)', ...
    'CEP1 (Septin7 K.D.)', ...
    'CEP1 (100uM FCF)'};

fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
br = bar(1:3,[f1_frac'; f2_frac'], 'stacked');
xticks(1:3);
xticklabels(x_tick_labels);
br(1).FaceColor = 'flat';
br(2).FaceColor = 'flat';
for uu = 1:3
    br(1).CData(uu,:) = [1 0 0];
    br(2).CData(uu,:) = [0 1 0];
end
saveas(fig, fullfile(save_fold, strcat('PCA_DBSCAN_refine_GMM_cluster_', num2str(cl_id_sel), '.png')));

