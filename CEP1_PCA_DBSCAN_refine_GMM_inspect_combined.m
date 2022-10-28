clear; clc;

save_fold = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures';

dbscan_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\DBSCAN.mat';
load(dbscan_res_file);

bLab{1} = {'001', '002', '005', '016'};
bLab{2} = {'001', '003', '004', '009', '010', '012'};
bLab{3} = {'007', '012', '013', '021', '022', '024'};

bLoc{1} = 'CEP1 Structures/Control CEP1 Structures/Ctrl_CEP1-EGFP_';
bLoc{2} = 'CEP1 Structures/Septin7-Knockdown CEP1 Structures/Sept7-KD_CEP1-EGFP_';
bLoc{3} = 'CEP1 Structures/100uM FCF CEP1 Structures/100uM FCF_CEP1-EGFP_';

phen_data = zeros(size(cl_idx));

cnt = 0;
cnt2 = 0;
for uu = 1:3
    ulab = bLab{uu}; 
    for w = 1:length(ulab)
        disp([uu,w]);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        ndat = size(features1,1);
        sel_indx = filt_vect(cnt+1:cnt+ndat);
        ndat2 = sum(sel_indx);
        
        phen_data(cnt2+1:cnt2+ndat2) = uu;
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end

cl_idx_comb = cl_idx_sel;

ref_cl_4_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\PCA_DBSCAN_refine_GMM_cluster_4.mat';
ref_cl_3_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\PCA_DBSCAN_refine_GMM_cluster_3.mat';
ref_cl_8_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\PCA_DBSCAN_refine_GMM_cluster_8.mat';
ref_cl_6_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\PCA_DBSCAN_refine_GMM_cluster_6.mat';
ref_cl_7_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\PCA_DBSCAN_refine_GMM_cluster_7.mat';

load(ref_cl_4_res_file);
cl_n = max(unique(sort(cl_idx_comb)));
cnt = 1;
for i = 1:length(cl_idx_comb)
    if cl_idx_comb(i) == cl_id_sel
        if ref_idx(cnt) == 2
            cl_idx_comb(i) = cl_n+1;
        end
        cnt = cnt+1;
    end
end

load(ref_cl_3_res_file);
cl_n = max(unique(sort(cl_idx_comb)));
cnt = 1;
for i = 1:length(cl_idx_comb)
    if cl_idx_comb(i) == cl_id_sel
        if ref_idx(cnt) == 2
            cl_idx_comb(i) = cl_n+1;
        end
        cnt = cnt+1;
    end
end

load(ref_cl_8_res_file);
cl_n = max(unique(sort(cl_idx_comb)));
cnt = 1;
for i = 1:length(cl_idx_comb)
    if cl_idx_comb(i) == cl_id_sel
        if ref_idx(cnt) == 2
            cl_idx_comb(i) = cl_n+1;
        end
        cnt = cnt+1;
    end
end

load(ref_cl_6_res_file);
cl_n = max(unique(sort(cl_idx_comb)));
cnt = 1;
for i = 1:length(cl_idx_comb)
    if cl_idx_comb(i) == cl_id_sel
        if ref_idx(cnt) == 2
            cl_idx_comb(i) = cl_n+1;
        end
        cnt = cnt+1;
    end
end

load(ref_cl_7_res_file);
cl_n = max(unique(sort(cl_idx_comb)));
cnt = 1;
for i = 1:length(cl_idx_comb)
    if cl_idx_comb(i) == cl_id_sel
        if ref_idx(cnt) == 2
            cl_idx_comb(i) = cl_n+1;
        end
        cnt = cnt+1;
    end
end

cell_data = zeros(size(cl_idx_comb));

cnt = 0;
cnt2 = 0;
cell_n = 1;
for uu = 1:3
    ulab = bLab{uu}; 
    for w = 1:length(ulab)
        disp([uu,w]);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        ndat = size(features1,1);
        sel_indx = filt_vect(cnt+1:cnt+ndat);
        ndat2 = sum(sel_indx);
        cell_data(cnt2+1:cnt2+ndat2) = cell_n;
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
        cell_n = cell_n + 1;
    end
end


cl_num = max(cl_idx_comb);
cell_num = max(cell_data);

dens_data = zeros(cell_num, cl_num);
leg_plts = [];
cm = turbo(3);

fig = figure('Position', [50 50 1500 700]);
hold on;
grid on;
box on;
for uu = 1:3
    ulab = bLab{uu}; 
    dens_vals = cell(1, cl_num);
    dens_mean = zeros(1, cl_num);
    dens_CI_1 = zeros(1, cl_num);
    dens_CI_2 = zeros(1, cl_num);
    
    for cl_id = 1:cl_num
        for cell_i = unique(cell_data(phen_data == uu))'
            struct_i_idx = (phen_data == uu) & (cl_idx_comb == cl_id) & (cell_data == cell_i);
            struct_all_idx = (phen_data == uu) & (cell_data == cell_i);
            dens_vals{cl_id} = [dens_vals{cl_id}, sum(struct_i_idx)/sum(struct_all_idx)];
        end
        
        dens_mean(cl_id) = mean(dens_vals{cl_id});
        ts = tinv([0.025  0.975], length(dens_vals{cl_id})-1);
        SEM = std(dens_vals{cl_id})/sqrt(length(dens_vals{cl_id}));
        dens_CI_1(cl_id) = dens_mean(cl_id) + ts(1)*SEM;
        dens_CI_2(cl_id) = dens_mean(cl_id) + ts(2)*SEM;
        
        eb = errorbar(cl_id+uu/15, dens_mean(cl_id), ...
        dens_CI_1(cl_id)-dens_mean(cl_id), dens_CI_2(cl_id)-dens_mean(cl_id), ...
        'Marker', '.', 'MarkerSize', 20, 'LineStyle', 'none', 'LineWidth', 1, 'Color', cm(uu,:));
        leg_plts(uu) = eb;
        drawnow;
        
    end
end
xticks([1:cl_num]+1/5);
xticklabels([1:cl_num]);
leg_names = {'CEP1 Structures Control', ...
    'CEP1 Structures Septin7-Knockdown', ...
    'CEP1 Structures 100uM FCF'};
legend(leg_plts, leg_names);
xlabel('Cluster id');
ylabel('Structure density');
set(gca, 'FontSize', 12);

saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_refine_GMM_comined_statistics.png'));

ylim([0, 0.2]);
saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_refine_GMM_comined_statistics2.png'));


cl_num = max(cl_idx_comb);
cm = jet(cl_num);

cnt = 0;
cnt2 = 0;
fig = figure('Position',get(0,'Screensize'));
for uu = 1:3
    ulab = bLab{uu};
    cell_ids = unique(cell_data(phen_data == uu));
    for w = 1:length(cell_ids)
        cell_i = cell_ids(w);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        mkdir([loc '_PCA_DBSCAN_refine_GMM_combined_plot']);
        ndat = size(features1,1);
        im_obj_ids = unique(sort(Lexcl(Lexcl>0)));
        sel_indx = filt_vect(cnt+1:cnt+ndat);
        im_obj_sel_ids = im_obj_ids(sel_indx);
        ndat2 = sum(sel_indx);
        struct_all_idx = (phen_data == uu) & (cell_data == cell_i);
        cl_cell_idx_comb = cl_idx_comb(struct_all_idx);
        IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
        for cl_id_i = 1:cl_num
            im_indx = ismember(Lexcl(:), im_obj_sel_ids(cl_cell_idx_comb == cl_id_i));
            
            IM_r = zeros(size(Lexcl));
            IM_g = zeros(size(Lexcl));
            IM_b = zeros(size(Lexcl));
            IM_r(im_indx) = cm(cl_id_i,1);
            IM_g(im_indx) = cm(cl_id_i,2);
            IM_b(im_indx) = cm(cl_id_i,3);
            %IM_rgb = cat(3, IM_r, IM_g, IM_b);
            IM_rgb(:,:,1) = IM_rgb(:,:,1) + IM_r;
            IM_rgb(:,:,2) = IM_rgb(:,:,2) + IM_g;
            IM_rgb(:,:,3) = IM_rgb(:,:,3) + IM_b;
            
            IM_rgb_i = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
            IM_rgb_i(:,:,1) = IM_r;
            IM_rgb_i(:,:,2) = IM_g;
            IM_rgb_i(:,:,3) = IM_b;
            %IM = zeros(size(Lexcl));
            %IM(im_indx) = 1;
            
            clf;
            hold on;
            axis ij;
            axis image;
            axis off;
            image(IM_rgb_i);
            title(strcat('cluster=', num2str(cl_id_i)));
            drawnow;
            saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_combined_plot', '\cluster=' num2str(cl_id_i) '.png']);
        end
        clf;
        hold on;
        axis ij;
        axis image;
        axis off;
        image(IM_rgb);
        title('all clusters');
        drawnow;
        saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_combined_plot', '\all clusters.png']);
        
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end



dens_data = zeros(cell_num, cl_num);
dens_phen = zeros(cell_num, 1);


cell_counter = 1;
for uu = 1:3
    ulab = bLab{uu};
    cell_ids = unique(cell_data(phen_data == uu));
    for w = 1:length(cell_ids)
        cell_i = cell_ids(w);
        ndat = size(features1,1);
        for cl_id = 1:cl_num
            struct_i_idx = (phen_data == uu) & (cl_idx_comb == cl_id) & (cell_data == cell_i);
            struct_all_idx = (phen_data == uu) & (cell_data == cell_i);
            dens_data(cell_counter, cl_id) = sum(struct_i_idx)/sum(struct_all_idx);
        end
        dens_phen(cell_counter) = uu;
        cell_counter = cell_counter + 1;
    end
end



f_vals = dens_data;
f_vals = (dens_data - mean(dens_data,1))./std(dens_data,1);

[pcaval, score, latent, ~, explained] = pca(f_vals, 'Algorithm','svd');

explained_cum = cumsum(explained);

f_N = size(f_vals, 2);
p_N = size(f_vals, 1);

pca_sel_n = sum(explained_cum  < 90);

fig = figure('Position', [50 50 1000 500]);
hold on;
grid on;
box on;
plot(1:f_N, explained_cum, 'Color', [1 0 0], 'LineWidth', 3);
plot([1, f_N], [90 90], 'Color', [0 0 0], 'LineWidth', 1);
xlim([1 f_N]);
xlabel('Principal component');
ylabel('Cumulative variance explained %');
title({strcat('number of PCs that explain more then 90% of variance:', num2str(pca_sel_n)), ...
    '(will be used for clustering with DBSCAN)'});
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_refine_GMM_comined_cells_structure_density_pca_variance_explained.png'));


c_N = 1000;
cm = turbo(c_N);
c_data = zeros(length(dens_phen), 3);
for i = 1:length(dens_phen)
    id = fix(dens_phen(i)/max(dens_phen)*(c_N-1)) + 1;
    c_data(i,:) = cm(id,:);
end

phen_ids_1 = find(dens_phen == 1);
phen_ids_2 = find(dens_phen == 2);
phen_ids_3 = find(dens_phen == 3);

sz = 500;

fig = figure('Position', [50 50 800 800]);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data, 'LineWidth', 2);
sc1 = scatter(score(phen_ids_1(1),1),  score(phen_ids_1(1),2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data(phen_ids_1(1),:), 'LineWidth', 2);
sc2 = scatter(score(phen_ids_2(1),1),  score(phen_ids_2(1),2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data(phen_ids_2(1),:), 'LineWidth', 2);
sc3 = scatter(score(phen_ids_3(1),1),  score(phen_ids_3(1),2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data(phen_ids_3(1),:), 'LineWidth', 2);
legend([sc1, sc2, sc3], leg_names, 'Location', 'NorthWest');
xlabel('PC1');
ylabel('PC2');
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_refine_GMM_comined_cells_structure_density_pca_first_two_components.png'));


save(fullfile(save_fold, 'PCA_DBSCAN_refine_GMM_comined.mat'), ...
    'phen_data' , 'cl_idx_comb', 'cell_data', 'cl_num', 'cell_num', ...
    'dens_data', 'dens_phen');
