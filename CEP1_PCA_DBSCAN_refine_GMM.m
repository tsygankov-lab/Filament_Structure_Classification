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

cl_id_sel = 4;
pca_comp_sel = 3;

cl_id_sel = 3;
pca_comp_sel = 1;

cl_id_sel = 8;
pca_comp_sel = 1;

cl_id_sel = 6;
pca_comp_sel = 4;

cl_id_sel = 7;
pca_comp_sel = 1;

cl_id = cl_id_sel;

sel_vect = (cl_idx_sel == cl_id);

f_vals = features1_comb_norm(sel_vect, :);

[pcaval, score, latent, ~, explained] = pca(f_vals, 'Algorithm','svd');

explained_cum = cumsum(explained);

f_N = size(features1_comb_norm, 2);
p_N = size(features1_comb_norm, 1);

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
saveas(fig, fullfile(save_fold, strcat('refine_GMM_cl_', num2str(cl_id_sel), '_PCA_explained_variance.png')));


c_N = 1000;
cm = jet(c_N);
phen_data_sel = phen_data(sel_vect);
c_data = zeros(length(phen_data_sel), 3);
for i = 1:length(phen_data_sel)
    id = fix(phen_data_sel(i)/max(phen_data_sel)*(c_N-1)) + 1;
    c_data(i,:) = cm(id,:);
end

phen_ids_1 = find(phen_data_sel == 1);
phen_ids_2 = find(phen_data_sel == 2);
phen_ids_3 = find(phen_data_sel == 3);

x_tick_labels = {'CEP1 (Control)', ...
    'CEP1 (Septin7 K.D.)', ...
    'CEP1 (100uM FCF)'};

fig = figure('Position', [50 50 1200 900]);
subplot(4,4,1);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,2), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC2');
title('PC1 and PC2');

subplot(4,4,2);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,3), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC3');
title('PC1 and PC3');

subplot(4,4,3);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,4), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC4');
title('PC1 and PC4');

subplot(4,4,4);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,5), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC5');
title('PC1 and PC5');

subplot(4,4,5);
hold on;
grid on;
box on;
scatter(score(:,2),  score(:,6), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC6');
title('PC1 and PC6');

subplot(4,4,6);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,7), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC7');
title('PC1 and PC7');


subplot(4,4,7);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,8), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC8');
title('PC1 and PC8');

subplot(4,4,8);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,9), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC9');
title('PC1 and PC9');

subplot(4,4,9);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,10), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC10');
title('PC1 and PC10');

subplot(4,4,10);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,11), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC11');
title('PC1 and PC11');

subplot(4,4,11);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,12), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC12');
title('PC1 and PC12');

subplot(4,4,12);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,13), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC13');
title('PC1 and PC13');

subplot(4,4,13);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,14), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC14');
title('PC1 and PC14');

subplot(4,4,14);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,15), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC15');
title('PC1 and PC15');

subplot(4,4,15);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,16), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
xlabel('PC1');
ylabel('PC16');
title('PC1 and PC16');

subplot(4,4,16);
hold on;
grid on;
box on;
axis off;

sc1 = scatter(score(phen_ids_1(1),1),  score(phen_ids_1(1),17), 'Marker', '.', ...
    'SizeData', 50, 'CData', c_data(phen_ids_1(1),:));
sc2 = scatter(score(phen_ids_2(1),1),  score(phen_ids_2(1),17), 'Marker', '.', ...
    'SizeData', 50, 'CData', c_data(phen_ids_2(1),:));
sc3 = scatter(score(phen_ids_3(1),1),  score(phen_ids_3(1),17), 'Marker', '.', ...
    'SizeData', 50, 'CData', c_data(phen_ids_3(1),:));
legend([sc1, sc2, sc3], x_tick_labels, ...
    'Position', get(gca, 'Position'), 'FontSize', 12);
drawnow;

saveas(fig, fullfile(save_fold, strcat('refine_GMM_cl_', num2str(cl_id_sel), '_PCA_components.png')));

score_sel = score(:,pca_comp_sel);


gmfit = fitgmdist(score_sel,2);
ref_idx = cluster(gmfit,score_sel);

if cl_id_sel == 8
    thr = 2;
    ref_idx = zeros(size(score_sel));
    ref_idx(score_sel <= thr) = 1;
    ref_idx(score_sel > thr) = 2;
end

if cl_id_sel == 6
    thr = -2.5;
    ref_idx = zeros(size(score_sel));
    ref_idx(score_sel >= thr) = 1;
    ref_idx(score_sel < thr) = 2;
end

if cl_id_sel == 7
    thr = 3;
    ref_idx = zeros(size(score_sel));
    ref_idx(score_sel <= thr) = 1;
    ref_idx(score_sel > thr) = 2;
end

if cl_id_sel == 4
    ref_idx(ref_idx == 1) = 3;
    ref_idx(ref_idx == 2) = 1;
    ref_idx(ref_idx == 3) = 2;
end

fig = figure('Position', [50 50 500 500]); 
hold on;
grid on;
box on;
histogram(score_sel, 'Normalization','pdf');
%plot([thr, thr], get(gca, 'YLim'), 'Color', [1 0 0], 'LineWidth', 2);
title(strcat('PC', num2str(pca_comp_sel)));
drawnow;
saveas(fig, fullfile(save_fold, strcat('refine_GMM_cl_', num2str(cl_id_sel), '_PC_', num2str(pca_comp_sel), '.png')));

b_IDX = cell(1,3);

cnt = 0;
cnt2 = 0;
cnt3 = 0;
cl_num = max(unique(cl_idx_sel));
for uu = 1:3
    ulab = bLab{uu};
    for w = 1:length(ulab)
        disp([uu,w]);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        ndat = size(features1,1);
        
        im_obj_ids = unique(sort(Lexcl(Lexcl>0)));
        filt_indx = filt_vect(cnt+1:cnt+ndat);
        im_obj_filt_ids = im_obj_ids(filt_indx);
        ndat2 = sum(filt_indx);
        sel_indx = sel_vect(cnt2+1:cnt2+ndat2);
        im_obj_sel_ids = im_obj_filt_ids(sel_indx);
        ndat3 = sum(sel_indx);
        indx = ref_idx(cnt3+1:cnt3+ndat3);
        
        b_IDX{uu}{w} = {};
        for cl_id = 1:2
            im_indx = ismember(Lexcl(:), im_obj_sel_ids(indx == cl_id));
            b_IDX{uu}{w}{cl_id} = im_indx;
        end
        
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
        cnt3 = cnt3 + ndat3;
    end
end

cm = [1, 0, 0; 0, 1, 0];

fig = figure('Position',get(0,'Screensize'));
for uu = 1:3
%for uu = 3:3
    ulab = bLab{uu}; 
    for w = 1:length(ulab)
        disp([uu,w]);
        loc = [bLoc{uu} ulab{w}];
        mkdir([loc '_PCA_DBSCAN_refine_GMM']);
        load([loc '.mat']);
        IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
        for cl_id = 1:2
            im_indx = b_IDX{uu}{w}{cl_id};
            IM_r = zeros(size(Lexcl));
            IM_g = zeros(size(Lexcl));
            IM_b = zeros(size(Lexcl));
            IM_r(im_indx) = cm(cl_id,1);
            IM_g(im_indx) = cm(cl_id,2);
            IM_b(im_indx) = cm(cl_id,3);
            %IM_rgb = cat(3, IM_r, IM_g, IM_b);
            IM_rgb(:,:,1) = IM_rgb(:,:,1) + IM_r;
            IM_rgb(:,:,2) = IM_rgb(:,:,2) + IM_g;
            IM_rgb(:,:,3) = IM_rgb(:,:,3) + IM_b;
        end
        clf;
        hold on;
        axis ij;
        axis image;
        axis off;
        image(IM_rgb);
        title(strcat('DBSCAN cluster=', num2str(cl_id_sel), ' tail threshold(PC=', num2str(pca_comp_sel), ')'));
        drawnow;
        saveas(fig, [loc '_PCA_DBSCAN_refine_GMM' ...
            strcat('/cluster_', num2str(cl_id_sel), '_tail_threshold(PC=', num2str(pca_comp_sel), ').png')]);
    end
end

save(fullfile(save_fold, strcat('PCA_DBSCAN_refine_GMM_cluster_', num2str(cl_id_sel), '.mat')), ...
    'phen_data', 'cl_id_sel', 'pca_comp_sel', 'sel_vect', 'f_vals', 'pcaval', ...
    'score', 'latent', 'explained', 'explained_cum', 'f_N', 'p_N', 'pca_sel_n', ...
    'score_sel', 'gmfit', 'ref_idx', 'b_IDX');


