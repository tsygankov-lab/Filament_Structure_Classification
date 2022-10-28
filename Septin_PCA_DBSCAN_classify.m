clear; clc;

save_fold = 'E:\Septin_structure_analysis\updated_pictures\data7\Septin_structures';
mkdir(save_fold);

% reading quantitative metrics from mat files 
bLab = cell(1,3);
bLoc = cell(1,3);
bNum = zeros(1,3);

bLab{1} = {'006', '109', '110', '112'};
bLab{2} = {'001', '010', '101', '102', '103'};
bLab{3} = {'003', '101', '108', '109', '110'};

bLoc{1} = 'Septin Structures/Control Septin Structures/Ctrl_Sept7_';
bLoc{2} = 'Septin Structures/CEP1 Knockdown Septin Structures/CEP1-KD_Sept7_';
bLoc{3} = 'Septin Structures/CEP1 Overexpression Septin Structures/CEP1-OE_Sept7_';


features1_comb = [];

for uu = 1:3   % six phenotypes
    ulab = bLab{uu}; 
    cnt = 0;
    for w = 1:length(ulab)   % images from each phenotype
        disp([uu,w]);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        cnt = cnt + size(features1,1);
        features1_comb = [features1_comb; features1];
    end
    bNum(uu) = cnt;
end

disp('Number of metrics from each phenotype');
disp(bNum);


areas = features1_comb(:,1); 
min_area_th = 10;
max_area_th = 400;
filt_vect = (areas >= min_area_th) & (areas <= max_area_th);
features1_comb_filt = features1_comb(filt_vect,:);

fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
histogram(features1_comb(:,1));
title('area');
drawnow;
saveas(fig, fullfile(save_fold, 'area_histogram.png'));

fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
histogram(features1_comb_filt(:,1));
title('area (filtered)');
saveas(fig, fullfile(save_fold, 'area_histogram_filtered.png'));

% normalizing all metrics
features1_comb_norm = features1_comb_filt;
    
for i = 1:size(features1_comb_filt,2)
    
    if std(features1_comb_filt(:,i))>0
        features1_comb_norm(:,i) = (features1_comb_filt(:,i) - mean(features1_comb_filt(:,i)))/std(features1_comb_filt(:,i));
    else
        features1_comb_norm(:,i) = features1_comb_filt(:,i) - mean(features1_comb_filt(:,i));
    end
end
 
% displaying normalized metrics
% displaying combined metrics
fig = figure('Position',get(0,'Screensize'));
subplot(2,1,1);
imagesc(log(features1_comb));
title('combimed features');
subplot(2,1,2);
tmp = features1_comb_norm;
tmp(tmp>10) = 10;
imagesc(tmp);
title('normalized combined features');
 
% PCA
[pcaval,score,latent,~,explained] = pca(features1_comb_norm,'Algorithm','svd');

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
saveas(fig, fullfile(save_fold, 'PCA_explained_variance.png'));


fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
scatter3(score(:,1),  score(:,2), score(:,3), 'Marker', '.', ...
    'SizeData', 20, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_first_two_components.png'));


score_sel = score(:,1:pca_sel_n);

epsilon = 2;
minpts = 3;

GR_thr_perc = 0.6;
tic;
cl_idx = dbscan(score_sel, epsilon, minpts);
toc;


c_N = 1000;
cm = jet(c_N);
c_data = zeros(length(cl_idx), 3);
for i = 1:length(cl_idx)
    if cl_idx(i) ~= -1
        id = fix(cl_idx(i)/max(cl_idx)*(c_N-1)) + 1;
        c_data(i,:) = cm(id,:);
    else
        c_data(i,:) = [0 0 0];
    end
end

fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
scatter3(score_sel(:,1),  score_sel(:,2), score_sel(:,3), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_first_two_components_DBSCAN_all.png'));

disp('number of clusters');
disp(max(cl_idx));

disp('percent of unclassified points');
disp(sum(cl_idx==-1)/length(cl_idx)*100);


[GR_counts, GR_ids] = groupcounts(cl_idx);
GR_percents = GR_counts/sum(GR_counts)*100;

[~, GR_sort_ids] = sort(GR_counts, 'descend');

x_tick_labels = cell(1,1);
for i = 1:length(GR_ids)
    x_tick_labels{i} = num2str(GR_ids(GR_sort_ids(i)));
end

cl_top_idx = GR_ids(GR_percents > GR_thr_perc & GR_ids ~= -1);
disp(length(cl_top_idx));
cl_idx_sel = zeros(length(cl_idx), 1);
for i = 1:length(cl_top_idx)
    cl_idx_sel(cl_idx == cl_top_idx(i)) = i;
end

c_data = zeros(length(cl_idx_sel), 3);
for i = 1:length(cl_idx_sel)
    if cl_idx(i) ~= -1
        id = fix(cl_idx_sel(i)/max(cl_idx_sel)*(c_N-1)) + 1;
        c_data(i,:) = cm(id,:);
    else
        c_data(i,:) = [0 0 0];
    end
end


fig = figure('Position', [50 50 1500 600]);
hold on;
grid on;
box on;
bar(1:length(GR_ids), GR_percents(GR_sort_ids), 'EdgeColor', [0 0 0], 'FaceColor', [1 0 0]);
plot([0, length(GR_ids)], GR_thr_perc*[1, 1], 'Color', [0 0 0], 'LineWidth', 2);
xlim([0,length(GR_ids)]);
xticks(1:length(GR_ids));
xticklabels(x_tick_labels);
xlabel('cluster id');
ylabel('percent of data points');
title({strcat('number of clusters:', num2str(max(cl_idx))), ...
    strcat('percent of unclassified points:', num2str(sum(cl_idx==-1)/length(cl_idx)*100), '%'), ...
    strcat('number of selected clusters:', num2str(length(cl_top_idx)), '(with more then', num2str(GR_thr_perc), '% of points)')});
drawnow;
saveas(fig, fullfile(save_fold, 'DBSCAN_clusters.png'));



fig = figure('Position', [50 50 1500 700]);
hold on;
grid on;
box on;
scatter3(score_sel(:,1),  score_sel(:,2), score_sel(:,3), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_first_two_components_DBSCAN_selected.png'));


fig = figure('Position', [50 50 1700 700]);
subplot(2,3,1);
hold on;
grid on;
box on;
scatter(score_sel(:,1),  score_sel(:,2), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
title('PC1 and PC2');

subplot(2,3,2);
hold on;
grid on;
box on;
scatter(score_sel(:,1),  score_sel(:,3), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
title('PC1 and PC3');

subplot(2,3,3);
hold on;
grid on;
box on;
scatter(score_sel(:,1),  score_sel(:,4), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
title('PC1 and PC4');

subplot(2,3,4);
hold on;
grid on;
box on;
scatter(score_sel(:,2),  score_sel(:,3), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
title('PC2 and PC3');

subplot(2,3,5);
hold on;
grid on;
box on;
scatter(score_sel(:,2),  score_sel(:,4), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
title('PC2 and PC4');

subplot(2,3,6);
hold on;
grid on;
box on;
scatter(score_sel(:,3),  score_sel(:,4), 'Marker', '.', ...
    'SizeData', 20, 'CData', c_data);
title('PC3 and PC4');

drawnow;
saveas(fig, fullfile(save_fold, 'PCA_different_components_DBSCAN_selected.png'));


cl_num = max(unique(cl_idx_sel));
for cl_id  = 1:cl_num
    
    c_data = zeros(length(cl_idx_sel), 3);
    for i = 1:length(cl_idx_sel)
        if cl_idx_sel(i) == cl_id
            c_data(i,:) = [1 0 0];
        end
    end
    
    fig = figure('Position', [50 50 1700 700]);
    subplot(2,3,1);
    hold on;
    grid on;
    box on;
    scatter(score_sel(:,1),  score_sel(:,2), 'Marker', '.', ...
        'SizeData', 20, 'CData', c_data);
    title('PC1 and PC2');
    
    subplot(2,3,2);
    hold on;
    grid on;
    box on;
    scatter(score_sel(:,1),  score_sel(:,3), 'Marker', '.', ...
        'SizeData', 20, 'CData', c_data);
    title('PC1 and PC3');
    
    subplot(2,3,3);
    hold on;
    grid on;
    box on;
    scatter(score_sel(:,1),  score_sel(:,4), 'Marker', '.', ...
        'SizeData', 20, 'CData', c_data);
    title('PC1 and PC4');
    
    subplot(2,3,4);
    hold on;
    grid on;
    box on;
    scatter(score_sel(:,2),  score_sel(:,3), 'Marker', '.', ...
        'SizeData', 20, 'CData', c_data);
    title('PC2 and PC3');
    
    subplot(2,3,5);
    hold on;
    grid on;
    box on;
    scatter(score_sel(:,2),  score_sel(:,4), 'Marker', '.', ...
        'SizeData', 20, 'CData', c_data);
    title('PC2 and PC4');
    
    subplot(2,3,6);
    hold on;
    grid on;
    box on;
    scatter(score_sel(:,3),  score_sel(:,4), 'Marker', '.', ...
        'SizeData', 20, 'CData', c_data);
    title('PC3 and PC4');
    
    sgtitle(strcat('selected cluster = ', num2str(cl_id)));
    
    drawnow;
    saveas(fig, fullfile(save_fold, strcat('PCA_different_components_DBSCAN_selected_cluster_', num2str(cl_id), '.png')));
    
    
end


save(fullfile(save_fold, 'DBSCAN.mat'), ...
    'epsilon', 'minpts', 'GR_thr_perc', 'explained_cum', 'f_N', 'p_N', ...
    'pca_sel_n', 'cl_idx', 'GR_counts', 'GR_ids', 'GR_percents', ...
    'GR_sort_ids', 'cl_top_idx', 'cl_idx_sel', 'score', 'score_sel', ...
    'pcaval', 'latent', 'explained', 'features1_comb_norm', 'features1_comb', ...
    'areas', 'min_area_th', 'max_area_th', 'filt_vect', 'features1_comb_filt');


