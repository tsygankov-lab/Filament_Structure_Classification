clear; clc;

save_fold = 'E:\Septin_structure_analysis\updated_pictures\data7\Septin_structures';

dbscan_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\Septin_structures\DBSCAN.mat';
load(dbscan_res_file);

bLab{1} = {'006', '109', '110', '112'};
bLab{2} = {'001', '010', '101', '102', '103'};
bLab{3} = {'003', '101', '108', '109', '110'};

bLoc{1} = 'Septin Structures/Control Septin Structures/Ctrl_Sept7_';
bLoc{2} = 'Septin Structures/CEP1 Knockdown Septin Structures/CEP1-KD_Sept7_';
bLoc{3} = 'Septin Structures/CEP1 Overexpression Septin Structures/CEP1-OE_Sept7_';

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


cl_1 = [4, 1, 5];
cl_2 = [2, 3, 6];
cl_3 = [7, 8, 9];

cl_idx_comb(ismember(cl_idx_sel, cl_1)) = 1;
cl_idx_comb(ismember(cl_idx_sel, cl_2)) = 2;
cl_idx_comb(ismember(cl_idx_sel, cl_3)) = 3;

cl_idx_join = cl_idx_comb;

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

uAreas = {};
uAreas_all = {};
cnt = 0;
cnt2 = 0;
for uu = 1:3
    ulab = bLab{uu};
    cell_ids = unique(cell_data(phen_data == uu));
    uAreas{uu} = {};
    uAreas_all{uu} = {};
    for w = 1:length(cell_ids)
        cell_i = cell_ids(w);
        loc = [bLoc{uu} ulab{w}];
        load([loc '.mat']);
        ndat = size(features1,1);
        im_obj_ids = unique(sort(Lexcl(Lexcl>0)));
        sel_indx = filt_vect(cnt+1:cnt+ndat);
        im_obj_sel_ids = im_obj_ids(sel_indx);
        ndat2 = sum(sel_indx);
        struct_all_idx = (phen_data == uu) & (cell_data == cell_i);
        cl_cell_idx_comb = cl_idx_comb(struct_all_idx);
        uAreas{uu}{w} = {};
        uAreas_all{uu}{w} = 0;
        for cl_id_i = 1:cl_num
            disp([uu, w, cl_id_i]);
            im_indx = ismember(Lexcl(:), im_obj_sel_ids(cl_cell_idx_comb == cl_id_i));
            
            im = zeros(size(Lexcl));
            im(im_indx) = 1;
            rp = regionprops(im>0, 'Area');
            areas_i = cat(1,rp.Area);
            uAreas{uu}{w}{cl_id_i} = areas_i;
            uAreas_all{uu}{w} = uAreas_all{uu}{w} + sum(areas_i);
            
        end
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end




leg_plts = [];

cm = turbo(5);
cm = cm(2:end-1,:);

fig = figure('Position', [50 50 900 900]);
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
        for w = 1:length(unique((cell_data(phen_data == uu))))
            dens_vals{cl_id} = [dens_vals{cl_id}, sum(uAreas{uu}{w}{cl_id})/uAreas_all{uu}{w}];
        end
        
        dens_mean(cl_id) = mean(dens_vals{cl_id});
        ts = tinv([0.025  0.975], length(dens_vals{cl_id})-1);
        SEM = std(dens_vals{cl_id})/sqrt(length(dens_vals{cl_id}));
        dens_CI_1(cl_id) = dens_mean(cl_id) + ts(1)*SEM;
        dens_CI_2(cl_id) = dens_mean(cl_id) + ts(2)*SEM;
        
        col = cm(uu,:);
        
        eb = errorbar(cl_id+uu/10, dens_mean(cl_id), ...
        dens_CI_1(cl_id)-dens_mean(cl_id), dens_CI_2(cl_id)-dens_mean(cl_id), ...
        'Marker', '.', 'MarkerSize', 20, 'LineStyle', 'none', 'LineWidth', 2, 'Color', col, 'CapSize', 10);
        leg_plts(uu) = eb;
        drawnow;
        
    end
end
xticks([1:cl_num]+1/5);
xticklabels([1:cl_num]);
leg_names = {'Septin Structures Control', ...
    'Septin Structures CEP1 Knockdown', ...
    'Septin Structures CEP1 Overexpression'};
legend(leg_plts, leg_names, 'Location', 'SouthWest');
xticklabels({'Class 1', 'Class 2', 'Class 3'});
ylabel('Structure density');
set(gca, 'FontSize', 20);

saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_join_statistics.png'));



cm = gray(3);

fig = figure('Position', [50 50 1500 900]);
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
        for w = 1:length(unique((cell_data(phen_data == uu))))
            dens_vals{cl_id} = [dens_vals{cl_id}, sum(uAreas{uu}{w}{cl_id})/uAreas_all{uu}{w}];
        end
        
        dens_mean(cl_id) = mean(dens_vals{cl_id});
        ts = tinv([0.025  0.975], length(dens_vals{cl_id})-1);
        SEM = std(dens_vals{cl_id})/sqrt(length(dens_vals{cl_id}));
        dens_CI_1(cl_id) = dens_mean(cl_id) + ts(1)*SEM;
        dens_CI_2(cl_id) = dens_mean(cl_id) + ts(2)*SEM;
        
        col = cm(uu,:);
        
        br = bar(cl_id+uu/10, dens_mean(cl_id), 'BarWidth', 0.1, 'EdgeColor', [0 0 0], 'FaceColor', col, 'LineWidth', 2);
        
        eb = errorbar(cl_id+uu/10, dens_mean(cl_id), ...
        dens_CI_1(cl_id)-dens_mean(cl_id), dens_CI_2(cl_id)-dens_mean(cl_id), ...
        'Marker', '.', 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 2, 'Color', [0 0 0], 'CapSize', 10);
        leg_plts(uu) = br;
        drawnow;
        
    end
end

xticks([1:cl_num]+1/5);
xticklabels([1:cl_num]);
leg_names = {'Septin Structures Control', ...
    'Septin Structures CEP1 Knockdown', ...
    'Septin Structures CEP1 Overexpression'};
%legend(leg_plts, leg_names, 'Location', 'NorthEast');
xticklabels({'Class 1', 'Class 2', 'Class 3'});
ylabel('Structure density');
set(gca, 'FontSize', 35);
set(gca, 'LineWidth', 2);
ylim([0, 0.65]);

saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_join_statistics_barplot.png'));


cell_num = max(cell_data);

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
            dens_data(cell_counter, cl_id) = sum(uAreas{uu}{w}{cl_id})/uAreas_all{uu}{w};
        end
        dens_phen(cell_counter) = uu;
        cell_counter = cell_counter + 1;
    end
end




f_vals = (dens_data - mean(dens_data,1))./std(dens_data,1);

[pcaval, score, latent, ~, explained] = pca(f_vals, 'Algorithm','svd');

explained_cum = cumsum(explained);

f_N = size(f_vals, 2);
p_N = size(f_vals, 1);


c_data = zeros(length(dens_phen), 3);
for i = 1:length(dens_phen)
    c_data(i,:) = cm(dens_phen(i),:);
end

phen_ids_1 = find(dens_phen == 1);
phen_ids_2 = find(dens_phen == 2);
phen_ids_3 = find(dens_phen == 3);

leg_names = {'Septin Structures Control', ...
    'Septin Structures CEP1 Knockdown', ...
    'Septin Structures CEP1 Overexpression'};

sz = 100;

fig = figure('Position', [50 50 900 900]);
hold on;
grid on;
box on;
scatter(score(:,1),  score(:,2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data, 'LineWidth', 2);
sc1 = scatter(score(phen_ids_1(1),1),  score(phen_ids_1(1),2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data(phen_ids_1(1),:), 'LineWidth', 2);
scatter(mean(score(phen_ids_1,1)),  mean(score(phen_ids_1,2)), 'Marker', 'o', ...
    'SizeData', 3*sz, 'MarkerEdgeColor', c_data(phen_ids_1(1),:), ...
    'MarkerFaceColor', c_data(phen_ids_1(1),:), 'LineWidth', 2);
sc2 = scatter(score(phen_ids_2(1),1),  score(phen_ids_2(1),2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data(phen_ids_2(1),:), 'LineWidth', 2);
scatter(mean(score(phen_ids_2,1)),  mean(score(phen_ids_2,2)), 'Marker', 'o', ...
    'SizeData', 3*sz, 'MarkerEdgeColor', c_data(phen_ids_2(1),:), ...
    'MarkerFaceColor', c_data(phen_ids_2(1),:), 'LineWidth', 2);
sc3 = scatter(score(phen_ids_3(1),1),  score(phen_ids_3(1),2), 'Marker', 'o', ...
    'SizeData', sz, 'CData', c_data(phen_ids_3(1),:), 'LineWidth', 2);
scatter(mean(score(phen_ids_3,1)),  mean(score(phen_ids_3,2)), 'Marker', 'o', ...
    'SizeData', 3*sz, 'MarkerEdgeColor', c_data(phen_ids_3(1),:), ...
    'MarkerFaceColor', c_data(phen_ids_3(1),:), 'LineWidth', 2);
legend([sc1, sc2, sc3], leg_names, 'Location', 'NorthEast');
xlabel('PC1');
ylabel('PC2');
set(gca, 'FontSize', 20);
drawnow;

saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_join_pca_first_two_components_area_density_centroids.png'));


vbls = {'Class 1', 'Class 2', 'Class 3'};

fig  = figure('Position', [50 50 900 900]);
hold on;
grid on;
box on;
axis equal;
h = biplot(pcaval(:,1:2), 'Scores', score(:,1:2), 'VarLabels', vbls, 'MarkerSize', 30);
hID = get(h, 'tag'); 
hPt = h(strcmp(hID,'obsmarker')); 
for i = 1:size(score,1)
    set(hPt(i), 'Color', c_data(i,:));
end
hMk = h(strcmp(hID,'varmarker')); 
for i = 1:length(hMk)
    set(hMk(i), 'Color', [0 0 0]);
end
hLn = h(strcmp(hID,'varline')); 
for i = 1:length(hLn)
    set(hLn(i), 'Color', [0 0 0]);
end
set(gca,'FontSize',18);
hLb = h(strcmp(hID,'varlabel')); 
for i = 1:length(hLn)
    set(hLb(i), 'FontSize',18);
end
legend([hPt(find(dens_phen==1,1)), hPt(find(dens_phen==2,1)), hPt(find(dens_phen==3,1))], ...
    leg_names, 'Location', 'northoutside');
set(gca,'linewidth',3);
set(gca, 'FontSize', 20);

saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_join_pca_first_two_components_area_density_biplot.png'));


save(fullfile(save_fold, 'PCA_DBSCAN_join.mat'), ...
    'phen_data', 'cl_1', 'cl_2', 'cl_3', 'cl_idx_join', 'cell_data', ...
    'cl_num', 'cell_num', 'dens_data', 'dens_phen', 'uAreas', 'uAreas_all', 'cm');




cm = zeros(length(dens_phen), 3);
for i = 1:length(dens_phen)
    if dens_phen(i) == 1
        cm(i,:) = [0 0 0];
    end
    if dens_phen(i) == 2
        cm(i,:) = [0.5 0.5 0.5];
    end
    if dens_phen(i) == 3
        cm(i,:) = [1 1 1];
    end
end

fig = figure('Position', [50 50 800 800]);
hold on;
grid on;
box on;
sc1 = scatter3(dens_data(dens_phen==1,1), dens_data(dens_phen==1,2), dens_data(dens_phen==1,3), ...
    'CData', cm(dens_phen==1,:), 'Marker', 'o', 'SizeData', 150, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
sc2 = scatter3(dens_data(dens_phen==2,1), dens_data(dens_phen==2,2), dens_data(dens_phen==2,3), ...
    'CData', cm(dens_phen==2,:), 'Marker', 'o', 'SizeData', 150, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.5 0.5 0.5]);
sc3 = scatter3(dens_data(dens_phen==3,1), dens_data(dens_phen==3,2), dens_data(dens_phen==3,3), ...
    'CData', cm(dens_phen==3,:), 'Marker', 'o', 'SizeData', 150, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 1]);
xlabel('Class 1');
ylabel('Class 2');
zlabel('Class 3');
view(45,30);
%legend([sc1, sc2, sc3], leg_names, 'Location', 'NorthEast');
set(gca, 'FontSize', 20);

saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_join_area_density_3D_scatter.png'));