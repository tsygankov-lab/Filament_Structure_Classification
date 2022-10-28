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

b_CL_struct_areas = cell(1,3);
b_CL_cell_areas = cell(1,3);

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
        im_obj_ids = unique(sort(Lexcl(Lexcl>0)));
        sel_indx = filt_vect(cnt+1:cnt+ndat);
        im_obj_sel_ids = im_obj_ids(sel_indx);
        ndat2 = sum(sel_indx);
        
        indx = cl_idx_sel(cnt2+1:cnt2+ndat2);
        indx_init = cl_idx(cnt2+1:cnt2+ndat2);
        
        b_CL_cell_areas{uu}{w} = sum(CellMask(:));
        b_CL_struct_areas{uu}{w} = {};
        for cl_id = 1:length(cl_top_idx)
            im_indx = ismember(Lexcl(:), im_obj_sel_ids(indx == cl_id));
            im = zeros(size(msk));
            im(im_indx) = 1;
            
            im_con = bwconncomp(im);
            b_CL_struct_areas{uu}{w}{cl_id} = [];
            for i = 1:im_con.NumObjects
                b_CL_struct_areas{uu}{w}{cl_id} = [b_CL_struct_areas{uu}{w}{cl_id}, length(im_con.PixelIdxList{i})];
            end
        end
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end

struct_all_areas = cell(1, cl_num);
for uu = 1:3
    ulab = bLab{uu}; 
    for w = 1:length(ulab)
        for cl_id = 1:length(cl_top_idx)
            struct_all_areas{cl_id} = [struct_all_areas{cl_id},  b_CL_struct_areas{uu}{w}{cl_id}];
        end
        cnt = cnt + ndat;
    end
end
struct_mean_areas = zeros(1, cl_num);
struct_mean_areas_CI_1 = zeros(1, cl_num);
struct_mean_areas_CI_2 = zeros(1, cl_num);
for cl_id = 1:cl_num
    struct_mean_areas(cl_id) = mean(struct_all_areas{cl_id});
    ts = tinv([0.025  0.975], length(struct_all_areas{cl_id})-1);
    SEM = std(struct_all_areas{cl_id})/sqrt(length(struct_all_areas{cl_id}));
    struct_mean_areas_CI_1(cl_id) = struct_mean_areas(cl_id) + ts(1)*SEM;
    struct_mean_areas_CI_2(cl_id) = struct_mean_areas(cl_id) + ts(2)*SEM;
end
 
[~, area_sort_idx] = sort(struct_mean_areas);


x_tick_labels = cell(1,1);
for i = 1:cl_num
    x_tick_labels{i} = num2str(cl_top_idx(i));
end

fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
bar(1:cl_num, struct_mean_areas(area_sort_idx), 'EdgeColor', [0 0 0], 'FaceColor', [1 0 0]);
xticks(1:cl_num);
xticklabels(x_tick_labels((area_sort_idx)));
xlabel('cluster id');
ylabel('average structure area');
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_clusters_by_average_structure_area.png'));


fig = figure('Position', [50 50 500 500]);
hold on;
grid on;
box on;
eb = errorbar(1:cl_num, struct_mean_areas(area_sort_idx), ...
    -struct_mean_areas_CI_1(area_sort_idx)+struct_mean_areas(area_sort_idx), ...
    -struct_mean_areas_CI_1(area_sort_idx)+struct_mean_areas(area_sort_idx), ...
    'Marker', '.', 'LineStyle', 'none', 'LineWidth', 3, 'Color', [1 0 0]);
xticks(1:cl_num);
xticklabels(x_tick_labels((area_sort_idx)));
xlabel('cluster id');
ylabel('average structure area');
drawnow;
saveas(fig, fullfile(save_fold, 'PCA_DBSCAN_clusters_by_average_structure_area_error_bar.png'));


for uu = 1:3
    ulab = bLab{uu}; 
    dens_vals = cell(1, cl_num);
    dens_mean = zeros(1, cl_num);
    dens_CI_1 = zeros(1, cl_num);
    dens_CI_2 = zeros(1, cl_num);
    
    for cl_id = 1:cl_num
        for w = 1:length(ulab)
            dens_vals{cl_id} = [dens_vals{cl_id}, length(b_CL_struct_areas{uu}{w}{cl_id})/b_CL_cell_areas{uu}{w}];
        end
        
        dens_mean(cl_id) = mean(dens_vals{cl_id});
        ts = tinv([0.025  0.975], length(dens_vals{cl_id})-1);
        SEM = std(dens_vals{cl_id})/sqrt(length(dens_vals{cl_id}));
        dens_CI_1(cl_id) = dens_mean(cl_id) + ts(1)*SEM;
        dens_CI_2(cl_id) = dens_mean(cl_id) + ts(2)*SEM;
    end
    
    fig = figure('Position', [50 50 1000 500]);
    hold on;
    grid on;
    box on;
    eb = errorbar(1:cl_num, dens_mean(area_sort_idx), ...
        dens_CI_1-dens_mean, dens_CI_2-dens_mean, ...
        'Marker', '.', 'MarkerSize', 20, 'LineStyle', 'none', 'LineWidth', 1, 'Color', [1 0 0]);
    xticks(1:cl_num);
    xticklabels(x_tick_labels((area_sort_idx)));
    xlabel('cluster id');
    ylabel('density of structures');
    %ylim([0 2.5*10^-3]);
    title(bLoc{uu});
    drawnow;
    saveas(fig, fullfile(save_fold, strcat('PCA_DBSCAN_uu_', num2str(uu), '.png')));
end
