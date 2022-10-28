clear; clc;

dbscan_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures\DBSCAN.mat';
load(dbscan_res_file);

save_fold = 'E:\Septin_structure_analysis\updated_pictures\data7\CEP1_structures';
load(fullfile(save_fold, 'PCA_DBSCAN_refine_GMM_comined.mat'), ...
    'phen_data' , 'cl_idx_comb', 'cell_data', 'cl_num', 'cell_num', ...
    'dens_data', 'dens_phen');

bLab{1} = {'001', '002', '005', '016'};
bLab{2} = {'001', '003', '004', '009', '010', '012'};
bLab{3} = {'007', '012', '013', '021', '022', '024'};

bLoc{1} = 'CEP1 Structures/Control CEP1 Structures/Ctrl_CEP1-EGFP_';
bLoc{2} = 'CEP1 Structures/Septin7-Knockdown CEP1 Structures/Sept7-KD_CEP1-EGFP_';
bLoc{3} = 'CEP1 Structures/100uM FCF CEP1 Structures/100uM FCF_CEP1-EGFP_';

b_IDX = cell(1,3);
b_IDX_rest = cell(1,3);
b_IDX_unclassified = cell(1,3);

cl_top_idx = 1:max(cl_idx_comb);
cl_idx_sel = cl_idx_comb;

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
        b_IDX{uu}{w} = {};
        for cl_id = 1:length(cl_top_idx)
            im_indx = ismember(Lexcl(:), im_obj_sel_ids(indx == cl_id));
            b_IDX{uu}{w}{cl_id} = im_indx;
        end
        
        b_IDX_unclassified{uu}{w} = ismember(Lexcl(:), im_obj_sel_ids(indx_init == -1));
        rest_indx = (~ismember(indx_init, cl_top_idx)) & (indx_init ~= -1);
        b_IDX_rest{uu}{w} = ismember(Lexcl(:), im_obj_sel_ids(rest_indx));
        
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end


fr_size = [500, 500];

cm = jet(cl_num);

fig = figure();

uu = 1;
ulab = bLab{uu};
w = 4;
loc = [bLoc{uu} ulab{w}];
mkdir([loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop']);
fr_loc = [1100 1300];
load([loc '.mat']);
ndat = size(features1,1);
IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
IM_rgb_crop = zeros(fr_size(1), fr_size(2), 3);
for cl_id_i = 1:length(cl_top_idx)
    cl_id = cl_top_idx(cl_id_i);
    im_indx = b_IDX{uu}{w}{cl_id_i};
    IM_r = zeros(size(Lexcl));
    IM_g = zeros(size(Lexcl));
    IM_b = zeros(size(Lexcl));
    IM_r(im_indx) = cm(cl_id,1);
%     IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1) = 1;
    IM_g(im_indx) = cm(cl_id,2);
    IM_b(im_indx) = cm(cl_id,3);
    IM_rgb(:,:,1) = IM_rgb(:,:,1) + IM_r;
    IM_rgb(:,:,2) = IM_rgb(:,:,2) + IM_g;
    IM_rgb(:,:,3) = IM_rgb(:,:,3) + IM_b;
    
    IM_rgb_i = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
    IM_rgb_i(:,:,1) = IM_r;
    IM_rgb_i(:,:,2) = IM_g;
    IM_rgb_i(:,:,3) = IM_b;
    
    IM_r_crop = IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_g_crop = IM_g(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_b_crop = IM_b(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_rgb_crop(:,:,1) = IM_rgb_crop(:,:,1) + IM_r_crop;
    IM_rgb_crop(:,:,2) = IM_rgb_crop(:,:,2) + IM_g_crop;
    IM_rgb_crop(:,:,3) = IM_rgb_crop(:,:,3) + IM_b_crop;
    
    IM_rgb_crop_i = zeros(fr_size(1), fr_size(2), 3);
    IM_rgb_crop_i(:,:,1) = IM_r_crop;
    IM_rgb_crop_i(:,:,2) = IM_g_crop;
    IM_rgb_crop_i(:,:,3) = IM_b_crop;
    
    clf;
    hold on;
    axis ij;
    axis off;
    image(IM_rgb_crop_i);
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    set(gca,'position',[0 0 1 1]);
    set(gca, 'XLim', [0 fr_size(1)], 'YLim', [0 fr_size(2)]);
    set(fig,'PaperUnits','points');
    set(fig,'PaperSize',[fr_size(1) fr_size(1)]);
    set(fig,'PaperPosition',[0 0 fr_size(1) fr_size(2)]);
    set(gca,'Color','k');
    fig.Color = 'black';
    fig.InvertHardcopy = 'off';
    drawnow;
    saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop', '\cluster=' num2str(cl_id) '.png']);

end
% set(fig, 'Position',get(0,'Screensize'));
% clf;
% hold on;
% axis ij;
% axis image;
% axis off;
% image(IM_rgb);
% title('all selected clusters');
% drawnow;

clf;
hold on;
axis ij;
axis off;
image(IM_rgb_crop);
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
set(gca,'position',[0 0 1 1]);
set(gca, 'XLim', [0 fr_size(1)], 'YLim', [0 fr_size(2)]);
set(fig,'PaperUnits','points');
set(fig,'PaperSize',[fr_size(1) fr_size(1)]);
set(fig,'PaperPosition',[0 0 fr_size(1) fr_size(2)]);
set(gca,'Color','k');
fig.Color = 'black';
fig.InvertHardcopy = 'off';
drawnow;
saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop', '\all selected clusters.png']);





uu = 2;
ulab = bLab{uu};
w = 3;
loc = [bLoc{uu} ulab{w}];
mkdir([loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop']);
fr_loc = [1100 850];
load([loc '.mat']);
ndat = size(features1,1);
IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
IM_rgb_crop = zeros(fr_size(1), fr_size(2), 3);
for cl_id_i = 1:length(cl_top_idx)
    cl_id = cl_top_idx(cl_id_i);
    im_indx = b_IDX{uu}{w}{cl_id_i};
    IM_r = zeros(size(Lexcl));
    IM_g = zeros(size(Lexcl));
    IM_b = zeros(size(Lexcl));
    IM_r(im_indx) = cm(cl_id,1);
%     IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1) = 1;
    IM_g(im_indx) = cm(cl_id,2);
    IM_b(im_indx) = cm(cl_id,3);
    IM_rgb(:,:,1) = IM_rgb(:,:,1) + IM_r;
    IM_rgb(:,:,2) = IM_rgb(:,:,2) + IM_g;
    IM_rgb(:,:,3) = IM_rgb(:,:,3) + IM_b;
    
    IM_rgb_i = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
    IM_rgb_i(:,:,1) = IM_r;
    IM_rgb_i(:,:,2) = IM_g;
    IM_rgb_i(:,:,3) = IM_b;
    
    IM_r_crop = IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_g_crop = IM_g(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_b_crop = IM_b(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_rgb_crop(:,:,1) = IM_rgb_crop(:,:,1) + IM_r_crop;
    IM_rgb_crop(:,:,2) = IM_rgb_crop(:,:,2) + IM_g_crop;
    IM_rgb_crop(:,:,3) = IM_rgb_crop(:,:,3) + IM_b_crop;
    
    IM_rgb_crop_i = zeros(fr_size(1), fr_size(2), 3);
    IM_rgb_crop_i(:,:,1) = IM_r_crop;
    IM_rgb_crop_i(:,:,2) = IM_g_crop;
    IM_rgb_crop_i(:,:,3) = IM_b_crop;
    
    clf;
    hold on;
    axis ij;
    axis off;
    image(IM_rgb_crop_i);
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    set(gca,'position',[0 0 1 1]);
    set(gca, 'XLim', [0 fr_size(1)], 'YLim', [0 fr_size(2)]);
    set(fig,'PaperUnits','points');
    set(fig,'PaperSize',[fr_size(1) fr_size(1)]);
    set(fig,'PaperPosition',[0 0 fr_size(1) fr_size(2)]);
    set(gca,'Color','k');
    fig.Color = 'black';
    fig.InvertHardcopy = 'off';
    drawnow;
    saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop', '\cluster=' num2str(cl_id) '.png']);

end
% set(fig, 'Position',get(0,'Screensize'));
% clf;
% hold on;
% axis ij;
% axis image;
% axis off;
% image(IM_rgb);
% title('all selected clusters');
% drawnow;

clf;
hold on;
axis ij;
axis off;
image(IM_rgb_crop);
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
set(gca,'position',[0 0 1 1]);
set(gca, 'XLim', [0 fr_size(1)], 'YLim', [0 fr_size(2)]);
set(fig,'PaperUnits','points');
set(fig,'PaperSize',[fr_size(1) fr_size(1)]);
set(fig,'PaperPosition',[0 0 fr_size(1) fr_size(2)]);
set(gca,'Color','k');
fig.Color = 'black';
fig.InvertHardcopy = 'off';
drawnow;
saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop', '\all selected clusters.png']);




uu = 3;
ulab = bLab{uu};
w = 2;
loc = [bLoc{uu} ulab{w}];
mkdir([loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop']);
fr_loc = [550 1250];
load([loc '.mat']);
ndat = size(features1,1);
IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
IM_rgb_crop = zeros(fr_size(1), fr_size(2), 3);
for cl_id_i = 1:length(cl_top_idx)
    cl_id = cl_top_idx(cl_id_i);
    im_indx = b_IDX{uu}{w}{cl_id_i};
    IM_r = zeros(size(Lexcl));
    IM_g = zeros(size(Lexcl));
    IM_b = zeros(size(Lexcl));
    IM_r(im_indx) = cm(cl_id,1);
%     IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1) = 1;
    IM_g(im_indx) = cm(cl_id,2);
    IM_b(im_indx) = cm(cl_id,3);
    IM_rgb(:,:,1) = IM_rgb(:,:,1) + IM_r;
    IM_rgb(:,:,2) = IM_rgb(:,:,2) + IM_g;
    IM_rgb(:,:,3) = IM_rgb(:,:,3) + IM_b;
    
    IM_rgb_i = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
    IM_rgb_i(:,:,1) = IM_r;
    IM_rgb_i(:,:,2) = IM_g;
    IM_rgb_i(:,:,3) = IM_b;
    
    IM_r_crop = IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_g_crop = IM_g(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_b_crop = IM_b(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1);
    IM_rgb_crop(:,:,1) = IM_rgb_crop(:,:,1) + IM_r_crop;
    IM_rgb_crop(:,:,2) = IM_rgb_crop(:,:,2) + IM_g_crop;
    IM_rgb_crop(:,:,3) = IM_rgb_crop(:,:,3) + IM_b_crop;
    
    IM_rgb_crop_i = zeros(fr_size(1), fr_size(2), 3);
    IM_rgb_crop_i(:,:,1) = IM_r_crop;
    IM_rgb_crop_i(:,:,2) = IM_g_crop;
    IM_rgb_crop_i(:,:,3) = IM_b_crop;
    
    clf;
    hold on;
    axis ij;
    axis off;
    image(IM_rgb_crop_i);
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    set(gca,'position',[0 0 1 1]);
    set(gca, 'XLim', [0 fr_size(1)], 'YLim', [0 fr_size(2)]);
    set(fig,'PaperUnits','points');
    set(fig,'PaperSize',[fr_size(1) fr_size(1)]);
    set(fig,'PaperPosition',[0 0 fr_size(1) fr_size(2)]);
    set(gca,'Color','k');
    fig.Color = 'black';
    fig.InvertHardcopy = 'off';
    drawnow;
    saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop', '\cluster=' num2str(cl_id) '.png']);

end
% set(fig, 'Position',get(0,'Screensize'));
% clf;
% hold on;
% axis ij;
% axis image;
% axis off;
% image(IM_rgb);
% title('all selected clusters');
% drawnow;

clf;
hold on;
axis ij;
axis off;
image(IM_rgb_crop);
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
set(gca,'position',[0 0 1 1]);
set(gca, 'XLim', [0 fr_size(1)], 'YLim', [0 fr_size(2)]);
set(fig,'PaperUnits','points');
set(fig,'PaperSize',[fr_size(1) fr_size(1)]);
set(fig,'PaperPosition',[0 0 fr_size(1) fr_size(2)]);
set(gca,'Color','k');
fig.Color = 'black';
fig.InvertHardcopy = 'off';
drawnow;
saveas(fig, [loc '_PCA_DBSCAN_refine_GMM_comb_plot_crop', '\all selected clusters.png']);

