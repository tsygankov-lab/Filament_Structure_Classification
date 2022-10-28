clear; clc;

dbscan_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\Septin_structures\DBSCAN.mat';
load(dbscan_res_file);

join_cl_res_file = 'E:\Septin_structure_analysis\updated_pictures\data7\Septin_structures\PCA_DBSCAN_join.mat';
load(join_cl_res_file);

bLab{1} = {'006', '109', '110', '112'};
bLab{2} = {'001', '010', '101', '102', '103'};
bLab{3} = {'003', '101', '108', '109', '110'};

bLoc{1} = 'Septin Structures/Control Septin Structures/Ctrl_Sept7_';
bLoc{2} = 'Septin Structures/CEP1 Knockdown Septin Structures/CEP1-KD_Sept7_';
bLoc{3} = 'Septin Structures/CEP1 Overexpression Septin Structures/CEP1-OE_Sept7_';

b_IDX = cell(1,3);

cnt = 0;
cnt2 = 0;
cl_num = max(unique(cl_idx_join));
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
        
        indx = cl_idx_join(cnt2+1:cnt2+ndat2);
        indx_init = cl_idx(cnt2+1:cnt2+ndat2);
        b_IDX{uu}{w} = {};
        for cl_id = 1:cl_num
            im_indx = ismember(Lexcl(:), im_obj_sel_ids(indx == cl_id));
            b_IDX{uu}{w}{cl_id} = im_indx;
        end
        
        cnt = cnt + ndat;
        cnt2 = cnt2 + ndat2;
    end
end

cm = jet(5);
cm = cm(2:end-1,:);
cm = [1, 0.3, 0.3; 0.3, 1, 0.3; 0.3, 0.3, 1];

fr_size = [150, 150];

fig = figure();

uu = 1;
ulab = bLab{uu};
w = 2;
loc = [bLoc{uu} ulab{w}];
mkdir([loc '_PCA_DBSCAN_plot_join_crop']);
fr_loc = [1720 1520];
load([loc '.mat']);
IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
IM_rgb_crop = zeros(fr_size(1), fr_size(2), 3);
for cl_id = 1:cl_num
    im_indx = b_IDX{uu}{w}{cl_id};
    IM_r = zeros(size(Lexcl));
    IM_g = zeros(size(Lexcl));
    IM_b = zeros(size(Lexcl));
    IM_r(im_indx) = cm(cl_id,1);
    %IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1) = 1;
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
    saveas(fig, [loc '_PCA_DBSCAN_plot_join_crop', '\cluster=' num2str(cl_id) '.png']);
    
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
saveas(fig, [loc '_PCA_DBSCAN_plot_join_crop', '\all selected clusters.png']);







uu = 2;
ulab = bLab{uu};
w = 3;
loc = [bLoc{uu} ulab{w}];
mkdir([loc '_PCA_DBSCAN_plot_join_crop']);
fr_loc = [1420 1510];
load([loc '.mat']);
IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
IM_rgb_crop = zeros(fr_size(1), fr_size(2), 3);
for cl_id = 1:cl_num
    im_indx = b_IDX{uu}{w}{cl_id};
    IM_r = zeros(size(Lexcl));
    IM_g = zeros(size(Lexcl));
    IM_b = zeros(size(Lexcl));
    IM_r(im_indx) = cm(cl_id,1);
    %IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1) = 1;
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
    saveas(fig, [loc '_PCA_DBSCAN_plot_join_crop', '\cluster=' num2str(cl_id) '.png']);

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
saveas(fig, [loc '_PCA_DBSCAN_plot_join_crop', '\all selected clusters.png']);






uu = 3;
ulab = bLab{uu};
w = 1;
loc = [bLoc{uu} ulab{w}];
 mkdir([loc '_PCA_DBSCAN_plot_join_crop']);
fr_loc = [1760 920];
load([loc '.mat']);
IM_rgb = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
IM_rgb_crop = zeros(fr_size(1), fr_size(2), 3);
for cl_id = 1:cl_num
    im_indx = b_IDX{uu}{w}{cl_id};
    IM_r = zeros(size(Lexcl));
    IM_g = zeros(size(Lexcl));
    IM_b = zeros(size(Lexcl));
    IM_r(im_indx) = cm(cl_id,1);
    %IM_r(fr_loc(1):fr_loc(1)+fr_size(1)-1, fr_loc(2):fr_loc(2)+fr_size(2)-1) = 1;
    IM_g(im_indx) = cm(cl_id,2);
    IM_b(im_indx) = cm(cl_id,3);
    IM_rgb(:,:,1) = IM_rgb(:,:,1) + IM_r;
    IM_rgb(:,:,2) = IM_rgb(:,:,2) + IM_g;
    IM_rgb(:,:,3) = IM_rgb(:,:,3) + IM_b;
    
% %     IM_rgb_i = zeros(size(Lexcl, 1), size(Lexcl, 2), 3);
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
    saveas(fig, [loc '_PCA_DBSCAN_plot_join_crop', '\cluster=' num2str(cl_id) '.png']);

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
saveas(fig, [loc '_PCA_DBSCAN_plot_join_crop', '\all selected clusters.png']);



