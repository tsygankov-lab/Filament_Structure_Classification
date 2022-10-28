clear; clc;

warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:nearlySingularMatrix');

for phn = 1:3

    if phn == 1
        bThr1 = [16000, 12000, 12000, 12000];
        bThr2 = [2000,  2000,  2000,  2000];
        bLab = { '006', '109', '110', '112'};
    elseif phn == 2
        bThr1 = [24000, 24000, 24000, 24000, 24000];
        bThr2 = [2000,  2000,  2000,  2000,  2000];
        bLab = { '001', '010', '101', '102', '103'};
    elseif phn == 3
        bThr1 = [24000, 20000, 20000, 20000, 20000];
        bThr2 = [2000,  2000,  2000,  2000,  2000];
        bLab = { '003', '101', '108', '109', '110'};
    elseif phn == 4
        bThr1 = [4650,  5500,  5500,  4000];
        bThr2 = [1500,  1500,  1500,  1500];
        bLab = { '001', '002', '005', '016'};
    elseif phn == 5
        bThr1 = [2500,  1,     1,     5000,  2000,  3000];
        bThr2 = [1500,  1500,  1500,  1500,  1500,  1500];
        bLab = { '001', '003', '004', '009', '010', '012'};
    elseif phn == 6
        bThr1 = [2314,  3087,  3856,  4625,  4500,  4500];
        bThr2 = [1500,  1500,  1500,  1500,  1500,  1500];
        bLab = { '007', '012', '013', '021', '022', '024'};
    end
    
    for w = 1:length(bLab)
        if phn == 1
            loc = ['Septin Structures/Control Septin Structures/Ctrl_Sept7_' bLab{w}];
        elseif phn == 2
            loc = ['Septin Structures/CEP1 Knockdown Septin Structures/CEP1-KD_Sept7_' bLab{w}];
        elseif phn == 3
            loc = ['Septin Structures/CEP1 Overexpression Septin Structures/CEP1-OE_Sept7_' bLab{w}];
        elseif phn == 4
            loc = ['CEP1 Structures/Control CEP1 Structures/Ctrl_CEP1-EGFP_' bLab{w}];
        elseif phn == 5
            loc = ['CEP1 Structures/Septin7-Knockdown CEP1 Structures/Sept7-KD_CEP1-EGFP_' bLab{w}];
        elseif phn == 6
            loc = ['CEP1 Structures/100uM FCF CEP1 Structures/100uM FCF_CEP1-EGFP_' bLab{w}];
        end
        disp(loc);
        disp(class(imread([loc '.tif'],1)));
        
        Rdisk = 7;
        Thr1 = bThr1(w);

        Amax = 5000;
        Thr2 = bThr2(w);
        M = 30;

        F0 = double(imread([loc '.tif'],1));
        %F0 = double(im2gray(im2uint16(imread([loc '.tif'],1))));
        F = F0;
        F(F < Thr1) = 0;
        F(F > 0) = 1;
        F = imfill(F,'holes');
        se = strel('disk',Rdisk,0);
        F = imerode(F,se);
        L = bwlabel(F);
        pr = regionprops(L,'Area');
        A = cat(1,pr.Area);
        for i = 1:length(A)
            if A(i)<Amax
                L(L == i) = 0;
            end
        end
        CellMask = imdilate(L>0,se);


        CR = F0;
        BG = imgaussfilt(CR,5);
        J = CR - BG;
        msk = imgaussfilt(J,1);
        msk(msk<Thr2) = 0;
        msk(msk>0) = 1;
        
%         msk_min_area = 20;
%         msk_L = bwlabel(msk);
%         pr = regionprops(msk_L,'Area');
%         A = cat(1,pr.Area);
%         for i = 1:length(A)
%             if A(i)<msk_min_area
%                 msk_L(msk_L == i) = 0;
%             end
%         end
%         msk = double(msk_L>0);

        fig = figure('Position', [50 50 1800 900]);
        clf;
        subplot(1,2,1);
        hold on;
        colormap turbo;
        axis off;
        axis ij;
        axis image;
        imagesc(F0);
        
        subplot(1,2,2);
        hold on;
        colormap turbo;
        axis off;
        axis ij;
        axis image;
        imagesc(0.1*CellMask + 0.5*CellMask.*msk + 0.3*(~CellMask).*msk);
        drawnow;
        
        saveas(fig,[loc '_filt.fig']);
        saveas(fig,[loc '_filt.png']);
        close(fig);
        
        [B,L,nob,~] = bwboundaries(CellMask.*msk,'noholes');

        %%{
        stats = regionprops(L,'Area','ConvexArea','Eccentricity','EquivDiameter','Extent','FilledArea',...
            'MajorAxisLength','MinorAxisLength','Perimeter','Solidity');
        Mstat = 14;
        AR = cat(1,stats.Area);
        CO = cat(1,stats.ConvexArea);
        EC = cat(1,stats.Eccentricity);
        EQ = cat(1,stats.EquivDiameter);
        EX = cat(1,stats.Extent);
        FI = cat(1,stats.FilledArea);
        MA = cat(1,stats.MajorAxisLength);
        MI = cat(1,stats.MinorAxisLength);
        PE = cat(1,stats.Perimeter);
        SO = cat(1,stats.Solidity);
        
        features = zeros(nob,Mstat+2*M);
        %}

        fig1 = figure('Position',get(0,'Screensize'));            
        hold on;
        %Jtmp = J;
        %Jtmp(J<-40) = -40;
        %Jtmp(J> 40) =  40;
        %imagesc(Jtmp);
        imagesc(CellMask.*F0);
        axis equal; axis ij; axis off;


        debg = zeros(size(L));
        err = ones(1,nob);
        mI = zeros(1,nob);
        sI = zeros(1,nob);

        Lexcl = L;

        %tINT = zeros(size(L));  %!!!
        for i = 1:nob

            disp([phn w i/nob]);
            bnd = B{i};
            nbnd = size(bnd,1);

            mI(i) = mean(F0(L==i));
            sI(i) = std(F0(L==i),1);
            %tINT(L==i) = mI(i);  %!!!

            tmp = L;
            tmp(L==i) = 0;
            tmp(tmp>0) = 1;
            dmap = bwdist(tmp);

            dval = zeros(nbnd,1);
            for j = 1:nbnd
                dval(j) = dmap(bnd(j,1),bnd(j,2));
                debg(bnd(j,1),bnd(j,2)) = dval(j);
            end

            [coeff,bnd_rec] = fourier_shape(bnd,dval,M,1);

            plot(bnd(:,2),bnd(:,1),'r');
            if ~isempty(coeff) && sum(isnan(coeff))==0 
                err(i) = sqrt(sum((bnd_rec(:,2) - bnd(:,2)).^2))/nbnd;
                
                features(i,1:Mstat) = [AR(i),CO(i),EC(i),EQ(i),EX(i),FI(i),MA(i),MI(i),PE(i),SO(i),mean(dval),std(dval),mI(i),sI(i)];
                features(i,(Mstat+1):(Mstat+length(coeff))) = coeff';
                if err(i) < 0.1
                    plot(bnd_rec(:,2),bnd_rec(:,1),'y');
                else
                    Lexcl(L==i) = 0;
                end
            else
                Lexcl(L==i) = 0;
            end
        end

        %%{
        features1 = features(err<0.1,:);
        features_norm = features1;

        for i = 1:(Mstat+2*M)

            if std(features1(:,i))>0
                features_norm(:,i) = (features1(:,i) - mean(features1(:,i)))/std(features1(:,i));
            else
                features_norm(:,i) = features1(:,i) - mean(features1(:,i));
            end
        end
        %}

        %fig2 = figure('Position',get(0,'Screensize'));            
        %tdebg = debg;
        %tdebg(debg>30) = 30;
        %imagesc(tdebg);
        %axis equal; axis ij; axis off;

        fig3 = figure('Position',get(0,'Screensize'));
        colormap(gray);
        hold on;
        imagesc(CellMask.*F0);
        axis ij; axis image; axis off;

        %%{
        saveas(fig3,[loc '_3.fig']);
        saveas(fig3,[loc '_3.png']);
        close(fig3);

        %saveas(fig2,[loc '2.fig']);
        %saveas(fig2,[loc '2.png']);
        %close(fig2);

        saveas(fig1,[loc '_1.fig']);
        saveas(fig1,[loc '_1.png']);
        close(fig1);

        save([loc '.mat'],'F0','msk','CellMask','Lexcl','debg','M','Rdisk','Amax','Thr1','Thr2','features','features1','features_norm','err');
        %}
    end
end


