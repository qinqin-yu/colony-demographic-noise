function [xplot, yplot] = spt_one_timestep_2(image_info, t, maxd, barrg)
%function [lub] = spt()
%%% Parameters
%%%% initial feature detections
%featsize=2;%4         % radius of the particule (px), dected by a mask of size around 2 or 3 px (init=3)
%masscut=25000; 	      % cutting intensity for particule detection (init=100000)
%%%% minimal distance beetwen beads (first image)
tolproxy=1;
%%%% other rejection processes
barcc=1.0;          % excentricity (0=circular, 1=linear), usually around 0.1 or 0.2 (init=0.1)
barrg=barrg;%5;%5;            % gyration radius (linked to the radius of particule in px) (init=5)
IdivRg=0;      % minimal integrated intensity required to detect a particule (init=0)
%%%% cutoff on max displacement
%maxd=15;%6;%4;%12;        % max displacement of a bead between 2 frames inside the mask (px)

barcc_s=0.2;

Nb_frames=1;

%Import files

%path=pwd;
%Ia=imread(fullfile(path,'0010.tif'));
%Ib=imread(fullfile(path,'0011.tif'));
%Ic=imread(fullfile(path,'0012.tif'));

% path_folder = ['/Users/qinqinyu/Documents/hallatschek_lab/data_changed/' ...
%     '20190129_keio_test/'];
% 
% %path_folder = ['../../../data_changed/' ...
% %    '20181118_ecoli_beads_high_mag/'];
% 
% oldpath = path;
% newpath = strcat(path_folder, '3595_focused/rfp_1-4');
% path(oldpath, newpath);
% 
% %Import cropped images (where best focus will be calculated) into a cell
% imagefiles = dir(strcat(path_folder, '3595_focused/rfp_1-4/*.tif'));      
% nfiles = length(imagefiles);    % Number of files found
% 
% image_info = cell(nfiles, 3);   
% %image matix | t | z | c
% 
% %Add information about the images from the file names
% %figure
% for ii=1:nfiles
% %     if mod(ii, 10) == 0
% %         ii
% %     end
%     ii
%    currentfilename = imagefiles(ii).name;
%    currentimage = imread(currentfilename);
%    image_info{ii, 2} = str2num(currentfilename(1:end-4));   %t
%    %image_info{ii, 3} = str2num(currentfilename(8:10));  %z
%    %image_info{ii, 4} = str2num(currentfilename(13:15)); %c
%    %image_info{ii, 5} = std(im2double(currentimage(:)));
%    
%    %Increase contrast
%    background = imopen(currentimage,strel('disk',12));
%    currentimage2=imsubtract(currentimage,background);
%    contrastedimage=double(contrast(currentimage2,0.1,0.1));
% 
%    image_info{ii, 1} = contrastedimage;    %image matrix
%    
%    M_PA = feature2D(contrastedimage,1,featsize,masscut);
%    image_info{ii, 3} = M_PA;
%    %subplot(1,nfiles,ii), hold on, imshow(contrastedimage,[]), plot(M_PA(:,1),M_PA(:,2),'+r','MarkerSize',8)
% end
% 
% path(oldpath)

tmat = cell2mat(image_info(:, 2));

%[~, a_ind] = min(t);
a_ind = find(tmat == t);

b_ind = find(tmat == t + 1);

aM_PA = image_info{a_ind, 3};
%b_ind
%image_info
bM_PA = image_info{b_ind, 3};

% % Increase contrast
%     
% background = imopen(Ia,strel('disk',12));
% Ia2=imsubtract(Ia,background);
% background = imopen(Ib,strel('disk',12));
% Ib2=imsubtract(Ib,background);
% background = imopen(Ic,strel('disk',12));
% Ic2=imsubtract(Ic,background);
% 
% a=double(contrast(Ia2,0.1,0.1));
% b=double(contrast(Ib2,0.1,0.1));
% c=double(contrast(Ic2,0.1,0.1));

% % Detect the features
% aM_PA = feature2D(a,1,featsize,masscut);
% bM_PA = feature2D(b,1,featsize,masscut);%*3/4); %if you want to be less restrictive on the 2nd image (in case you impose a hard cutoff on the distance)
% cM_PA = feature2D(c,1,featsize,masscut);
% 
% figure, 
% subplot(1,3,1), hold on, imshow(a,[]), plot(aM_PA(:,1),aM_PA(:,2),'+r','MarkerSize',8)
% subplot(1,3,2), hold on, imshow(b,[]), plot(bM_PA(:,1),bM_PA(:,2),'+r','MarkerSize',8) %fig beads detect 2%
% subplot(1,3,3), hold on, imshow(c,[]), plot(cM_PA(:,1),cM_PA(:,2),'+r','MarkerSize',8) %fig beads detect 2%

%%
% Feature rejection process on first image
if isempty(aM_PA)==0 %MC 11juin2014
    
    %%%% Exclude beads that are closer to one another than tolproxy
    tolerance=tolproxy;
    distMatrix=(tolerance+1)*ones(size(aM_PA,1),size(aM_PA,1));
    for iMat=1:size(aM_PA,1)
        for jMat=1:iMat
            if iMat==jMat
                distMatrix(iMat,jMat)=tolerance+1;
            else
                distMatrix(iMat,jMat)=sqrt((aM_PA(iMat,1)-aM_PA(jMat,1))^2+(aM_PA(iMat,2)-aM_PA(jMat,2))^2);
            end
        end
    end
    
    idx_line_col=[];
    for iline=1:size(aM_PA,1)
        idxcol=find(distMatrix(iline,:)<tolerance);
        if isempty(idxcol)==0
            idx_line_col=vertcat(idx_line_col,[iline*ones(size(idxcol,2),1), idxcol.']);
        end
    end
    
    duplicate=idx_line_col(:);
    %                 duplicate=[];
    %                 for i_idxlinecol=1:size(idx_line_col,1)
    %                     dipl=zeros(1,2);
    %                     dipl(:)=dX(idx_line_col(i_idxlinecol,:)).^2+dY(idx_line_col(i_idxlinecol,:)).^2;
    %                     [maxDipl, idx_maxDipl]=max(dipl);
    %                     duplicate=vertcat(duplicate,idx_line_col(i_idxlinecol,idx_maxDipl));
    %                 end
    duplicate=unique(duplicate);
    aM_PA(duplicate,:)=[];
    clear duplicate distMatrix tolerance
    %subplot(1,nfiles,1), hold on, plot(aM_PA(:,1),aM_PA(:,2),'xb','MarkerSize',8) %fig beads detect3%
    % Modif sorting MC Sept 2014
    
    %%%% Other rejection processes
    X_PA= aM_PA(:,5)>barcc;
    aM_PA(X_PA,1:5)=0;
    X_PA= aM_PA(:,4)>barrg;
    aM_PA(X_PA,1:5)=0;
    X_PA= aM_PA(:,3)./aM_PA(:,4)<IdivRg;
    aM_PA(X_PA,1:5)=0;
    aM_PA=aM_PA(aM_PA(:,1)~=0,:);
    %subplot(1,nfiles,1), hold on, plot(aM_PA(:,1),aM_PA(:,2),'og','MarkerSize',8) %fig beads detect4%
    image_info{a_ind, 3} = aM_PA;
    
    %%%% Processing second image, can be less strict on cutoffs
    if isempty(bM_PA)==0 %MC 11juin2014

        %Rejection processes
        %X_PA= bM_PA(:,5)>barcc_s;  % uncommented in original 
        %bM_PA(X_PA,1:5)=0; % uncommented in original 
        %                     X_PA= sM_PA(:,4)>barrg;
        %                     sM_PA(X_PA,1:5)=0;
        %                     X_PA= sM_PA(:,3)./sM_PA(:,4)<IdivRg;
        %                     sM_PA(X_PA,1:5)=0;
        %                     sM_PA=sM_PA(sM_PA(:,1)~=0,:);
        
        %%%% Links the 2 images
        %Mtrack_PA = [];
        %for jj = 1:nfiles
            %jj
            %M_PA = image_info{jj, 3};
            
            %Mtrack_PA = [Mtrack_PA; M_PA(:,1:2),jj*ones(size(M_PA,1),1),jj*ones(size(M_PA,1),1)];
        Mtrack_PA=[aM_PA(:,1:2),ones(size(aM_PA,1),1),ones(size(aM_PA,1),1); bM_PA(:,1:2),2*ones(size(bM_PA,1),1),2*ones(size(bM_PA,1),1)];...
        %    cM_PA(:,1:2),3*ones(size(cM_PA,1),1),3*ones(size(cM_PA,1),1)];
        %                     var_Mtrack=genvarname(['Mtrack_' num2str(ff,'%03d') '_' num2str(i) num2str(j)]);
        %                     eval(['Mtrack_' num2str(ff,'%02d') '_' num2str(i) num2str(j) '= Mtrack_PA ;']);%%%
        %end
        [lub] = trackmem(Mtrack_PA,maxd,2,2,0);
        %figure, 
        %subplot(1,3,1), hold on, imshow(a,[]), plot(aM_PA(:,1),aM_PA(:,2),'+r','MarkerSize',8)
        %subplot(1,3,2), hold on, imshow(b,[]), plot(bM_PA(:,1),bM_PA(:,2),'+r','MarkerSize',8) %fig beads detect 2%
        %subplot(1,3,3), hold on, imshow(c,[]), plot(cM_PA(:,1),cM_PA(:,2),'+r','MarkerSize',8) %fig beads detect 2%
    else %MC 11juin2014
        lub=-1;%[]; %MC 11juin2014
    end %MC 11juin2014
else
    lub = -1;
end

%         Mtrack_PA = [];
%         for jj = 1:30
%             %jj
%             M_PA = image_info{jj, 3};
%             
%             Mtrack_PA = [Mtrack_PA; M_PA(:,1:2),jj*ones(size(M_PA,1),1),jj*ones(size(M_PA,1),1)];
%         %Mtrack_PA=[aM_PA(:,1:2),ones(size(aM_PA,1),1),ones(size(aM_PA,1),1); bM_PA(:,1:2),2*ones(size(bM_PA,1),1),2*ones(size(bM_PA,1),1);...
%         %    cM_PA(:,1:2),3*ones(size(cM_PA,1),1),3*ones(size(cM_PA,1),1)];
%         %                     var_Mtrack=genvarname(['Mtrack_' num2str(ff,'%03d') '_' num2str(i) num2str(j)]);
%         %                     eval(['Mtrack_' num2str(ff,'%02d') '_' num2str(i) num2str(j) '= Mtrack_PA ;']);%%%
%         end
% [lub] = trackmem(Mtrack_PA,maxd,2,2,0);
%end
%size(aM_PA)
%size(bM_PA)
%size(lub)
if lub~=-1
    ind = lub(:,5);
    im = lub(:,4);

    x = lub(:,1);
    y = lub(:,2);

    xplot = zeros(max(im), max(ind));
    yplot = zeros(max(im), max(ind));
    %size(xplot)
    for i = 1:length(ind)
        xplot(im(i), ind(i)) = x(i);
        yplot(im(i), ind(i)) = y(i);
    end
    xplot(xplot == 0) = NaN;
    yplot(yplot == 0) = NaN;
else
    xplot = nan(2,1);
    yplot = nan(2,1);
end

end