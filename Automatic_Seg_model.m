% Define the path for each folder used in your project
Input_folder='D:\...\Input';
Input_folder_ph='D:\...\Images';
Current_folder='D:\...\Model';
Output_folder='D:\...\Output';
%The model will analyze each file within the Input_folder having the
%specified extension e.g. jpg
cd(Input_folder)
file_images=dir('*.jpg');
cd(Input_folder_ph)
file_images_ph=dir('*.jpg');
for n=16:size(file_images,1)
    cd(Input_folder)
    filename_img=file_images(n).name;
    [filepath, name, ext]=fileparts(filename_img);
    picture=imread(filename_img);
    
    %Transformation in a monochromatic image.
    img=rgb2gray(picture);
    board_new=img;
    %----------------------------------------------------------------------
    %enhancement of the contrast of the image. Increase the highest intensity
    %and decrease of the lowest ones.
    img_adj=imadjust(img,[100 120]/255, [0 255]/255);%
    figure
    imshow(img_adj);    
    cd(Input_folder)
    N=size(img_adj);
    n_row=N(1);
    n_col=N(2);
    
    %----------------------------------------------------------------------
    %application of the 2-D first derivative Sobel operator to the smoothed image
    %to highlight regions of the image with high first spatial derivatives.
    % This operator performs the spatial gradient evaluation in both x and y directions
    % using two 3x3 convolution.Since the edges give rise to high local
    %image contrast, peaks in the image gradient along the edges are expected
    
    s=[1 2 1;0 0 0; -1 -2 -1];
    
    %%convlolution of the 2 matrices.
    %The convolution process works as follows: for each pixel of the grayscale image the
    % code computes the weighted average of the values of the pixels surrounding it
    % using as weights the values of a square matrix that is centred on the pixel of
    % interest.
    %The resulting image represents its gradient in the vertical direction.
    
    img_gradient=conv2(double(img_adj),double(s), 'same');
    
    %%Initialization of the matrices that contain the positive and the
    %%negative gradient of the image intensity.
    
    img_gradient_pos= zeros(n_row, n_col);
    img_gradient_neg= zeros(n_row, n_col);
    board_up_edge= zeros(n_row,n_col);
    board_down_edge= zeros(n_row,n_col);
    board_left_edge= zeros(n_row,n_col);
    board_right_edge= zeros(n_row,n_col);
    
    %%Initialization of the vectors that contain the maximum value of
    %%gradient and its position (row) along each column of the image.
    
    max_val= zeros(n_col,1);
    pos_idx=zeros(n_col,1);
    row=zeros(n_col,1);
    col=zeros(n_col,1);
    
    %all the pixel values lower than zero are neglected. Only the positive
    %gradient is considered since we want to highlight the upper edge of
    %the white board.
    for j=1:n_col
        for i=1:n_row
            if img_gradient(i,j)>0
                img_gradient_pos(i,j)=img_gradient(i,j);
            else
                img_gradient_pos(i,j)=0;
            end
        end
    end
    
    %the white clusters of pixels (connected components) made by less than 500
    %pixles are removed from the image.
    
    img_gradient_pos_new=bwareaopen(img_gradient_pos, 800);
    %----------------------------------------------------------------------
    %The resulting image represents the remaining clusters of white pixels with values equal to 1.
    %The image is multiplied by the image containing the positive gradient
    %values. This allows to characterize each white pixel with the positive gradient values.
    
    immult=immultiply(img_gradient_pos, img_gradient_pos_new);
    
    %Search for the the maximum positive gradient value in each column of the "immult" image.
    
    for j=1:n_col
        [max_val(j),pos_idx(j)] = max(immult(1:n_row,j));
        pos_idx(j) =pos_idx(j);
        [row(j,1), col(j,1)]=ind2sub(size(immult(:,j)), pos_idx(j));
    end
    
    %% Filtration of the pixels correspondent to the maximum gradient that don't belong to the white board edge.
    
    %Step 1: If none of the pixel in the column j has positive gradient
    %but along the column j-1 a pixel, characterized as belonging to the edge, has been identified,
    %the code selects as edge the pixel located at the same row of of the pixel correspondent to the column j-1
    
    for j=2:(n_col)
        if  size(immult(immult(:,j)==0,1),1)==size(immult,1)
            if j>1 & max_val(j-1)~=0 & row(j-1)< size(img_gradient,1)/2-300
                row(j)= row(j-1);
                max_val(j)=1;
            end
            row(j)= n_row;
            
            %Step 2: If at least one of the pixel has a positive gradient value and
            %if the pixel with the maximum positive gradient in the vertical direction is
            % located in the upper half of the image and if its distance from the edge
            % identified in column j-1 is less than 10 pixels, then this will be the next edge
            % pixel and the code will move to column j+1.
            
        elseif  max_val(j)>0 & max_val(j-1)>0 & row(j-1)< size(img_gradient,1)/2-300 & row(j)< size(img_gradient,1)/2-300
            if row(j)>row(j-1)+10
                row(j)=row(j-1);
            end
        elseif  max_val(j)>0 & max_val(j-1)>0 & row(j-1)< size(img_gradient,1)/2-300  & row(j)>= size(img_gradient,1)/2-300
            [max_val(j), pos_idx(j)]=max(immult(1:size(img_gradient_pos,1)/2-300 ,j));
            [row(j), col(j)]=ind2sub(size(immult(:,j)), pos_idx(j));
            if row(j)==1 | row(j)>10+row(j-1)
                row(j)=row(j-1);
                if max_val(j)==0
                    max_val(j)=max_val(j-1);
                end
            elseif row(j)> 1 & row(j)<=10+row(j-1)
                row(j)=row(j);
            end
            
        end
    end
    
    for j=1:n_col
        for i=1:n_row
            if i<row(j)
                board_up_edge(i,j)=0;
            else
                board_up_edge(i,j)=1;
            end
        end
    end
    
    
    %% same of the previous loop, but for defining the lower limit
    m_val= zeros(n_col,1);
    pos_idx_m= zeros(n_col,1);
    row_m=zeros(n_col,1);
    col_m=zeros(n_col,1);
    for j=1:n_col
        for i= 1:n_row
            if img_gradient(i,j)<0
                img_gradient_neg(i,j)=-img_gradient(i,j);
            else
                img_gradient_neg(i,j)=0;
            end
        end
    end
    
    img_gradient_neg_new=bwareaopen(img_gradient_neg, 800);
    immult=immultiply(double(img_gradient_neg), double(img_gradient_neg_new));
    immult(immult==0)=inf;
    for j=1:n_col
        [m_val(j), pos_idx_m(j)] = min(immult(1:(n_row),j));
        pos_idx_m(j)=pos_idx_m(j);
        [row_m(j), col_m(j)]=ind2sub(size(immult(:,j)), pos_idx_m(j));
    end
    
    for j=2:(n_col)
        if m_val(j)==inf & m_val(j-1)==inf
            row_m(j)=1;
            elseif m_val(j-1)>0 & m_val(j)>0 & row_m(j-1)> size(img_gradient,1)/2+300 & row_m(j)> size(img_gradient,1)/2+300
            if row_m(j)<row_m(j-1)-10
                row_m(j)=row_m(j-1);
            else
                row_m(j)=row_m(j);
            end
        elseif m_val(j-1)>0 & m_val(j)>0 & row_m(j)<= size(img_gradient,1)/2+300
            [m_val(j), pos_idx_m(j)]=min(immult(size(img_gradient_neg,1)/2-300:n_row ,j));
            pos_idx_m(j) =pos_idx_m(j);
            [row_m(j), col_m(j)]=ind2sub(size(immult(:,j)), pos_idx_m(j));
            if row_m(j)==n_row | row_m(j)<row_m(j-1)-10
                row_m(j)=row_m(j-1);
                if m_val(j)==inf
                    m_val(j)=m_val(j-1);
                end
            else
                row_m(j)=row_m(j);
            end
        end
    end
   
    for j=1:n_col
        for i=1:n_row
            if (i>row_m(j))
                board_down_edge(i,j)=0;
            else
                board_down_edge(i,j)=1;
            end
        end
    end
    
    
    %% same algorithm for defining the lateral limits of the board
    
    s1=s.';
    
    img_gradient_y=conv2(double(img_adj),double(s1),'same');
    img_gradient_y_pos= zeros(n_row, n_col);
    img_gradient_y_neg= zeros(n_row, n_col);
    max_val_y=zeros(n_row,1);
    pos_idx_y=zeros(n_row,1);
    row_y=zeros(n_row,1);
    col_y=zeros(n_row,1);
    
    for i=1:n_row
        for j=1 : n_col
            if img_gradient_y(i,j)>0
                img_gradient_y_pos(i,j)=img_gradient_y(i,j);
            else
                img_gradient_y_pos(i,j)=0;
            end
        end
    end
    
    img_gradient_y_pos_new=bwareaopen(img_gradient_y_pos, 500);
    immult=immultiply(img_gradient_y_pos, img_gradient_y_pos_new);
    
    for i=1:n_row
        [max_val_y(i), pos_idx_y(i)] = max((immult(i, 1:n_col)));
        
        [row_y(i), col_y(i)]=ind2sub(size(immult(i,:)),pos_idx_y(i));
    end
    for i= 2:(n_row)
        if max_val_y(i)==0
            col_y(i)=size(img_gradient_y,2);
        elseif max_val_y(i-1)>0 & max_val_y(i)>0 & col_y(i-1)< size(img_gradient_y,2)/2-300 & col_y(i)< size(img_gradient_y,2)/2-300
            if col_y(i)>10+col_y(i-1)
                col_y(i)=col_y(i-1);
            else col_y(i)=col_y(i);
            end
        elseif  max_val_y(i-1)>0 & max_val_y(i)>0 & col_y(i)>= size(img_gradient,2)/2-300
            [max_val_y(i), pos_idx_y(i)]=max(immult(i,1:size(img_gradient_pos,2)/2-300));
            pos_idx_y(i) =pos_idx_y(i);
            [row_y(i), col_y(i)]=ind2sub(size(immult(i,:)), pos_idx_y(i));
            if col_y(i)==1 | col_y(i)>10+col_y(i-1)
                col_y(i)=col_y(i-1);
                if max_val_y(i)==0
                    max_val_y(i)=max_val_y(i-1);
                end
            else col_y(i)=col_y(i);
            end
        end
    end
    
    for i=1:n_row
        for j=1:n_col
            if (j<col_y(i))
                board_left_edge(i,j)=0;
            else
                board_left_edge(i,j)=1;
            end
        end
    end
    
    m_val_my=zeros(n_row,1);
    pos_idx_my=zeros(n_row,1);
    row_my=zeros(n_row,1);
    col_my=zeros(n_row,1);
    for i=1:n_row
        for j= 1:n_col
            if img_gradient_y(i,j)<0
                img_gradient_y_neg(i,j)=-img_gradient_y(i,j);
            else
                img_gradient_y_neg(i,j)=0;
            end
        end
    end
    img_gradient_y_neg_new=bwareaopen(img_gradient_y_neg, 500);
    immult=immultiply(img_gradient_y_neg, img_gradient_y_neg_new);
    immult(immult<=0)=inf;
    for i=1:n_row
        [m_val_my(i), pos_idx_my(i)] = min(immult(i, 1:n_col));
        pos_idx_my(i)=pos_idx_my(i);
        [row_my(i), col_my(i)]=ind2sub(size(immult(i,:)), pos_idx_my(i));
    end
    
    for i= 2:(n_row)
        if m_val_my(i)==inf 
            col_my(i)=1;
            elseif m_val_my(i-1)>0 & m_val_my(i)>0 & col_my(i-1)> size(img_gradient_y,2)/2-300 & col_my(i)> size(img_gradient_y,2)/2-300
            if col_my(i)<col_my(i-1)-10
                col_my(i)=col_my(i-1);
            else col_my(i)=col_my(i);
            end
        elseif m_val_my(i-1)>0 & m_val_my(i)>0 & col_y(i)<= size(img_gradient,2)/2-300
            [m_val_my(i), pos_idx_my(i)]=min(immult(i,size(img_gradient_pos,2)/2-300:n_col));
            pos_idx_my(i) =pos_idx_my(i);
            [row_my(i), col_my(i)]=ind2sub(size(immult(i,:)), pos_idx_my(i));
            if col_my(i)==n_col | col_my(i)<col_my(i-1)-10
                col_my(i)=col_my(i-1);
                if m_val_my(i)==inf
                    m_val_my(i)=m_val_my(i-1);
                end
            else col_my(i)=col_my(i);
            end
        end
    end
    
    for i=1:n_row
        for j= 1:n_col
            if (j>col_my(i))
                board_right_edge(i,j)=0;
            else
                board_right_edge(i,j)=1;
            end
        end
    end
    
    
    %% definition of a new matrix for representing the white board
    board=zeros(n_row, n_col);
    
    for  i=1:n_row
        for j=1:n_col
            if board_up_edge(i,j)==0 | board_down_edge(i,j)==0 | board_left_edge(i,j)==0 | board_right_edge(i,j)==0
                board(i,j)=0;
            else board(i,j)=1;
            end
        end
    end
    
    c_board=bwlabel(board,8);
    c_board_bounding=regionprops(c_board, 'Area');
    c_board_bounding=struct2cell(c_board_bounding)';
    c_board_bounding=cell2mat(c_board_bounding);
    
    AA=zeros(size(board));
    for i=1:size(board,1)
        for j=1:size(board,2)
            for r=1:size(c_board_bounding,1)
                if c_board(i,j)==r
                    AA(i,j)=c_board_bounding(r,1);
                end
            end
            
        end
    end
    for i=1:size(board,1)
        for j=1:size(board,2)
            if AA(i,j)~= max(c_board_bounding)
                AA(i,j)=0;
            end
        end
    end
    board=AA;
    
    %% new image representing only the board ith the pebbles on the top
    
    for i=1:size(board,1)
        for j=1:size(board,2)
            if board(i,j)==1 && img(i,j)>130 || board(i,j)==0 && img(i,j)>130
                board(i,j)=1;
            elseif board(i,j)==0 && img(i,j)<100
                board(i,j)=0;
            elseif board(i,j)==1 && img(i,j)<100
                board(i,j)=0;
            end
        end
    end
    choice=menu('Choose', 'accept', 'reject');
    if choice==1
        board=bwareaopen(board, 50000);
        board=imfill(board, 'holes');
        for i=1:size(board,1)
            for j=1:size(board,2)
                if board(i,j)==1;
                    img(i,j)=img(i,j);
                    picture(i,j,:)=picture(i,j,:);
                else
                    img(i,j)=0;
                    picture(i,j,:)=[180,180,180];
                end
            end
        end
        figure
        imshow(picture)
        img_mono=rgb2gray(picture);
        img_mono_bin=imbinarize(img_mono);
        board_adj=imadjust(img,[60 200]/255, [0 1],1.2);
        board_bin=imbinarize(board_adj);
        figure
        imshow(board_bin)
        figure
        imshow(img_mono_bin)
        board_bin_opposite=1-board_bin;
    end
    
    multiply=immultiply(board_bin_opposite, board); 
    board_new=bwareaopen(multiply, 500);
    
    figure
    imshow(board_new)
    f=1;
    cc=bwlabel(board_new, 8);
    cc_rgb=label2rgb(cc);
    props=regionprops(cc, 'Area');
    props=struct2cell(props);
    props=cell2mat(props);
    choice=menu('Choose', 'cc_ok', 'cc_wrong');
    while choice==2
        board_new=imerode(board_new, strel('disk',2,8));
        figure
        imshow(board_new);
        g=0;
        choice=menu('Choose', 'accept', 'reject');
    end
    cc_new=bwlabel(board_new, 8);
    cc_rgb_new=label2rgb(cc_new);
    %
    BoundingBox=[];
    centroids=[];
    P=[];
    P=regionprops(cc_new, 'centroid', 'BoundingBox');
    centroids = cat(1,P.Centroid);
    for q=1:size(P, 1)
        BoundingBox(q, :)=P(q).BoundingBox;
    end
    centroids=centroids(BoundingBox(:,3)>30 & BoundingBox(:,3)<2000 & BoundingBox(:,4)>30 & BoundingBox(:,4)<2000, :);
    BoundingBox=BoundingBox(BoundingBox(:,3)>30 & BoundingBox(:,3)<2000 & BoundingBox(:,4)>30 & BoundingBox(:,4)<2000, :);
    cd 'C:\Users\agu005\Desktop\Master Thesis Gurini\Image_Processing'
    figure
    imshow(cc_new)
    hold on
    plot(centroids(:,1),centroids(:,2),'b*')
    hold off
    print('Centroidi', '-dtiff', '-r1000');
    cd(Input_folder_ph)
    filename_img_ph=file_images_ph(n).name;
    [filepath_ph,name_ph, ext_ph]=fileparts(filename_img_ph);
    im_ph=imread(filename_img_ph);
    im_ph_reverse=~rgb2gray(im_ph);
    im_ph_filled=imfill(im_ph_reverse, 'holes');
    cc_ph=bwlabel(im_ph_filled, 8);
    props_ph=regionprops(cc_ph, 'Centroid', 'BoundingBox');
    cent_ph=[];
    BB=[];
    cent_ph=cat(1, props_ph.Centroid);
    figure
    imshow(im_ph)
    hold on
    plot(cent_ph(:,1),cent_ph(:,2),'b*')
    hold off
    for q=1:size(props_ph, 1)
        BB(q, :)=props_ph(q).BoundingBox;
    end
    cent_ph=cent_ph(BB(:,3)>30 & BB(:,4)>30, :);
    BB=BB(BB(:,3)>30 & BB(:,4)>30, :);
    
    for w =31:size(centroids, 1)
        r_cent=centroids(w,2);
        c_cent=centroids(w,1);
        rect_1(1, :)=BoundingBox(w,:);
        if rect_1(1,1)>rect_1(1,3)*0.3 & rect_1(1,2)>rect_1(1,4)*0.3 & rect_1(1,1)<size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)<size(board,1)-rect_1(1,4)*1.6
            rect(1,1)=rect_1(1,1)-rect_1(1,3)*0.3;
            rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
            
            rect(1,3)=rect_1(1,3)*1.6;
            rect(1,4)=rect_1(1,4)*1.6;
        elseif rect_1(1,1)<rect_1(1,3)*0.3 & rect_1(1,2)>rect_1(1,4)*0.3 & rect_1(1,1)<size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)<size(board,1)-rect_1(1,4)*1.6
            if rect_1(1,1)>10
                rect(1,1)=rect_1(1,1)-10;
                rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
                rect(1,3)=rect_1(1,3)*1.3+10;
                rect(1,4)=rect_1(1,4)*1.6;
            else
                rect(1,1)=0.5;
                rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
                rect(1,3)=rect_1(1,3)*1.3+rect_1(1,1)-0.5;
                rect(1,4)=rect_1(1,4)*1.6;
            end
        elseif rect_1(1,1)<rect_1(1,3)*0.3 & rect_1(1,2)>rect_1(1,4)*0.3 & rect_1(1,1)<size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)>size(board,1)-rect_1(1,4)*1.6
            if rect_1(1,1)>10
                rect(1,1)=rect_1(1,1)-10;
                rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
                rect(1,3)=rect_1(1,3)*1.3+10;
                rect(1,4)=rect_1(1,4)*1.6;
            else
                rect(1,1)=0.5;
                rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
                rect(1,3)=rect_1(1,3)*1.3+rect_1(1,1)-0.5;
                rect(1,4)=rect_1(1,4)*1.6;
            end
        elseif rect_1(1,1)>rect_1(1,3)*0.3 & rect_1(1,2)<rect_1(1,4)*0.3 & rect_1(1,1)<size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)<size(board,1)-rect_1(1,4)*1.6
            if rect_1(1,2)>10
                rect(1,1)=rect_1(1,1)-rect_1(1,3)*0.3;
                rect(1,2)=rect_1(1,2)-10;
                rect(1,3)=rect_1(1,3)*1.6;
                rect(1,4)=rect_1(1,4)*1.3+10;
            else
                rect(1,1)=rect_1(1,1)-rect_1(1,3)*0.3;
                rect(1,2)=0.5;
                rect(1,3)=rect_1(1,3)*1.6;
                rect(1,4)=rect_1(1,4)*1.3+rect_1(1,2)-0.5;
            end
            
        elseif rect_1(1,1)<rect_1(1,3)*0.3 & rect_1(1,2)<rect_1(1,4)*0.3 & rect_1(1,1)<size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)<size(board,1)-rect_1(1,4)*1.6
            if rect_1(1,1)>10 & rect_1(1,1)>10
                rect(1,1)=rect_1(1,1)-10;
                rect(1,2)=rect_1(1,2)-10;
                
                rect(1,3)=rect_1(1,3)*1.3+10;
                rect(1,4)=rect_1(1,4)*1.3+10;
            else
                rect(1,1)=0.5;
                rect(1,2)=0.5;
                
                rect(1,3)=rect_1(1,3)*1.3+rect_1(1,1)-0.5;
                rect(1,4)=rect_1(1,4)*1.3+rect_1(1,1)-0.5;
            end
        elseif rect_1(1,1)>rect_1(1,3)*0.3 & rect_1(1,2)>rect_1(1,4)*0.3 & rect_1(1,1)<size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)>size(board,1)-rect_1(1,4)*1.6
            rect(1,1)=rect_1(1,1)-rect_1(1,3)*0.3;
            rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
            
            rect(1,3)=rect_1(1,3)*1.6;
            rect(1,4)=size(board,1)-rect(1,2);
        elseif rect_1(1,1)>rect_1(1,3)*0.3 & rect_1(1,2)>rect_1(1,4)*0.3 & rect_1(1,1)>size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)<size(board,1)-rect_1(1,4)*1.6
            rect(1,1)=rect_1(1,1)-rect_1(1,3)*0.3;
            rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
            
            rect(1,3)=size(board,2)-rect(1,1);
            rect(1,4)=rect_1(1,4)*1.6;
        elseif rect_1(1,1)>rect_1(1,3)*0.3 & rect_1(1,2)>rect_1(1,4)*0.3 & rect_1(1,1)>size(board,2)-rect_1(1,3)*1.6 & rect_1(1,2)>size(board,1)-rect_1(1,4)*1.6
            rect(1,1)=rect_1(1,1)-rect_1(1,3)*0.3;
            rect(1,2)=rect_1(1,2)-rect_1(1,4)*0.3;
            
            rect(1,3)=size(board,2)-rect(1,1);
            rect(1,4)=size(board,1)-rect(1,2);
        end
        aa=imcrop(picture, rect);
        imshow(aa);
        choice=menu('Choose', 'accept', 'reject');
        if choice==1
            aa_mask=imcrop(board,rect);
            
            cd(Current_folder)
            [Shape,  threshOut]=Segmentation_function(aa, aa_mask, rect, w); 
            
            choice=menu('Choose', 'accept', 'reject');
            if choice==1
                sample{w,1}=Shape;
                dist_cent=[];
                for q=1:size(cent_ph,1)
                    dist_cent(q,1)=sqrt((r_cent-cent_ph(q,2))^2+(c_cent-cent_ph(q,1))^2);
                end
                cd(Output_folder)
                imwrite(1-imfill(Shape, 'holes'), strcat(name, '.', num2str(w), '.png'))
                BB_new=BB(dist_cent(:,1)==min(dist_cent(:,1)), :);
                rect_ph=[BB_new(1)-40, BB_new(2)-40, BB_new(3)+80, BB_new(4)+80];
                ph=imcrop(rgb2gray(im_ph), rect_ph);
                ph_bin=imbinarize(ph);
                ph_bin_rev=1-ph_bin;
                con_comp_ph=bwlabel(ph_bin_rev);
                con_comp_ph_prop=regionprops(con_comp_ph, 'Area');
                con_comp_ph_prop=struct2cell(con_comp_ph_prop);
                con_comp_ph_prop=cell2mat(con_comp_ph_prop);
                if size(con_comp_ph_prop,2)>=2
                    area_con=zeros(size(ph,1), size(ph,2));
                    for i=1:size(ph,1)
                        for j=1:size(ph,2)
                            for p=1:size(con_comp_ph_prop,2)
                                if con_comp_ph(i,j)==p
                                    area_con(i,j)=con_comp_ph_prop(1,p);
                                end
                            end
                        end
                    end
                    for i=1:size(ph,1)
                        for j=1:size(ph,2)
                            if area_con(i,j)~=max(con_comp_ph_prop)
                                ph_bin(i,j, :)=[1];
                            end
                        end
                    end
                end
                cd(Output_folder)
                imwrite(ph_bin, strcat(name_ph, '.',num2str(w), '.png'))
                BB_new=[];
                
                close all
            end
        end
        
        w=w+1;
    end
    
    close all
    
end
