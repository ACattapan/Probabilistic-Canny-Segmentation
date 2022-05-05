function [r_e, c_e, r_c, c_c, BW, re_1, re_2, ce_1, ce_2, angle, threshOut]= Endpoints(img,img2)

%import image and application of the morphological operators for
%reducing the internal edges linked to the intra-granular color and texture
%variations
gray=rgb2gray(img);
gray=imgaussfilt(gray, 0.8);
N=size(gray);
n_row=N(1);
n_col=N(2);
se = strel('disk',15,8);
Ge = imerode(gray,se);
Gobr = imreconstruct(Ge,gray);
se1=strel('disk',10,8);
Gobrd = imdilate(Gobr,se1);
Gobrcbr = imreconstruct(imcomplement(Gobrd),imcomplement(Gobr));
Gobrcbr = imcomplement(Gobrcbr);
%application of the Canny edge detector without imposing any parameters
[E, threshOut]=edge(Gobrcbr,'Canny');
%identification of the branchpoints
E=bwmorph(E, 'skel');
BR=bwmorph(E,'branchpoints');
%break of the branches of the connected component
for i=1:n_row
    for j=1:n_col
        if BR(i,j)==1
            E(i,j)=0;
        end
    end
end

cc=bwlabel(E,8);
cc_rgb=label2rgb(cc);
%classification of the connected components based on the area (so the
%number of pixels that compose the connected component)
props=regionprops(cc(:,:), 'Area');
props=struct2cell(props);
props=cell2mat(props);
N_props=size(props);
n_p=N_props(2);

Area=zeros(n_row,n_col);
for i=1:n_row
    for j=1:n_col
        for r=1:n_p
            if cc(i,j)==r
                Area(i,j)=props(1,r);
            end
        end
    end
end
%Area=imdilate(Area, strel('disk',2,8));
E2_rgb=label2rgb(Area);
E2_GRAY= rgb2gray(E2_rgb);

Area_max=zeros(n_row,n_col);
for i=1:n_row
    for j=1:n_col
        if Area(i,j)==max(props(:,:))
            Area_max(i,j)=Area(i,j);
        else
            Area_max(i,j)=0;
        end
    end
end
% Area_max_dil=imdilate(Area_max,strel('disk',2,8));
BW=imbinarize(Area_max);
idx=find(BW~=0);
[row,col]=ind2sub(N,idx);
figure
imshow(BW)

%endpoint of the longest component
s = regionprops(BW,'centroid');
centroids = cat(1,s.Centroid);

c_c=round(centroids(1));
r_c=round(centroids(2));

endpoints=bwmorph(BW, 'endpoint');
endpoints=find(endpoints~=0);
%Result from the application of Canny algorithm

[r_e, c_e]=ind2sub(N, endpoints);
h=1;
figure

a=imshow(gray);
hold on
a=imshow(E2_rgb);
hold off
alpha(a, 255-E2_GRAY);
%print('Canny_risultato','-dtiff', '-r1000')%, '-djpg', '-r100)
choice=menu('Choose', 'Long_con_pebble', 'Long_con_board', 'completely wrong');
if choice==3
    r_e=[];
    c_e=[];
    re_1=[];
    re_2=[];
    ce_1=[];
    ce_2=[];
    
elseif choice==2
    while  h==1
        for i=1:size(r_e,1)
            if r_e(i)-0.5>4 & c_e(i)-0.5>4 & size(img2,1)-r_e(i)>4 & size(img2,2)-c_e(i)>4
                Area_new=Area-Area_max;
                Area_new(Area_new~=max(max(Area_new)))=0;
                Area_max=Area_new;
                BW=imbinarize(Area_new);
                endpoints=bwmorph(BW, 'endpoint');
                endpoints=find(endpoints~=0);
                [r_e,c_e]=ind2sub(N, endpoints);
                figure
                imshow(BW)
                choice=menu('Choose', 'Long_con_pebble', 'Long_con_board');
                if choice==2
                    i=1;
                    h=1;
                    Area_max=Area_max+Area_new;
                else h=0;
                end
                
            elseif r_e(i)-0.5<4 & c_e(i)-0.5>4 &  size(img2,1)-r_e(i)>4 & size(img2,2)-c_e(i)>4
                Area_new=Area-Area_max;
                Area_new(Area_new~=max(max(Area_new)))=0;
                BW=imbinarize(Area_new);
                endpoints=bwmorph(BW, 'endpoint');
                endpoints=find(endpoints~=0);
                [r_e, c_e]=ind2sub(N, endpoints);
                figure
                imshow(BW)
                choice=menu('Choose', 'Long_con_pebble', 'Long_con_board');
                if choice==2
                    i=1;
                    h=1;
                    Area_max=Area_max+Area_new;
                else h=0; break
                end
                
            elseif r_e(i)-0.5>4 & c_e(i)-0.5<4 & size(img2,1)-r_e(i)>4 & size(img2,2)-c_e(i)>4
                
                Area_new=Area-Area_max;
                Area_new(Area_new~=max(max(Area_new)))=0;
                BW=imbinarize(Area_new);
                endpoints=bwmorph(BW, 'endpoint');
                endpoints=find(endpoints~=0);
                [r_e,c_e]=ind2sub(N, endpoints);
                figure
                imshow(BW)
                choice=menu('Choose', 'Long_con_pebble', 'Long_con_board');
                if choice==2
                    i=1;
                    h=1;
                    Area_max=Area_max+Area_new;
                else h=0; break
                end
                
            elseif r_e(i)-0.5>4 & c_e(i)-0.5>4 & size(img2,1)-r_e(i)<4 & size(img2,2)-c_e(i)>4
                
                Area_new=Area-Area_max;
                Area_new(Area_new~=max(max(Area_new)))=0;
                BW=imbinarize(Area_new);
                endpoints=bwmorph(BW, 'endpoint');
                endpoints=find(endpoints~=0);
                [r_e,c_e]=ind2sub(N, endpoints);
                figure
                imshow(BW)
                choice=menu('Choose', 'Long_con_pebble', 'Long_con_board');
                if choice==2
                    i=1;
                    h=1;
                    Area_max=Area_max+Area_new;
                else h=0; break
                end
                
            elseif r_e(i)-0.5>4 & c_e(i)-0.5>4 & size(img2,1)-r_e(i)<4 & size(img2,2)-c_e(i)<4
                
                Area_new=Area-Area_max;
                Area_new(Area_new~=max(max(Area_new)))=0;
                BW=imbinarize(Area_new);
                endpoints=bwmorph(BW, 'endpoint');
                endpoints=find(endpoints~=0);
                [r_e,c_e]=ind2sub(N, endpoints);
                figure
                imshow(BW)
                choice=menu('Choose', 'Long_con_pebble', 'Long_con_board');
                if choice==2
                    i=1;
                    h=1;
                    Area_max=Area_max+Area_new;
                else h=0; break
                end
                
            end
        end
    end
end
if ~isempty(endpoints)
    Br=bwmorph(BW, 'branchpoints');
    idx=find(Br~=0);
    [r_b, c_b]=ind2sub(N, idx);
    d_e_b=zeros(size(r_e,1),1);
    check=0;
    skel= bwmorph(BW,'skel',Inf);
    
    skel= bwmorph(BW,'skel',Inf);
    
    Dmask = false(size(skel));
    if ~isempty(r_e)
        for k = 1:numel(c_e)
            D = bwdistgeodesic(skel,c_e(k),r_e(k));
            
        end
        D=D+BW;
    end
    o=1;
    if ~isempty(r_b) & ~isempty(r_e) & size(r_e,1)>=2
        
        x_b=r_b(o);
        y_b=c_b(o);
        check=0;
        while check==0
            x_b=r_b(o);
            y_b=c_b(o);
            check=0;
            matrix=[];
            points=[];
            d=[];
            skel= bwmorph(BW,'skel',Inf);
            Dmask = false(size(skel));
            for k = 1:numel(c_e)
                D = bwdistgeodesic(skel,c_e(k),r_e(k));
            end
            D=D+BW;
            for j=1:size(r_e,1)
                if ~isempty(x_b)
                    d_e_b(j,1)=sqrt((r_e(j)-x_b)^2+(c_e(j)-y_b)^2);
                    j=j+1;
                else check=1;
                end
                
            end
            if any(d_e_b<sqrt(32))
                BW(r_e(d_e_b(:,1)<sqrt(32),1), c_e(d_e_b(:,1)<sqrt(32),1))=0;
                props=regionprops(BW, 'Area');
                BW=bwareaopen(BW, max(cell2mat(struct2cell(props)))-1);
                endpoints=bwmorph(BW, 'endpoint');
                endpoints=find(endpoints~=0);
                [r_e,c_e]=ind2sub(N, endpoints);
                Br=bwmorph(BW, 'branchpoints');
                idx=find(Br~=0);
                [c_b, r_b]=ind2sub(N, idx);
            else
                r=1;
                for n=x_b-4:x_b+4
                    for m=y_b-4:y_b+4
                        if n==x_b-4 & ~isnan(D(n,m)) | n==x_b+4 & ~isnan(D(n,m))| m==y_b+4 & ~isnan(D(n,m)) |  m==y_b-4 & ~isnan(D(n,m))
                            d(r,1)=abs(D(n,m)-D(x_b, y_b));
                            r=r+1;
                        end
                    end
                end
                d_min=min(d);
                
                mat=find(D==(D(x_b, y_b)+d_min) | D==(D(x_b, y_b)-d_min));
                mat(:,2)= D(mat(:,1));
                [mat(:,3), mat(:,4)]=ind2sub(N, mat(:,1));
                l=1;
                for h=1:size(mat,1)
                    if mat(h,3)<=x_b+4 & mat(h,4)<=y_b+4 & mat(h,3)>=x_b-4 & mat(h,4)>=y_b-4
                        matrix(l, :)=mat(h, :);
                        l=l+1;
                    end
                end
                for t=1:size(matrix,1)
                    if matrix(t,3)>=x_b && matrix(t,4)>y_b
                        matrix(t,5)=(atan((matrix(t,3)-x_b)/(matrix(t,4)-y_b))*180/pi)+270;
                    elseif matrix(t,3)>=x_b && matrix(t,4)<y_b
                        matrix(t,5)=(-atan((matrix(t,3)-x_b)/(matrix(t,4)-y_b))*180/pi)+180;
                    elseif matrix(t,3)<x_b && matrix(t,4)<y_b
                        matrix(t,5)=(-atan((matrix(t,3)-x_b)/(matrix(t,4)-y_b))*180/pi)+180;
                    elseif matrix(t,3)<x_b && matrix(t,4)>=y_b
                        matrix(t,5)=(-atan((matrix(t,3)-x_b)/(matrix(t,4)-y_b))*180/pi);
                        
                    end
                end
                for i=1:size(matrix,1)
                    matrix(i,6)=(matrix(t,3)-x_b)/(matrix(t,4)-y_b);
                end
                if size(matrix(matrix(:,2)>D(x_b,y_b), 2),1)>=2
                    if size(matrix(matrix(:,2)<D(x_b,y_b)),1)>=2
                        
                        for i=1:size(matrix,1)
                            matrix(i,7)=sqrt((matrix(i,3)-x_b)^2+(matrix(i,4)-y_b)^2);
                        end
                        m_comp=matrix(matrix(:,2)<D(x_b,y_b) & matrix(:,7)==min(matrix(matrix(:,2)<D(x_b,y_b),7)), 6);
                        
                        for i=1:size(matrix,1)
                            matrix(i,8)=abs(matrix(i,6))-abs(m_comp);
                            matrix(i,9)=abs(matrix(i,5)-matrix(matrix(:,2)<D(x_b,y_b)& matrix(:,7)==min(matrix(matrix(:,2)<D(x_b,y_b),7)), 5));
                        end
                    else
                        for i=1:size(matrix,1)
                            matrix(i,7)=sqrt((matrix(i,3)-x_b)^2+(matrix(i,4)-y_b)^2);
                        end
                        for i=1:size(matrix,1)
                            matrix(i,8)=abs(matrix(i,6))-abs(matrix(matrix(:,2)<D(x_b,y_b),6));
                            matrix(i,9)=abs(matrix(i,5)-matrix(matrix(:,2)<D(x_b,y_b), 5));
                        end
                        
                    end
                    points=matrix(matrix(:,2)>D(x_b,y_b), :);
                    
                    BW(points(points(:,8)~=min(points(:,8)) | abs(points(:,9)-180)~=min(abs(points(:,9)-180)),3), points(points(:,8)~=min(points(:,8)) | abs(points(:,9)-180)~=min(abs(points(:,9)-180)),4))=0;
                    
                    props=regionprops(BW, 'Area');
                    BW=bwareaopen(BW, max(cell2mat(struct2cell(props)))-1);
                    endpoints=bwmorph(BW, 'endpoint');
                    endpoints=find(endpoints~=0);
                    [r_e,c_e]=ind2sub(N, endpoints);
                    Br=bwmorph(BW, 'branchpoints');
                    idx=find(Br~=0);
                    [r_b, c_b]=ind2sub(N, idx);
                    
                elseif size(matrix(matrix(:,2)<D(x_b,y_b), 2),1)>=2
                    if size(matrix(matrix(:,2)>D(x_b,y_b)),1)>=2
                        
                        for i=1:size(matrix,1)
                            matrix(i,7)=sqrt((matrix(i,3)-x_b)^2+(matrix(i,4)-y_b)^2);
                        end
                        m_comp=matrix(matrix(:,2)>D(x_b,y_b) & matrix(:,7)==min(matrix(matrix(:,2)>D(x_b,y_b),7)), 6);
                        
                        for i=1:size(matrix,1)
                            matrix(i,8)=abs(matrix(i,6))-abs(m_comp);
                            matrix(i,9)=abs(matrix(i,5)-matrix(matrix(:,2)>D(x_b,y_b)& matrix(:,7)==min(matrix(matrix(:,2)>D(x_b,y_b),7)), 5));
                        end
                    else
                        for i=1:size(matrix,1)
                            matrix(i,7)=sqrt((matrix(i,3)-x_b)^2+(matrix(i,4)-y_b)^2);
                        end
                        for i=1:size(matrix,1)
                            matrix(i,8)=abs(matrix(i,6))-abs(matrix(matrix(:,2)>D(x_b,y_b),6));
                            matrix(i,9)=abs(matrix(i,5)-matrix(matrix(:,2)>D(x_b,y_b), 5));
                        end
                    end
                    points=matrix(matrix(:,2)<D(x_b,y_b), :);
                    
                    BW(points(points(:,8)~=min(points(:,8)) & abs(points(:,9)-180)~=min(abs(points(:,9)-180)),3), points(points(:,8)~=min(points(:,8)) & abs(points(:,9)-180)~=min(abs(points(:,9)-180)),4))=0;
                    
                    props=regionprops(BW, 'Area');
                    BW=bwareaopen(BW, max(cell2mat(struct2cell(props)))-1);
                    endpoints=bwmorph(BW, 'endpoint');
                    endpoints=find(endpoints~=0);
                    [r_e,c_e]=ind2sub(N, endpoints);
                    Br=bwmorph(BW, 'branchpoints');
                    idx=find(Br~=0);
                    [r_b, c_b]=ind2sub(N, idx);
                end
            end
            if size(r_e, 1)>2 | ~isempty(r_b)
                Br=bwmorph(BW, 'branchpoints');
                idx=find(Br~=0);
                [r_b, c_b]=ind2sub(N, idx);
                if any(r_b==x_b) && any(c_b==y_b)
                    for j=1:size(r_e,1)
                        d_e_b(j,1)=sqrt((r_e(j)-x_b)^2+(c_e(j)-y_b)^2);
                        j=j+1;
                        
                    end
                    
                    if any(d_e_b(:,1)<=sqrt(32))
                        check=0;
                        d_e_b=[];
                        mat=[];
                        matrix=[];
                        
                    else
                        
                        check=0;
                        o=o;
                        
                        
                    end
                else
                    
                    check=0;
                    o=1;
                    
                end
            else
                check=1;
                
                
            end
        end
    elseif  ~isempty(r_b) & ~isempty(r_e)
        if D(r_e,c_e)==1 | D(r_e,c_e)==max(max(D))
            x_b=r_b(o);
            y_b=c_b(o);
            check=0;
            while check==0
                x_b=r_b(o);
                y_b=c_b(o);
                check=0;
                matrix=[];
                points=[];
                d=[];
                skel= bwmorph(BW,'skel',Inf);
                Dmask = false(size(skel));
                for k = 1:numel(c_e)
                    D = bwdistgeodesic(skel,c_e(k),r_e(k))+BW;
                    
                end
                for j=1:size(r_e,1)
                    if ~isempty(x_b)
                        d_e_b(j,1)=sqrt((r_e(j)-x_b)^2+(c_e(j)-y_b)^2);
                        j=j+1;
                    else check=1;
                    end
                end
                if any(d_e_b<sqrt(32))
                    BW(r_e(d_e_b(:,1)<sqrt(32),1), c_e(d_e_b(:,1)<sqrt(32),1))=0;
                    props=regionprops(BW, 'Area');
                    BW=bwareaopen(BW, max(cell2mat(struct2cell(props)))-1);
                    endpoints=bwmorph(BW, 'endpoint');
                    endpoints=find(endpoints~=0);
                    [r_e,c_e]=ind2sub(N, endpoints);
                    Br=bwmorph(BW, 'branchpoints');
                    idx=find(Br~=0);
                    [r_b, c_b]=ind2sub(N, idx);
                else
                    r=1;
                    for n=x_b-4:x_b+4
                        for m=y_b-4:y_b+4
                            if n==x_b-4 & ~isnan(D(n,m)) | n==x_b+4 & ~isnan(D(n,m))| m==y_b+4 & ~isnan(D(n,m)) |  m==y_b-4 & ~isnan(D(n,m))
                                d(r,1)=abs(D(n,m)-D(x_b, y_b));
                                d(r,2)=abs(D(n,m)-D(r_e,c_e));
                                r=r+1;
                            end
                        end
                    end
                    d_min=min(d(:,2));
                    [x_delete, y_delete]=ind2sub(size(D), find(D==D(r_e,c_e)+ d_min | D==D(r_e,c_e)- d_min));
                    BW(x_delete,y_delete)=0;
                    props=regionprops(BW, 'Area');
                    BW=bwareaopen(BW, max(cell2mat(struct2cell(props)))-1);
                    endpoints=bwmorph(BW, 'endpoint');
                    endpoints=find(endpoints~=0);
                    [r_e,c_e]=ind2sub(N, endpoints);
                    Br=bwmorph(BW, 'branchpoints');
                    idx=find(Br~=0);
                    [r_b, c_b]=ind2sub(N, idx);
                                       
                end
                if size(r_e, 1)>2 | ~isempty(r_b)
                    Br=bwmorph(BW, 'branchpoints');
                    idx=find(Br~=0);
                    [r_b, c_b]=ind2sub(N, idx);
                    if any(r_b==x_b) && any(c_b==y_b)
                        for j=1:size(r_e,1)
                            d_e_b(j,1)=sqrt((r_e(j)-x_b)^2+(c_e(j)-y_b)^2);
                            j=j+1;
                            
                        end
                        
                        if any(d_e_b(:,1)<=sqrt(32))
                            check=0;
                            d_e_b=[];
                            mat=[];
                            matrix=[];
                            
                        else
                            
                            check=0;
                            o=o;
                            
                            
                        end
                    else
                       
                        check=0;
                        o=1;
                        
                    end
                else
                    check=1; 
                end
            end
        end
    end
    
    angle=zeros(size(r_e,1),1);
    alfa=zeros(size(r_e,1),1);
    
    for t=1:size(r_e,1)
        if r_e(t,1)>=r_c && c_e(t,1)>=c_c
            angle(t,1)= (-atan((r_e(t,1)-r_c)/(c_e(t,1)-c_c))*180/pi)+360;
        elseif r_e(t,1)>=r_c && c_e(t,1)<c_c
            angle(t,1)=(-atan((r_e(t,1)-r_c)/(c_e(t,1)-c_c))*180/pi)+180;
        elseif r_e(t,1)<r_c && c_e(t,1)<c_c
            angle(t,1)= ((-atan((r_e(t,1)-r_c)/(c_e(t,1)-c_c))*180/pi))+180;
        elseif r_e(t,1)<r_c && c_e(t,1)>=c_c
            angle(t,1)=(-atan((r_e(t,1)-r_c)/(c_e(t,1)-c_c))*180/pi);
            
        end
    end
    if size(r_e,1)==2
        if r_e(1)>r_c && r_e(2)==r_c && angle(2)==360;
            angle(2)=0;
        elseif r_e(1)<r_c && r_e(2)==r_c && angle(2)==360;
            angle(2)=360;
        end
        [angle, idx]=sort(angle(:,1), 'descend');
        r_e=r_e(idx,:);
        c_e=c_e(idx, :);
        
        re_1=r_e(1);
        re_2=r_e(2);
        ce_1=c_e(1);
        ce_2=c_e(2);
        
    else
        r_e=[];
        c_e=[];
        
        re_1=[];
        re_2=[];
        ce_1=[];
        ce_2=[];
    end
else
    
    r_e=[];
    c_e=[];
    re_1=[];
    re_2=[];
    ce_1=[];
    ce_2=[];
    angle=[];
    threshout=[];
end
%