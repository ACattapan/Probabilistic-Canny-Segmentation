function [Shape, threshOut]=Segmentation_function_else(img,img2, rect, w) 

[r_e, c_e, r_c, c_c, BW, re_1, re_2, ce_1, ce_2, angle, threshOut]= Endpoints(img, img2);

gray=rgb2gray(img);

if ~isempty(r_e) & ~isempty(c_e) &any(r_e>2) & any(c_e>2)
    delta_t=6;
    Dmask = false(size(BW));
    for k = 1:numel(c_e)
        D = bwdistgeodesic(BW,c_e(k),r_e(k));
       
    end
    i=1;
    v=0;
    connected_points=zeros(25,5);
    while v==0
        d=1
        
        for y=r_e(i)-4:r_e(i)+4
            for g=c_e(i)-4:c_e(i)+4
                if y==r_e(i)-4 & ~isnan(D(y,g)) | y==r_e(i)+4 & ~isnan(D(y,g)) | g==c_e(i)+4 & ~isnan(D(y,g))  |  g==c_e(i)-4 & ~isnan(D(y,g))
                    connected_points(d,1)=y;
                    connected_points(d,2)=g;
                    connected_points(d,3)=D(y,g);
                    connected_points(d,4)=D(r_e(i),c_e(i))-D(y,g);
                    connected_points(d,5)=sqrt((y-r_e(i))^2+(g-c_e(i))^2);
                    d=d+1;
                end
            end
        end
        
        con= connected_points(connected_points(:,4)>10 & connected_points(:,3)~=0, :);
        con=con(con(:,5)==min(con(:,5)),:);
        
        if ~isempty(con)
            X0=r_e(i);
            Y0=c_e(i);
            X1=con(con(:,4)==min(con(:,4)),1);
            Y1=con(con(:,4)==min(con(:,4)),2);
            BW=Drawline2D(BW,X0,Y0,X1,Y1);
            endpoints=bwmorph(BW, 'endpoint');
            endpoints=find(endpoints~=0);
            [r_e,c_e]=ind2sub(size(img), endpoints);
            Dmask = false(size(BW));
            for k = 1:numel(c_e)
                D = bwdistgeodesic(BW,c_e(k),r_e(k));
                
            end
            
        end
        if ~isempty(r_e) & ~isempty(con)
            connected_points=[];
            con=[];
            i=size(r_e,1);
            v=0;
        else v=1;
        end
    end
    
    if   size(r_e,1)==2
        
        % imshow(gray)
        N=size(gray);
        n_row=N(1);
        n_col=N(2);
        se = strel('disk',15,8);
        Ge = imerode(gray,se);
        Gobr = imreconstruct(Ge,gray);
        % imshow(Gobr)
        se1=strel('disk',10,8);
        Gobrd = imdilate(Gobr,se1);
        Gobrcbr = imreconstruct(imcomplement(Gobrd),imcomplement(Gobr));
        Gobrcbr = imcomplement(Gobrcbr);
        %COMPUTE A SERIES OF EDGE DETENCTION CHANGING PARAMETERS: ST.
        %DEVIATION, THRESHOLD INTERVAL
        %parameters for st. deviation: Lower and upper limits
        %two thresholds UL and LL. edge disregards all edges with edge strength below the lower threshold,
        %and preserves all edges with edge strength above the higher threshold. The
        %positive thing of this algorithm is that is able to detect not only the
        %stronger edges, but also the weaker connected to the stronger ones. The
        %algorithm defines indeed two threshold, one for the stronger and the other
        %for the weaker.
        %FROM MATLAB: The Canny method differs from the other edge-detection methods in that
        %it uses two different thresholds (to detect strong and weak edges), and includes
        %the weak edges in the output only if they are connected to strong edges.
        %This method is therefore less likely than the others to be affected by noise, and
        %more likely to detect true weak edges.
        
        sd1= [0.2]%,sqrt(2)];
        sd2= [5]%, sqrt(2)];
        Deltasd=0.1;
        lowm= 0.1;
        lowM= [0.2]
        highm= [0.3]
        highM= [0.3]
        Deltalim=0.1;
        Dk= [8];
        Dh= [8];
        A=[];
        C_coord=[];
        
        
        i=1;
        for lowindx = lowm:Deltalim:lowM
            if lowindx==0.0
                lowindx=[];
            end
            for highindx = highm:Deltalim:highM
                if highindx==0.0
                    highindx=[];
                end
                
                if or(isempty(highindx),isempty(lowindx))
                    
                    for sd = sd:Deltasd:sd2
                        E_s(:,:,i)= edge(gray,'Canny',[lowindx highindx],sd);
                        
                        i=i+1;
                        
                    end
                    
                elseif lowindx<highindx
                    
                    for sd = sd1:Deltasd:sd2
                        E_s(:,:,i)= edge(Gobrcbr,'Canny',[lowindx highindx],sd);
                        
                        i=i+1;
                    end
                    
                end
            end
        end
        
        %%---------------------------------------
        %Last image obtained by the application of Canny edge without
        %imposing the parameters
        E_s(:,:,i)=edge(Gobrcbr, 'Canny');
        %-------------------------------------------------------
        %%External function MeanEdges. From the edges (value in each pixel==1) detected by the function
        %%edge(Canny) for each combination of the parameters on which the function is based, the MeanEdges,
        %%calculates the mean value in each pixel. The output is an image where in
        %%each pixel, is saved a value from 0 to 1 which is a function of the
        %%frequency with which the edge in that pixel is identified.
        
        E_k(:,:)= mean(E_s(:,:,:),[3]);
        %E_k_1=imdilate(E_k, strel('disk', 1,8));
        % cd 'C:\Users\agu005\Desktop\Master Thesis Gurini\Image_Processing'
        %figure
        %imshow(E_k_1)
       % print('Canny_probabilistico', '-dtiff', '-r1000');
        %-------------------------------------------------------------------
        %Transformation of the cartesian coordinates into polar ones
        for t=1:size(r_e,1)
            dist(t,1)=round(sqrt((r_e(t,1)-r_c)^2+ (c_e(t,1)-c_c)^2));
        end
        %polar coordinates matrix
        p_coordinate=[angle, dist];
        %-----------------------------
        %Definition of the angle width between the two endpoints and linear
        %distance
        if p_coordinate(1,1)- p_coordinate(2,1)<180
            delta=p_coordinate(1,1)- p_coordinate(2,1);
        else
            delta=360 - p_coordinate(1,1)+ p_coordinate(2,1);
        end
        d=sqrt((r_e(1)-r_e(2))^2+(c_e(1)-c_e(2))^2);
        %%------------------initialization of the variables
        %delta theshold decides when the loop ends
        delta_t=6;
        k=r_e(1);
        h=c_e(1);
        %q and p are two indices needed for saving data in vectors
        q=1;
        p=1;
        check=0;
        BW1=BW;
        
        while check==0
            %definition of the matrix of the possible nonzero values around the current
            %endpoint analyzed, with which the longest connected component can
            %merge
            if (k-Dk)<=0 | (k+Dk)>=size(img, 1) | (h-Dh)<=0 | (h+Dh)>=size(img,2)
                check=1;
            else
                E_w= E_k(k-Dk:k+Dk, h-Dh:h+Dh);
                %first condition for selecting which endpoint can be considered. If
                %there is no point with which the endpoint can be connected, it is
                %necessary to pass to the other endpoint.
                
                if any(E_w(:,:),'all')
                    r=1;
                    for n=k-Dk:k+Dk
                        for m=h-Dh:h+Dh
                            if E_k(n,m)~=0
                                %% definition of a new matrix A in which are registered information
                                %%about row and column of the pixel inside the window
                                %%+- Dk, +-Dh, its E_k value, the polar coordinates
                                %%with respect to the centroid, angular distance to the
                                %%other endpoint,slope of the function E_k from the
                                %%endpoint (k,h) and the neighboor (n,m)
                                A(r,1)=n; %row
                                A(r,2)=m; %column
                                A(r,3)=E_k(n,m); %value of the edge strength
                                %angular distance from the centroid
                                if n>=r_c && m>=c_c
                                    A(r,4)= (-atan((n-r_c)/(m-c_c))*180/pi)+360;
                                elseif n>=r_c && m<c_c
                                    A(r,4)=(-atan((n-r_c)/(m-c_c))*180/pi)+180;
                                elseif n<r_c && m<c_c
                                    A(r,4)= (-atan((n-r_c)/(m-c_c))*180/pi)+180;
                                elseif n<r_c && m>=c_c
                                    A(r,4)= (-atan((n-r_c)/(m-c_c))*180/pi);
                                    
                                end
                                %distance from the actual endpoint
                                A(r,5)=sqrt((n-r_e(q))^2+(m-c_e(q))^2);
                                %delta angle between the pixel (n,m) and the other
                                %endpoint. It is necessary to consider two different
                                %cases dependent on which endpoint it is considered.
                                %Indeed, since the endpoints have been ordered
                                %according to the angular position with respect to the
                                %centroid (alfa of the first endpoint is lower than the
                                %angle of the second), if the first e.p is considered
                                %the angular distance is A(r,4)+360-alfa(2), otherwise
                                %it is  360-A(r,4)+alfa(1)
                                
                                if q~=size(r_e,1) && p_coordinate(q+1,1)> A(r,4) && p_coordinate(q+1,1)-A(r,4) <180
                                    A(r,6)= p_coordinate(q+1,1)-A(r,4);
                                elseif q~=size(r_e,1) && p_coordinate(q+1,1)> A(r,4) && p_coordinate(q+1,1)-A(r,4)>180
                                    A(r,6) = 360+ A(r,4) - p_coordinate(q+1,1);
                                elseif q~=size(r_e,1)&& p_coordinate(q+1,1)< A(r,4)&& A(r,4)- p_coordinate(q+1,1) <180
                                    A(r,6)= A(r,4) - p_coordinate(q+1,1);
                                elseif q~=size(r_e,1) && p_coordinate(q+1,1)< A(r,4)&& A(r,4)- p_coordinate(q+1,1)>180
                                    A(r,6) = 360-A(r,4)+ p_coordinate(q+1,1);
                                elseif q==size(r_e,1) && p_coordinate(q-1,1)> A(r,4) && p_coordinate(q-1,1)-A(r,4) <180
                                    A(r,6)= p_coordinate(q-1,1)-A(r,4);
                                elseif q==size(r_e,1) && p_coordinate(q-1,1)> A(r,4) && p_coordinate(q-1,1)-A(r,4)>180
                                    A(r,6) = 360+ A(r,4) - p_coordinate(q-1,1);
                                elseif q==size(r_e,1)&& p_coordinate(q-1,1)< A(r,4)&& A(r,4)- p_coordinate(q-1,1) <180
                                    A(r,6)= A(r,4) - p_coordinate(q-1,1);
                                elseif q==size(r_e,1) && p_coordinate(q-1,1)< A(r,4)&& A(r,4)- p_coordinate(q-1,1)>180
                                    A(r,6) = 360-A(r,4)+ p_coordinate(q-1,1);
                                end
                                
                                %Difference between the values of edge strength assumed
                                %in the current endpoint and the pixel n, m
                                A(r,7)=E_k(n,m)-E_k(r_e(q),c_e(q));
                                %Definition of the gradient of the function E_k
                                A(r,8)=A(r,7)/A(r,5);
                                
                                %Distance of the pixel n,m from the centroid
                                A(r,9)=sqrt((n-r_c)^2+(m-c_c)^2);
                                A(r,10)=abs(A(r,9)-p_coordinate(q,2));
                                r=r+1;
                                
                            end
                        end
                    end
                    %CREATION of new matrices (B, M) which contain the pixels that
                    %soddisfy the conditions required of delta lower than the current
                    %value, distance from the actual endpoint ~=0, distance from the
                    %centroid similar to the distance value assumed in the current
                    %endpoint
                    
                    
                    B=A(A(:,6)<delta & A(:,5)~=0,:);
                    if ~isempty(B)  %&(A(:,6)<delta) | A(:,5)~=0
                        %MAtrix M represents the pixel that becomes the new endpoint. Condition of maximum gradient
                        %minimum distance from the endpoint, and minimum angular
                        %distance from the other endpoint
                        
                        M=B(B(:,8)==max(B(:,8)),:);
                        M=M(M(:,5)==min(M(:,5)),:);
                        M=M(M(:,6)==min(M(:,6)),:);
                        M=M(M(:,10)==min(M(:,10)), :);
                        
                        %interpolation of the previous endpoint with the new one though
                        %a line defined with the function Drawline2D
                        X0=r_e(q);
                        Y0=c_e(q);
                        X1=M(:,1);
                        Y1=M(:,2);
                        BW1=Drawline2D(BW1,X0,Y0,X1,Y1);
                        r_e(q)=X1;
                        c_e(q)=Y1;
                        p_coordinate(q,1)=M(:,4);
                        %Updating of the matrices that store the polar coordinates of
                        %the endpoints.
                        p_coordinate(q,2)=sqrt((r_e(q)-r_c)^2+(c_e(q)-c_c)^2);
                        
                        delta=M(:,6);
                    elseif any(A(:,6)==delta & A(:,5)~=0)
                        B=A(A(:,6)==delta & A(:,5)~=0,:); %& A(:,9)<p_coordinate(q,2) & A(:,9)>p_coordinate(q,2)-8
                        
                        %MAtrix M represents the pixel that becomes the new endpoint. Condition of maximum gradient
                        %minimum distance from the endpoint, and minimum angular
                        %distance from the other endpoint
                        
                        M=B(B(:,8)==max(B(:,8)),:);
                        M=M(M(:,5)==min(M(:,5)),:);
                        M=M(M(:,6)==min(M(:,6)),:);
                        M=M(M(:,10)==min(M(:,10)), :);
                        
                        %interpolation of the previous endpoint with the new one though
                        %a line defined with the function Drawline2D
                        X0=r_e(q);
                        Y0=c_e(q);
                        X1=M(:,1);
                        Y1=M(:,2);
                        BW1=Drawline2D(BW1,X0,Y0,X1,Y1);
                        r_e(q)=X1;
                        c_e(q)=Y1;
                        p_coordinate(q,1)=M(:,4);
                        p_coordinate(q,2)=sqrt((r_e(q)-r_c)^2+(c_e(q)-c_c)^2);
                        
                        delta=M(:,6);
                    elseif  any(A(:,5)~=0)
                        B=A((A(:,6)-delta==min(A(A(:,1)~=r_e(q) & A(:,2)~=c_e(q),6)-delta)),:); %& A(:,9)<p_coordinate(q,2) & A(:,9)>p_coordinate(q,2)-8
                        
                        %MAtrix M represents the pixel that becomes the new endpoint. Condition of maximum gradient
                        %minimum distance from the endpoint, and minimum angular
                        %distance from the other endpoint
                        
                        M=B(B(:,8)==max(B(:,8)),:);
                        M=M(M(:,5)==min(M(:,5)),:);
                        M=M(M(:,6)==min(M(:,6)),:);
                        M=M(M(:,10)==min(M(:,10)), :);
                        
                        %interpolation of the previous endpoint with the new one though
                        %a line defined with the function Drawline2D
                        X0=r_e(q);
                        Y0=c_e(q);
                        X1=M(:,1);
                        Y1=M(:,2);
                        BW1=Drawline2D(BW1,X0,Y0,X1,Y1);
                        r_e(q)=X1;
                        c_e(q)=Y1;
                        p_coordinate(q,1)=M(:,4);
                        p_coordinate(q,2)=sqrt((r_e(q)-r_c)^2+(c_e(q)-c_c)^2);
                        delta=M(:,6);
                    end
                    
                end
                
                %Mtrix C_coord which stores the value of row and column of the endpoint
                %at each step. It has been created in order to check whether the both
                %endpoints can't be connected with other neighboor pixels inside the
                %window. When this occur, the loop ends and the two endpoint can be
                %connected with a line.
                C_coord(p,1)=r_e(q);
                C_coord(p,2)=c_e(q);
                if p<4 && delta>delta_t && d>sqrt(Dk^2+Dh^2)
                    if q==size(r_e,1)
                        q=q-1;
                        k=r_e(q);
                        h=c_e(q);
                    else
                        q=q+1;
                        k=r_e(q);
                        h=c_e(q);
                    end
                    delta=delta;
                    d=sqrt((r_e(1)-r_e(2))^2+(c_e(1)-c_e(2))^2);
                    A=[];
                    p=p+1;
                    check=0;
                elseif p<4 && delta<delta_t && d<sqrt(Dk^2+Dh^2)
                    check=1; %break
                elseif p>=4 & size(find(C_coord(:,1)==C_coord(p-1,1) & C_coord(:,2)==C_coord(p-1,2)),1)~=1 & size(find(C_coord(:,1)==C_coord(p,1) & C_coord(:,2)==C_coord(p,2)),1)~=1 |  delta<=delta_t | d<=sqrt(Dk^2+Dh^2) %C_coord(p-1,1)==C_coord(C_coord(p-1,2)==C_coord(1:p-3,2),1) & C_coord(p,1)==C_coord(C_coord(p, 2)==C_coord(1:p-2,2),1) |  delta<=delta_t | d<=sqrt(Dk^2+Dh^2) %any(C_coord(p-1,:)==C_coord(1:p-3,:))& any(C_coord(p,1)==C_coord(1:p-2,1))
                    check=1; %break
                else                    
                    if q==size(r_e,1)
                        q=q-1;
                        k=r_e(q);
                        h=c_e(q);
                    else
                        q=q+1;
                        k=r_e(q);
                        h=c_e(q);
                    end
                    delta=delta;
                    d=sqrt((r_e(1)-r_e(2))^2+(c_e(1)-c_e(2))^2);
                    A=[];
                    p=p+1;
                    check=0;
                    
                end
                
            end
        end
        X0=r_e(1);
        Y0=c_e(1);
        X1=r_e(2);
        Y1=c_e(2);
        
        BW1=Drawline2D(BW1,X0,Y0,X1,Y1);
        Shape=BW1;
        %         Shape=imdilate(Shape, strel('disk',4,8));
        Shape_dilate=imdilate(Shape, strel('disk', 1,8));
        Shape_rgb=label2rgb(Shape_dilate);
        Shape_gray=rgb2gray(Shape_rgb);
         
        figure
        a=imshow(img);
        hold on
        a=imshow(Shape_rgb);
        hold off
        alpha(a, 255-Shape_gray);
        print('Result', '-dtiff','-r1000');
        
    else
        Shape=BW;
        Shape_dilate=imdilate(Shape, strel('disk', 1,8));
    end
else
    Shape=BW;
    Shape_dilate=imdilate(Shape, strel('disk', 1,8));
end

