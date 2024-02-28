%Loic Daumail, 01-31-2024
%% Random Lines
iw = 1280/2; %devide by 2 since the screen will be split in half for each mirror of the stereoscope, in pixels
ih = 1024;
id = round(sqrt(iw^2+ih^2)); %we compute the diagonal distance of the screen: as we will apply a rotation, the matrix that we will initially fill with line segments will be id x id in dimensions
% pixelpitch = .252; %same for horizontal and vertical in mm
viewdist = 46; %in cm
iwd = 47;%asind( iw*pixelpitch/(viewdist*10)); %screen dimensions in degrees
ihd = 35.25;%asind( ih*pixelpitch/(viewdist*10));
idd = round(sqrt(iwd^2+ihd^2));%asind( id*pixelpitch/(viewdist*10));

%we want this density of line segments:
linedensity = 11; %11 line segments per square degree
nlines = round(linedensity*(iwd*ihd)); %total number of lines needed for this density
linelen = 16;
cont =1;
ori = 45;%[45, 135;135,45]; %orientation of background and figure

origBitmap = ones(id+10+linelen, id+10+linelen);%preallocate image of random line segments, add 10 to make lines start at indices diffrent from 1, to avoid having a "start effect"
%add linelen as well to remove black part on right side of the image

lindex = randi(id+10+linelen,nlines,1);
colindex = randi(id+10+linelen,nlines,1);
for i =1:nlines
    origBitmap(lindex(i),colindex(i):colindex(i)+linelen-1) = 0;
end

img = origBitmap(10:end-linelen-1,10:id+10-1);

% figure();
%imshow(img)
%saveDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\Saliency_study\experiment_design\line_stims\';
%saveas(gcf,strcat(saveDir, sprintf('figGround_random_line_segm_ori%d_cont%d_len%d_xlocation%d_ylocation%d.png', 100*orientations(d), 100*cont(c), 10000*scalingF(s), locsComb(1,a),locsComb(2,a))));
%saveas(gcf,strcat(saveDir, 'figGround_random_line_segm.png'));


%% rotate image at 45 dva and 135 dva

%store img in a uint8 format for functions to work
A = uint8(ceil((img.*cont).*255)); %put intensity scale from 0 to 255 and convert into uint8 %trim image to size of screenWidth

%Mid point of the image
midx=ceil((size(A,1)+1)/2);
midy=ceil((size(A,2)+1)/2);

x2=zeros([size(A,1) size(A,2)]);
y2=zeros([size(A,1) size(A,2)]);

%ori = [45, 135]; %background 45 degrees, figure 135 degrees
%ori = [135, 45]; %background 135 degrees, figure 45 degrees
B = uint8(zeros([size(A),length(ori)])); % preallocate for rotated images
for s = 1:length(ori)
    for i=1:size(A,1)
        for j=1:size(A,2)
            %Cartesian to Polar co-ordinates
            [theta1,rho1]=cart2pol(i-midx,j-midy);
            phi=theta1+ori(s)*pi/180;
            %Polar to Cartesian co-ordinates
            [l,m]=pol2cart(phi,rho1);
            x2(i,j)=ceil(l)+midx;
            y2(i,j)=ceil(m)+midy;
        end
    end
    %The result may produce value lesser than 1 or greater than the image size.
    x2=max(x2,1);
    x2=min(x2,size(A,1));
    
    y2=max(y2,1);
    y2=min(y2,size(A,2));
    
    for i=1:size(A,1)
        for j=1:size(A,2)
            B(i,j,s)=A(x2(i,j),y2(i,j));
        end
    end
end
img2 =  B(round((id-ih)/2):round((id-ih)/2)+ih-1, round((id-iw)/2):round((id-iw)/2)+iw-1,:);

background = img2(:,:,1);
% fig = img2(:,:,2);

figure();
imshow(background)

%% Grid lines

ex.lineColor = [0 0 0];
ex.lineHdeg = 1/2;
ex.lineThdeg = 0.1;
ex.lineHdegSq = 5;
ex.lineH = ex.lineHdeg*ex.ppd;
ex.lineW = ex.lineHdeg*ex.ppd;
ex.lineTh = ex.lineThdeg*ex.ppd;
ex.lineHsq = ex.lineHdegSq*ex.ppd;
%% RIGHT SQUARE
% Grid Line coordinates
% 45 deg line coordinates
xRforLl =  3 * xc / 2 - ex.lineW;
xRforLr =  3 * xc / 2 + ex.lineW;
yRforLBot = yc + ex.lineH;
yRforLTop = yc - ex.lineH;

% 135 deg line coordinates
xRbackRl = 3 * xc / 2 - ex.lineW;
xRbackRr = 3 * xc / 2 + ex.lineW;
yRbackRTop = yc - ex.lineH;
yRbackRBot = yc + ex.lineH;
%            Grid lines
            for c = 1:3
                xshift = 2*ex.lineH*(c-1);
                for r=1:3
                    yshift = 2*ex.lineW*(r-1);
                    %lower right cross/grid
                    Screen('DrawLines', w, [xRforLl+xshift, xRforLr+xshift; yRforLBot+yshift, yRforLTop+yshift], ex.lineTh, ex.lineColor);
                    Screen('DrawLines', w, [xRbackRl+xshift, xRbackRr+xshift; yRbackRTop+yshift, yRbackRBot+yshift], ex.lineTh, ex.lineColor);
                    %Upper right grid
                    Screen('DrawLines', w, [xRforLl+xshift, xRforLr+xshift; yRforLBot-yshift, yRforLTop-yshift], ex.lineTh, ex.lineColor);
                    Screen('DrawLines', w, [xRbackRl+xshift, xRbackRr+xshift; yRbackRTop-yshift, yRbackRBot-yshift], ex.lineTh, ex.lineColor);
                    
                end
            end
            for c = 1:3
                xshift = 2*ex.lineH*(c-1);
                for r=1:3
                    yshift = 2*ex.lineW*(r-1);
                    %lower left cross/grid
                    Screen('DrawLines', w, [xRforLl-xshift, xRforLr-xshift; yRforLBot+yshift, yRforLTop+yshift], ex.lineTh, ex.lineColor);
                    Screen('DrawLines', w, [xRbackRl-xshift, xRbackRr-xshift; yRbackRTop+yshift, yRbackRBot+yshift], ex.lineTh, ex.lineColor);
                    
                    %Upper left grid
                    Screen('DrawLines', w, [xRforLl-xshift, xRforLr-xshift; yRforLBot-yshift, yRforLTop-yshift], ex.lineTh, ex.lineColor);
                    Screen('DrawLines', w, [xRbackRl-xshift, xRbackRr-xshift; yRbackRTop-yshift, yRbackRBot-yshift], ex.lineTh, ex.lineColor);
                    
                end
            end
            %     %% square
            % Define the square parameters (position and size)
            squareRectRight = [3 * xc / 2-ex.lineHsq/2, yc - ex.lineHsq/2, 3 * xc / 2+ex.lineHsq/2, yc + ex.lineHsq/2];
            % Draw the square
            Screen('FrameRect', w, [0 0 0], squareRectRight, ex.lineTh);