%% CEX 1

close all
run vl_setup.m;
ima=imread("a.jpg");
imb=imread("b.jpg");

figure()
imagesc(ima);
figure()
imagesc(imb);

%SIFT
[fA dA ] = vl_sift( single( rgb2gray( ima )) ); %fA is 4x947, so 947 features found in ima
[fB dB ] = vl_sift( single( rgb2gray( imb )) ); %fB is 4x865, so 865 features found in ima
matches = vl_ubcmatch( dA , dB ); %2x204 so 204 matches
xA = fA(1:2 , matches (1 ,:));
xB = fB(1:2 , matches (2 ,:));

%Now we can find homography H, transforming xA<->xB using DLT
%Ransac to find good correspondences.
threshold=5; %threshold for inliers = 5 pixels
min=4; %minimum nmbr of correspondences needed
%random index
perm = randperm (size(matches,2));

%Normalize it
xAmean=mean(xA(1:2,:),2); %(x0,y0)
sxA=std(xA(1:2,:),0,2); %(sx,sy)
%Normalization matrix
NA=[1/sxA(1), 0, -xAmean(1)/sxA(1);0, 1/sxA(2), -xAmean(2)/sxA(2);0, 0, 1];
xAn=NA*[xA;ones(1,length(xA))];
%Normalize it
xBmean=mean(xB(1:2,:),2); %(x0,y0)
sxB=std(xB(1:2,:),0,2); %(sx,sy)
%Normalization matrix
NB=[1/sxB(1), 0, -xBmean(1)/sxB(1);0, 1/sxB(2), -xBmean(2)/sxB(2);0, 0, 1];
xBn=NB*[xB;ones(1,length(xB))];

maxnmbr=0;
Hbest=[];
iters=31;
for k=1:1:iters
    M=[];
    x1=[xAn(:,perm(4*k-3:4*k))];
    x2=[xBn(:,perm(4*k-3:4*k))];

    m=[];
    for i=1:1:4
        l=[[x1(:,i)',zeros(1,6)];[zeros(1,3),x1(:,i)',zeros(1,3)];[zeros(1,6),x1(:,i)']];
        r=zeros(3,4);
        r(1:3,i)=-x2(:,i);
        m=[l,r];
        M=[M;m];
    end

    [U,S,V]=svd(M);
    v=V(1:9,end);
    %lams=V(end,10:end);
    norm(M*V(:,end))

    Hn=reshape(v,[3,3])';
    H=inv(NB)*Hn*NA;
    xprojB=pflat(H*[xA;ones(1,length(xA))]);

    
    good_points = (sqrt(sum((xB(1:2,:) - xprojB(1:2,:)).^2))< threshold);
    nmbrgoodpoints=sum(good_points);

    if nmbrgoodpoints>maxnmbr
        maxnmbr=nmbrgoodpoints;
        Hbest=H;
        maxgoodpoints=good_points;
    end
end
xprojB=pflat(Hbest*[xA;ones(1,length(xA))]);

figure()
imagesc(imb)
hold on
plot(xprojB(1,maxgoodpoints),xprojB(2,maxgoodpoints),'co')
hold on
plot(xB(1,maxgoodpoints),xprojB(2,maxgoodpoints),'r*')
hold off

Htform = projective2d(Hbest');

Rout = imref2d( size( ima ) ,[ -200 800] ,[ -400 600]);
% Sets the size and output bounds of the new image .
[Atransf] = imwarp(ima , Htform ,  "OutputView" , Rout);
% Transforms the image
Idtform = projective2d ( eye (3));
[ Btransf ] = imwarp (imb , Idtform ,  "OutputView", Rout );
% Creates a larger version of the second image
AB = Btransf;
AB ( Btransf < Atransf ) = Atransf ( Btransf < Atransf );
% Writes both images in the new image . %( A somewhat hacky solution is needed
% since pixels outside the valid image area are not always zero ...)
figure()
imagesc ( Rout.XWorldLimits , Rout . YWorldLimits , AB );

% Plots the new image with the correct axes

%%


