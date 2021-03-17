close all; clc; clear
% Two methods to view the video
% Method 1: use MATLAB'd built-in video player
% load('cam1_1.mat')
% implay(vidFrames1_1)

% Method 2: create a loop and display each frame
% numFrames = size(vidFrames1_1, 4);
% for j = 1:numFrames
%     X = vidFrames1_1(:,:,:,j);
%     imshow(X); drawnow
% end

% Test 1
% Cam1-1
% -- Load Movie
% load ('cam1_1.mat');
% [height width rgb nf1_1] = size(vidFrames1_1);
% %  height, width, rgb, num_frames
% -- Watch Movie
% X11 = [];
% Y11 = [];
% for j = 1:nf1_1
%     X = vidFrames1_1(:,:,:,j);
% %   consider displacement of mass in z direction 
% %     imshow(X); drawnow
% % clean out the noises by just blackening 
% % the background, surrounding images
% % Take the cam1_1 as an example, X here is 480*640*3,
% % and the spring-mass is oscilating around pixels around
% % (:,240:480) and (120:460,:)
%     X(:,1:240) = 0;
%     X(:,480:end) = 0;
%     X(1:40,:) = 0;
%     X(460:end, :) = 0;
%     [Max Ind] = max(X(:));
%     [y11 x11] = ind2sub(size(X), Ind);
%     X11 = [X11 x11];
%     Y11 = [Y11 y11];
%     plot(X11,'b','Linewidth',[1]);
%     hold on
%     plot(Y11,'r','Linewidth',[1]);
%     legend('Position on x', 'Position on y','Location','Best')
%     title('Cam1-1')
%     xlabel('Number of Frames')
%     ylabel('Positions')
%     set(gca, 'Fontsize', [12])
% end

% Cam2-1
% -- Load Movie
% load ('cam2_1.mat');
% % height, width, rgb, num_frames
% [height width rgb nf2_1] = size(vidFrames2_1);
% % -- Watch Movie
% X21 = [];
% Y21 = [];
% for j = 1:nf2_1
%     X = vidFrames2_1(:,:,:,j);
% %     imshow(X); drawnow
%     X(:,1:240) = 0;
%     X(:,400:end) = 0;
% %     imshow(X)
%     [Max Ind] = max(X(:));
%     [y21 x21] = ind2sub(size(X), Ind);
%     X21 = [X21 x21];
%     Y21 = [Y21 y21];
%     plot(X21,'b','Linewidth',[1]);
%     hold on
%     plot(Y21,'r','Linewidth',[1]);
%     legend('Position on x', 'Position on y','Location','Best')
%     title('Cam2-1')
%     xlabel('Number of Frames')
%     ylabel('Positions')
%     set(gca, 'Fontsize', [12])
% end

% Cam3-1
% -- Load Movie
% load ('cam3_1.mat');
% % height, width, rgb, num_frames
% [height width rgb nf3_1] = size(vidFrames3_1);
% % -- Watch Movie
% X31 = [];
% Y31 = [];
% for j = 1:nf3_1
%     X = vidFrames3_1(:,:,:,j);
% %     imshow(X); drawnow
% %     keep around(200:360,:) and (:,1:510)
%     X(1:230,:) = 0;
%     X(340:end,:) = 0;
%     X(:,1:60) = 0;
%     X(:,500:end) = 0;
% %     imshow(X)
%     [Max Ind] = max(X(:));
%     [y31 x31] = ind2sub(size(X), Ind);
%     X31 = [X31 x31];
%     Y31 = [Y31 y31];
%     plot(X31,'b','Linewidth',[1]);
%     hold on
%     plot(Y31,'r','Linewidth',[1]);
%     legend('Position on x', 'Position on y','Location','Best')
%     title('Cam3-1')
%     xlabel('Number of Frames')
%     ylabel('Positions')
%     set(gca, 'Fontsize', [12])
% end
% 
% % To make the num_frames to be the same that equals to the shortest,
% % we need to find the minimum length by using:
% Xmin = min([length(X11) length(X21) length(X31)]);
% X11 = X11(1:Xmin);
% X21 = X21(1:Xmin);
% X31 = X31(1:Xmin);
% Y11 = Y11(1:Xmin);
% Y21 = Y21(1:Xmin);
% Y31 = Y31(1:Xmin);
% matrix = [X11;Y11;X21;Y21;X31;Y31];
% [m,n] = size(matrix);
% % Find the mean
% miu = mean(matrix,2);
% % miu = sum(matrix)/(h*w);
% new_matrix = zeros(m,n);
% % Subtract the mean from all entries in matrix we have
% for i = 1:n
%     for j = 1:m
%         new_matrix(j,i) = matrix(j,i) - miu(j,1);
%     end
% end
% % Find covariance by the formula
% Cx = (1/(n-1)) * new_matrix * new_matrix';
% [V,D] = eig(Cx);
% lambda = diag(D);
% Y = V'*new_matrix;
% plot(Y(1,:),'Linewidth',[1])
% title('Test 1 - Ideal Case')
% xlabel('Number of Frames')
% ylabel('Positions')
% set(gca, 'Fontsize', [12])

% % Test2 Noisy case
% % Cam1-2
% load ('cam1_2.mat');
% [height width rgb nf2] = size(vidFrames1_2);
% X12 = [];
% Y12 = [];
% % -- Watch Movie
% for j = 1:nf2
%     X = vidFrames1_2(:,:,:,j);
% %     imshow(X); drawnow
%     X(:,1:270) = 0;
%     X(:,480:end) = 0;
% % %      imshow(X)
%     [Max Ind] = max(X(:));
%     [y12 x12] = ind2sub(size(X), Ind);
%     X12 = [X12 x12];
%     Y12 = [Y12 y12];
% end
% %     plot(X12,'b','Linewidth',[1]);
% %     hold on
% %     plot(Y12,'r','Linewidth',[1]);
% %     legend('Position on x', 'Position on y','Location','Best')
% %     title('Cam1-2')
% %     xlabel('Number of Frames')
% %     ylabel('Positions')
% %     set(gca, 'Fontsize', [12])
% % end
% 
% % Cam2-2
% load ('cam2_2.mat');
% [height width rgb nf2_2] = size(vidFrames2_2);
% X22 = [];
% Y22 = [];
% % -- Watch Movie
% for j = 1:nf2_2
%     X = vidFrames2_2(:,:,:,j);
% %     imshow(X); drawnow
%     X(:,1:180) = 0;
%     X(:,430:end) = 0;
% % %     imshow(X)
%     [Max Ind] = max(X(:));
%     [y22 x22] = ind2sub(size(X), Ind);
%     X22 = [X22 x22];
%     Y22 = [Y22 y22];
% end
% % %     plot(X22,'b','Linewidth',[1]);
% % %     hold on
% % %     plot(Y22,'r','Linewidth',[1]);
% % %     legend('Position on x', 'Position on y','Location','Best')
% % %     title('Cam2-2')
% % %     xlabel('Number of Frames')
% % %     ylabel('Positions')
% % %     set(gca, 'Fontsize', [12])
% % end
% 
% % Cam3-2
% load ('cam3_2.mat');
% [height width rgb nf3_2] = size(vidFrames3_2);
% X32 = [];
% Y32 = [];
% % -- Watch Movie
% for j = 1:nf3_2
%     X = vidFrames3_2(:,:,:,j);
% %      imshow(X); drawnow
%     X(1:190,:) = 0;
%     X(340:end,:) = 0;
%     X(:,1:60) = 0;
%     X(:,500:end) = 0;
% %      imshow(X)
%     [Max Ind] = max(X(:));
%     [y32 x32] = ind2sub(size(X), Ind);
%     X32 = [X32 x32];
%     Y32 = [Y32 y32];
% %     plot(X32,'b','Linewidth',[1]);
% %     hold on
% %     plot(Y32,'r','Linewidth',[1]);
% %     legend('Position on x', 'Position on y','Location','Best')
% %     title('Cam3-2')
% %     xlabel('Number of Frames')
% %     ylabel('Positions')
% %     set(gca, 'Fontsize', [12])
% end
% % PCA
% Xmin = min([length(X12) length(X22) length(X32)]);
% X12 = X12(1:Xmin);
% X22 = X22(1:Xmin);
% X32 = X32(1:Xmin);
% Y12 = Y12(1:Xmin);
% Y22 = Y22(1:Xmin);
% Y32 = Y32(1:Xmin);
% matrix = [X12;Y12;X22;Y22;X32;Y32];
% [m,n] = size(matrix);
% % Find the mean
% miu = mean(matrix,2);
% % miu = sum(matrix)/(h*w);
% new_matrix = zeros(m,n);
% % Subtract the mean from all entries in matrix we have
% for i = 1:n
%     for j = 1:m
%         new_matrix(j,i) = matrix(j,i) - miu(j,1);
%     end
% end
% % Find covariance by the formula
% Cx = (1/(n-1)) * new_matrix * new_matrix';
% [V,D] = eig(Cx);
% lambda = diag(D);
% Y = V'*new_matrix;
% plot(Y(1,:),'Linewidth',[1])
% title('Test 2 - Noisy Case')
% xlabel('Number of Frames')
% ylabel('Positions')
% set(gca, 'Fontsize', [12])

% % Test3
% % Cam1-3
% load ('cam1_3.mat');
% [height width rgb nf3] = size(vidFrames1_3);
% X13 = [];
% Y13 = [];
% % -- Watch Movie
% for j = 1:nf3
%     X = vidFrames1_3(:,:,:,j);
% %     imshow(X); drawnow
%     X(:,1:250) = 0;
%     X(:,420:end) = 0;
% %     X(1:40,:) = 0;
% %     X(460:end, :) = 0;
% %      imshow(X)
%     [Max Ind] = max(X(:));
%     [y13 x13] = ind2sub(size(X), Ind);
%     X13 = [X13 x13];
%     Y13 = [Y13 y13];
% %     plot(X13,'b','Linewidth',[1]);
% %     hold on
% %     plot(Y13,'r','Linewidth',[1]);
% %     legend('Position on x', 'Position on y','Location','Best')
% %     title('Cam1-3')
% %     xlabel('Number of Frames')
% %     ylabel('Positions')
% %     set(gca, 'Fontsize', [12])
% end
% % % Cam2-3
% load ('cam2_3.mat');
% [height width rgb nf2_3] = size(vidFrames2_3);
% X23 = [];
% Y23 = [];
% % -- Watch Movie
% for j = 1:nf2_3
%     X = vidFrames2_3(:,:,:,j);
% %     imshow(X); drawnow
%     X(:,1:200) = 0;
%     X(:,430:end) = 0;
% %     imshow(X)
%     [Max Ind] = max(X(:));
%     [y23 x23] = ind2sub(size(X), Ind);
%     X23 = [X23 x23];
%     Y23 = [Y23 y23];
% % end
% %     plot(X23,'b','Linewidth',[1]);
% %     hold on
% %     plot(Y23,'r','Linewidth',[1]);
% %     legend('Position on x', 'Position on y','Location','Best')
% %     title('Cam2-3')
% %     xlabel('Number of Frames')
% %     ylabel('Positions')
% %     set(gca, 'Fontsize', [12])
% end
% % 
% % Cam3-3
% load ('cam3_3.mat');
% [height width rgb nf3_3] = size(vidFrames3_3);
% X33 = [];
% Y33 = [];
% % -- Watch Movie
% for j = 1:nf3_3
%     X = vidFrames3_3(:,:,:,j);
% %      imshow(X); drawnow
%     X(1:170,:) = 0;
%     X(340:end,:) = 0;
%     X(:,1:60) = 0;
%     X(:,500:end) = 0;
% %       imshow(X)
%     [Max Ind] = max(X(:));
%     [y33 x33] = ind2sub(size(X), Ind);
%     X33 = [X33 x33];
%     Y33 = [Y33 y33];
% %     plot(X33,'b','Linewidth',[1]);
% %     hold on
% %     plot(Y33,'r','Linewidth',[1]);
% %     legend('Position on x', 'Position on y','Location','Best')
% %     title('Cam3-3')
% %     xlabel('Number of Frames')
% %     ylabel('Positions')
% %     set(gca, 'Fontsize', [12])
% end
% % PCA
% Xmin = min([length(X13) length(X23) length(X33)]);
% X13 = X13(1:Xmin);
% X23 = X23(1:Xmin);
% X33 = X33(1:Xmin);
% Y13 = Y13(1:Xmin);
% Y23 = Y23(1:Xmin);
% Y33 = Y33(1:Xmin);
% matrix = [X13;Y13;X23;Y23;X33;Y33];
% [m,n] = size(matrix);
% % Find the mean
% miu = mean(matrix,2);
% new_matrix = zeros(m,n);
% % Subtract the mean from all entries in matrix we have
% for i = 1:n
%     for j = 1:m
%         new_matrix(j,i) = matrix(j,i) - miu(j,1);
%     end
% end
% % Find covariance by the formula
% Cx = (1/(n-1)) * new_matrix * new_matrix';
% [V,D] = eig(Cx);
% Y = V'*new_matrix;
% plot(Y(1,:),'Linewidth',[1])
% title('Test 3 - Horizontal Displacement')
% xlabel('Number of Frames')
% ylabel('Positions')
% set(gca, 'Fontsize', [12])

% Test4
% Cam1-4 Rotating
load ('cam1_4.mat');
[height width rgb nf1_4] = size(vidFrames1_4);
X14 = [];
Y14 = [];
% -- Watch Movie
for j = 1:nf1_4
    X = vidFrames1_4(:,:,:,j);
%     imshow(X); drawnow
    X(:,1:280) = 0;
    X(:,460:end) = 0;
%      imshow(X)
    [Max Ind] = max(X(:));
    [y14 x14] = ind2sub(size(X), Ind);
    X14 = [X14 x14];
    Y14 = [Y14 y14];
%     plot(X14,'b','Linewidth',[1]);
%     hold on
%     plot(Y14,'r','Linewidth',[1]);
%     legend('Position on x', 'Position on y','Location','Best')
%     title('Cam1-4')
%     xlabel('Number of Frames')
%     ylabel('Positions')
%     set(gca, 'Fontsize', [12])
end

% Cam2-4
load ('cam2_4.mat');
[height width rgb nf2_4] = size(vidFrames2_4);
X24 = [];
Y24 = [];
% -- Watch Movie
for j = 1:nf2_4
    X = vidFrames2_4(:,:,:,j);
%     imshow(X); drawnow
    X(:,1:200) = 0;
    X(:,430:end) = 0;
%     imshow(X)
    [Max Ind] = max(X(:));
    [y24 x24] = ind2sub(size(X), Ind);
    X24 = [X24 x24];
    Y24 = [Y24 y24];
%     plot(X24,'b','Linewidth',[1]);
%     hold on
%     plot(Y24,'r','Linewidth',[1]);
%     legend('Position on x', 'Position on y','Location','Best')
%     title('Cam2-4')
%     xlabel('Number of Frames')
%     ylabel('Positions')
%     set(gca, 'Fontsize', [12])
end
% % 
% Cam3-4
load ('cam3_4.mat');
[height width rgb nf3_4] = size(vidFrames3_4);
X34 = [];
Y34 = [];
% -- Watch Movie
for j = 1:nf3_4
    X = vidFrames3_4(:,:,:,j);
%      imshow(X); drawnow
    X(1:150,:) = 0;
    X(310:end,:) = 0;
    X(:,1:60) = 0;
    X(:,520:end) = 0;
%       imshow(X)
    [Max Ind] = max(X(:));
    [y34 x34] = ind2sub(size(X), Ind);
    X34 = [X34 x34];
    Y34 = [Y34 y34];
%     plot(X34,'b','Linewidth',[1]);
%     hold on
%     plot(Y34,'r','Linewidth',[1]);
%     legend('Position on x', 'Position on y','Location','Best')
%     title('Cam3-4')
%     xlabel('Number of Frames')
%     ylabel('Positions')
%     set(gca, 'Fontsize', [12])
end

% PCA
Xmin = min([length(X14) length(X24) length(X34)]);
X14 = X14(1:Xmin);
X24 = X24(1:Xmin);
X34 = X34(1:Xmin);
Y14 = Y14(1:Xmin);
Y24 = Y24(1:Xmin);
Y34 = Y34(1:Xmin);
matrix = [X14;Y14;X24;Y24;X34;Y34];
[m,n] = size(matrix);
% Find the mean
miu = mean(matrix,2);
new_matrix = zeros(m,n);
% Subtract the mean from all entries in matrix we have
for i = 1:n
    for j = 1:m
        new_matrix(j,i) = matrix(j,i) - miu(j,1);
    end
end
% Find covariance by the formula
Cx = (1/(n-1)) * new_matrix * new_matrix';
[V,D] = eig(Cx);
Y = V'*new_matrix;
plot(Y(1,:),'Linewidth',[1])
title('Test 4 - Horizontal Displacement and Rotation')
xlabel('Number of Frames')
ylabel('Positions')
set(gca, 'Fontsize', [12])