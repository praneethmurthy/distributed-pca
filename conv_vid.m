clear
close all

vidObj = VideoReader('../modis-expts/modis_images/output_new.mp4');

vidFrames = read(vidObj);
vidlength =size(vidFrames, 4);

framesize = read(vidObj, 1);
%vidlength = vidObj.NumFrames;


for ii = 1 : vidlength -1
    tmp = double(rgb2gray(vidFrames(:,:,:,ii)));
    tmp = imresize(tmp, 0.5);
    [nrow, ncol] = size(tmp);
    L(:,ii) = reshape(tmp, [], 1);
end

% size(framesize)
% vidlength
% 
% vid1 = zeros(size(framesize, 1) * size(framesize, 2), vidlength);
% vid2 = zeros(size(framesize, 1) * size(framesize, 2), vidlength);
% vid3 = zeros(size(framesize, 1) * size(framesize, 2), vidlength);
% 
% % i = 1;
% % while hasFrame(vidObj)
% %     frame = readFrame(vidObj);
% %     
% %     i = i +1;
% % end
% 
% i=1;
% while hasFrame(vidObj)
%     tmp = readFrame(vidObj);
%     vid1(:, i) = reshape(tmp(:, :, 1), [], 1);
%     vid2(:, i) = reshape(tmp(:, :, 2), [], 1);
%     vid3(:, i) = reshape(tmp(:, :, 3), [], 1);
%     i = i+1;
% end