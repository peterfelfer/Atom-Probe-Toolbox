function vol = volumeDownSample(vol,factor)
% function to downsample a volume represented by a 3D matrix by a power of
% 2. The input volume will be cropped such that the JxKxL downsampled volume
% is fac x (JxKxL). The cropping is symmetric, with right side preference
% if there is a odd number of elements. The downsampling is isotropic.

sz = size(vol);
szEnd = floor(sz/factor);

% preparation of volume to fit the downsapling
removeEls = mod(size(vol),factor); % determine the number of elements to remove on each side
vol = vol(... % cropping of volume for downsampling
    1+ceil(removeEls(1)/2):end-floor(removeEls(1)/2),...
    1+ceil(removeEls(2)/2):end-floor(removeEls(2)/2),...
    1+ceil(removeEls(3)/2):end-floor(removeEls(3)/2));

% actual downsampling. Based on this post: https://de.mathworks.com/matlabcentral/answers/1766320-how-to-downsample-the-3d-matrix-in-matlab
vol = squeeze(sum(reshape(vol, factor, szEnd(1), factor, szEnd(2), factor, szEnd(3)), [1,3,5])) / factor^3;
end