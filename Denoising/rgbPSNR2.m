function psnrvalue = rgbPSNR2(image1, image2)
% image1和image2大小相等
image1 = double(image1);
image2 = double(image2);

[row, col, ~] = size(image1); % 图像的长和宽
MSE = sum(bsxfun(@minus, image1, image2).^2, 'all') / (row * col * 3); % 计算MSE

MAX = double(intmax(class(image1))); % 图像的灰度级数
B = ceil(log2(MAX + 1)); % 编码一个像素所用二进制位数
psnrvalue = 20 * log10(MAX) - 10 * log10(MSE); % 计算PSNR
end
