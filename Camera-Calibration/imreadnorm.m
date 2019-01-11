function I = imreadnorm(file)

I = im2double(imread(file));

if (size(I,3) > 1)
  I = rgb2gray(I);
end

I = I - min(I(:));
I = I / max(I(:));
