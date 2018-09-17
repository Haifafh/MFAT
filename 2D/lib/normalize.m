function im = normalize(im)
im = double(im); im = (im - min(im(:))) / (max(im(:)) - min(im(:)));
end