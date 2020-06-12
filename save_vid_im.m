imshow(reshape(BG(:, 10)/255, imSize))
export_fig('figures/norst_10.pdf', '-pdf', '-transparent');
close;
imshow(reshape(BG(:, 100)/255, imSize))
export_fig('figures/norst_100.pdf', '-pdf', '-transparent');
close;
imshow(reshape(BG(:, 400)/255, imSize))
export_fig('figures/norst_400.pdf', '-pdf', '-transparent');
close;