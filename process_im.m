
% code called by "label_oatsarenodes_graph.m" to label plasmodium in images (1=plasmodium, 0=background)

function [im_bw,im_skel] = process_im(im, strel_size, THRESHOLD)

    sub_bg = imtophat(im, strel('disk', strel_size));
    
    im_bw = sub_bg > mean(sub_bg(sub_bg>0)) + THRESHOLD*std(sub_bg(sub_bg>0));
    im_skel = bwskel(logical(im_bw));
    
end
