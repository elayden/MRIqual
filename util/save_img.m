function save_img(use_header,img,fpath)
    if use_header
        try save_nii(img,fpath)
        catch; save_untouch_nii(img,fpath)
        end
    else save_untouch_nii(img,fpath)
    end
end