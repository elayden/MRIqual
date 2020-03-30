function [roidat,numrois,roi_dim] = load_ROI(roi,background_dim,input_type)
    if ~isempty(roi)
        switch class(roi)
            case 'char'
                try 
                    roi_img = load_nii(roi); 
                catch
                    try
                        roi_img = load_untouch_nii(roi);
                        warning(['Non-orthogonal shearing detected in affine matrix of ',...
                            input_type,' image. Loaded successfully without applying affine.'])
                    catch
                        error(['Error:  failed to load ',input_type,' image'])
                    end
                end
                roidat = roi_img.img; roi_dim = size(roidat); numrois = max(roidat(:));
            case 'struct'
                try
                    roidat = roi.img;
                    roi_dim = size(roidat); 
                    numrois = max(roidat(:));
                catch
                    error(['If input ',input_type,' is a structure, it should contain a field ''.img''.'])
                end
            case 'cell'
                try
                    roidat = zeros(background_dim);
                    for i = 1:length(roi)
                        for j = 1:size(roi{i},1)
                            roidat(roi{i}(j,1),roi{i}(j,2),roi{i}(j,3)) = i;
                        end
                    end
                catch
                    error(['Check ',input_type,' cell array input.'])
                end
                numrois = length(roi); roi_dim = background_dim;
            otherwise
                if isnumeric(roi) && ~isempty(roi)
                    roi_dim = size(roi);
                    if all(roi_dim==background_dim)
                        roidat = roi; 
                        numrois = max(roidat(:));
                        if numrois==0
                            error(['Input ',input_type,' matrix contains no non-zero entries.'])
                        end
                    else
                        error(['Input ',input_type,' matrix dimensions do not match background.'])
                    end
                else numrois = 0; roidat = 0;
                end
        end
    else numrois = 0; roi_dim = []; roidat = [];
    end
    if ~all(roi_dim==background_dim)
        error(['Input ',input_type,' matrix dimensions do not match background.'])
    end
end