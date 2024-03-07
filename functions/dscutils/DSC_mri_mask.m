function mask = DSC_mri_mask(volumes, options, app)

% Computes masks for DSC-MRI examinations.
%
% Input parameters: volumes (4D matrix) containing the time course of the
% DSC signal for all voxels.
% Options is the struct containing method options, the significant ones are:
%
% options.mask.npixel: represents the minimum number of pixels of a
% connected component used as a threshold
% to exclude scalp and adjacent areas outside the brain
%
% options.display - level 1 Shows processing progress
% - level 2 Shows masks and information on thresholding
% and image intensities to be masked
%
% Output parameters: mask structure, containing
% - aif: optimized mask for arterial input function search
% - data: optimized mask for whole brain masking
% - threshold: calculated threshold provided as output

nbin=101;

volume_sum = sum(volumes,4);
mask.data = zeros(size(volume_sum));
mask.aif = zeros(size(volume_sum));

gust = volume_sum(1:options.nR*options.nC*options.nS);

[prob,intensity] = hist(gust,nbin); %#ok<HIST>

prob = prob(3:end);
intensity = intensity(3:end);

[~,ind_max]=max(prob);
if (prob(ind_max+1)/prob(ind_max)) > 0.3 
    ind_max=ind_max-1;
end

temp=prob(ind_max+1:end);
clear prob; prob=temp;

temp=intensity(ind_max+1:end);
clear intensity; intensity=temp;

f = fittype('gauss2');
g1 = inline('a1.*exp(-((x-b1)./c1).^2)'); %#ok<*DINLN>

opt_fit = fitoptions('gauss2');
gfit = fit(double(intensity)',double(prob)',f,opt_fit);

[mask.threshold,~] = curveintersect(intensity,g1(gfit.a1,gfit.b1,gfit.c1,intensity), ... 
                                  intensity,g1(gfit.a2,gfit.b2,gfit.c2,intensity));

mask.aif = volume_sum>mask.threshold;
hf_mask = zeros(options.nS,1);

for s=1:options.nS
    %copro eventuali "buchi" creati dalla sogliatura
    temp = mask.aif(:,:,s);
    
    %elimino le componenti connesse minori e lascio intatte quelle maggiori
    %di #options.mask.pixel
    
    CC = bwconncomp(temp);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    
    for j=1:size(numPixels,2)
        
        if numPixels(idx) > numPixels(j) && numPixels(j) < options.mask.npixel
            
            temp(CC.PixelIdxList{j}) = 0;
        end
    end
    
    mask.data(:,:,s)=imfill(temp,'holes');
    
    if (options.display > 2)||((options.display >1) && (s == round(0.5*options.nS)))
        hf_mask(s)=figure();
        
        subplot(121)
        imagesc(volume_sum(:,:,s))
        colormap('gray')
        title(['Slice ' num2str(s) '/' num2str(options.nS) ' - Masked Data for AIF selection'])
        
        B = bwboundaries(temp,4);
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
        end
        
        subplot(122)
        imagesc(volume_sum(:,:,s))
        
        B = bwboundaries(mask.data(:,:,s),4);
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
        end
        colormap('gray')
        title(['Slice ' num2str(s) '/' num2str(options.nS) ' - Masked Data'])
        
    end
end

mask.aif = mask.aif.*mask.data;
mask.gfit = gfit;
mask.intensity = intensity;

end

