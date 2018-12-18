

%-Wiktor Olszowy


path_manage = fgetl(fopen('path_manage.txt'));

addpath('/applications/spm/spm12_7219');

cd([path_manage '/brain_parcellation']);

P = {'mask.nii'; 'labels_Neuromorphometrics.nii'};

%-for Peter Zeidman's routine (from SPM's mailing list), interpolation with splines is fine, as only one region is resliced, but not the entire brain!
%-I do not employ any splines, just nearest neighbour (NN)

spm_reslice(P, struct('which', 1, 'interp', 0));

VOIs  = [200 201 108 109 182 183 192 193];
%-pauses needed as otherwise I think next command can be started without finishing the former one
for i = 1:length(VOIs)
   system(['fslmaths rlabels_Neuromorphometrics    -thr  ' num2str(VOIs(i)) ' mask_VOI_' num2str(VOIs(i))]);
   pause(5);
   system(['fslmaths mask_VOI_' num2str(VOIs(i)) ' -uthr ' num2str(VOIs(i)) ' -bin mask_VOI_' num2str(VOIs(i))]);
   pause(5);
   system(['gunzip mask_VOI_'   num2str(VOIs(i)) '.nii.gz']);
   pause(5);
end

cd(path_manage)
