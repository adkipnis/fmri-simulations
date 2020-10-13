%% Safety checks
% Check how the AR(1) model holds
voxel = randi(size(y,2))
hold on
legend
plot(y(:,voxel))
plot(y_ar(:,voxel))
plot(eps(:,voxel))
hold off

% Check how permutation changes the timeseries of one voxel
hold on
legend
plot(y(:,voxel))
plot(y_perm(:,voxel))
hold off


a = sum(sum(y));
b = sum(sum(y_ar));
c = sum(sum(eps));
d = sum(sum(y_perm));
e = sum(sum(eps_perm));

% Check that the total residual amplitudes stay unchanged 
% 1. amp(y) = amp(y_pred + eps)
a
b+c

% 2. amp(y) = amp(y_perm)
a
d 

% 3. % amp(eps) = amp(eps_perm) (this assumes that n_permutations = 1)
c
e


% Plot residuals before and after permutations
tp = randi(178) % for the first tp the two should look identical

X1 = E(:,:,:,tp);
X2 = E_perm(:,:,:,tp,1);

sl = randi(size(E,3)) %axial slices
hold on
subplot(2,1,1)
imagesc(X1(:,:,sl))
subplot(2,1,2)
imagesc(X2(:,:,sl))
hold off

nifti_utils.vol_viewer_4D(E);
nifti_utils.vol_viewer_4D(E_perm(:,:,:,:,1));

% Check that the permutations don't change the total noise amplitude
sum(sum(sum(nansum(E))))
sum(sum(sum(nansum(E_perm(:,:,:,:,1)))))

% Add Signal and permuted noise
S = repaste_to_5d(res_mat, y_glm, mask_vec, Opts);
D = S + E_perm;