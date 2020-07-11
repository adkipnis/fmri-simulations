import os
import pickle
import numpy as np
import mri_config as config
import bdpy
from itertools import product

import brainiak



def main():
    # Settings ---------------------------------------------------------

    # Data settings
    subjects = config.subjects
    rois = config.rois
    num_voxel = config.num_voxel


    #n_iter = 200

    #results_dir = config.results_dir


    # Load data --------------------------------------------------------
    print('----------------------------------------')
    print('Loading data')

    data_all = {}
    for sbj in subjects:
        if len(subjects[sbj]) == 1:
            data_all[sbj] = bdpy.BData(subjects[sbj][0])
        else:
            # Concatenate data
            suc_cols = ['Run', 'Block']
            data_all[sbj] = concat_dataset([bdpy.BData(f) for f in subjects[sbj]],
                                           successive=suc_cols)


    for sbj, roi, feat in product(subjects, rois, features):
        print('--------------------')
        print('Subject:    %s' % sbj)
        print('ROI:        %s' % roi)
        print('Num voxels: %d' % num_voxel[roi])
        print('Feature:    %s' % feat)


        # Prepare data
        print('Preparing data')
        dat = data_all[sbj]

        x = dat.select(rois[roi])           # Brain data
        datatype = dat.select('DataType')   # Data type
        labels = dat.select('stimulus_id')  # Image labels in brain data

        y = data_feature.select(feat)             # Image features
        y_label = data_feature.select('ImageID')  # Image labels

        # For quick demo, reduce the number of units from 1000 to 100
        y = y[:, :100]

        y_sorted = get_refdata(y, y_label, labels)  # Image features corresponding to brain data

        # Get training and test dataset
        i_train = (datatype == 1).flatten()    # Index for training
        i_test_pt = (datatype == 2).flatten()  # Index for perception test
        i_test_im = (datatype == 3).flatten()  # Index for imagery test
        i_test = i_test_pt + i_test_im

        x_train = x[i_train, :]
        x_test = x[i_test, :]

        y_train = y_sorted[i_train, :]
        y_test = y_sorted[i_test, :]

        # Feature prediction
        pred_y, true_y = feature_prediction(x_train, y_train,
                                            x_test, y_test,
                                            n_voxel=num_voxel[roi],
                                            n_iter=n_iter)

        # Separate results for perception and imagery tests
        i_pt = i_test_pt[i_test]  # Index for perception test within test
        i_im = i_test_im[i_test]  # Index for imagery test within test

        pred_y_pt = pred_y[i_pt, :]
        pred_y_im = pred_y[i_im, :]

        true_y_pt = true_y[i_pt, :]
        true_y_im = true_y[i_im, :]

        # Get averaged predicted feature
        test_label_pt = labels[i_test_pt, :].flatten()
        test_label_im = labels[i_test_im, :].flatten()

        pred_y_pt_av, true_y_pt_av, test_label_set_pt \
            = get_averaged_feature(pred_y_pt, true_y_pt, test_label_pt)
        pred_y_im_av, true_y_im_av, test_label_set_im \
            = get_averaged_feature(pred_y_im, true_y_im, test_label_im)

        # Get category averaged features
        catlabels_pt = np.vstack([int(n) for n in test_label_pt])  # Category labels (perception test)
        catlabels_im = np.vstack([int(n) for n in test_label_im])  # Category labels (imagery test)
        catlabels_set_pt = np.unique(catlabels_pt)                 # Category label set (perception test)
        catlabels_set_im = np.unique(catlabels_im)                 # Category label set (imagery test)

        y_catlabels = data_feature.select('CatID')   # Category labels in image features
        ind_catave = (data_feature.select('FeatureType') == 3).flatten()

        y_catave_pt = get_refdata(y[ind_catave, :], y_catlabels[ind_catave, :], catlabels_set_pt)
        y_catave_im = get_refdata(y[ind_catave, :], y_catlabels[ind_catave, :], catlabels_set_im)

        # Prepare result dataframe
        results = pd.DataFrame({'subject' : [sbj, sbj],
                                'roi' : [roi, roi],
                                'feature' : [feat, feat],
                                'test_type' : ['perception', 'imagery'],
                                'true_feature': [true_y_pt, true_y_im],
                                'predicted_feature': [pred_y_pt, pred_y_im],
                                'test_label' : [test_label_pt, test_label_im],
                                'test_label_set' : [test_label_set_pt, test_label_set_im],
                                'true_feature_averaged' : [true_y_pt_av, true_y_im_av],
                                'predicted_feature_averaged' : [pred_y_pt_av, pred_y_im_av],
                                'category_label_set' : [catlabels_set_pt, catlabels_set_im],
                                'category_feature_averaged' : [y_catave_pt, y_catave_im]})

        # Save results
        makedir_ifnot(os.path.dirname(results_file))
        with open(results_file, 'wb') as f:
            pickle.dump(results, f)

        print('Saved %s' % results_file)

        dist.unlock()
