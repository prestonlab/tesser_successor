"""Dataset handling and measures for use with pymvpa2."""

import os
import numpy as np
from mvpa2.datasets.mri import fmri_dataset
from mvpa2.base.dataset import vstack
from mvpa2.mappers.zscore import zscore


def load_struct_timeseries(
    study_dir, subject, mask, feature_mask=None, verbose=1, zscore_run=True
):
    """Load structure learning timeseries as a dataset."""
    # functional timeseries for each run
    subject_dir = os.path.join(study_dir, f'tesser_{subject}')
    runs = list(range(1, 7))
    bold_images = []
    for run in runs:
        bold = os.path.join(
            subject_dir, 'BOLD', 'antsreg', 'data',
            f'functional_run_{run}_bold_mcf_brain_corr_notemp.feat',
            'filtered_func_data.nii.gz'
        )
        if not os.path.exists(bold):
            raise IOError(f'BOLD file does not exist: {bold}')
        bold_images.append(bold)

    # mask image to select voxels to load
    mask_dir = os.path.join(
        subject_dir, 'anatomy', 'antsreg', 'data', 'funcunwarpspace', 'rois', 'mni'
    )
    mask_file = os.path.join(mask_dir, f'{mask}.nii.gz')
    if not os.path.exists(mask_file):
        raise IOError(f'Mask file does not exist: {mask_file}')
    if verbose:
        print(f'Masking with: {mask_file}')

    # feature mask, if specified
    if feature_mask is not None:
        feature_file = os.path.join(mask_dir, 'f{feature_mask}.nii.gz')
        if not os.path.exists(feature_file):
            raise IOError(f'Feature mask does not exist: {feature_file}')
        if verbose:
            print(f'Using features within: {feature_file}')
        add_fa = {'include': feature_file}
    else:
        add_fa = None

    # load images and concatenate
    ds_list = []
    for run, bold_image in zip(runs, bold_images):
        ds_run = fmri_dataset(bold_image, mask=mask_file, add_fa=add_fa)
        ds_run.sa['run'] = run
        ds_list.append(ds_run)
    ds = vstack(ds_list)

    # normalize within run
    if zscore_run:
        zscore(ds, chunks_attr='run')

    # set the include feature attribute
    if feature_mask is not None:
        ds.fa['include'] = ds.fa['include'].astype(bool)
    else:
        ds.fa['include'] = np.ones(ds.shape[1], dtype=bool)
    return ds
