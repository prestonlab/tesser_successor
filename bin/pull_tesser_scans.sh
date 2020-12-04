#!/bin/bash
#
# Pull tesser scan data for local analysis.

dest="$1"
shift 1

rsync -azvu lonestar:/corral-repl/utexas/prestonlab/tesser/ "$dest" \
    --include="batch/" \
    --include="batch/analysis/" \
    --include="batch/analysis/rsa_beta/" \
    --include="batch/analysis/rsa_beta/rsa_event_textfiles/" \
    --include="tesser_*_run?_info.txt" \
    --include="tesser_*/" \
    --include="tesser_*/BOLD/" \
    --include="tesser_*/BOLD/functional_run_?/" \
    --include="tesser_*/BOLD/functional_run_?/QA/" \
    --include="confound.txt" \
    --include="tesser_*/BOLD/antsreg/" \
    --include="tesser_*/BOLD/antsreg/data/" \
    --include="tesser_*/BOLD/antsreg/data/functional_run_?_bold_mcf_brain_corr_notemp.feat/" \
    --include="filtered_func_data.nii.gz" \
    --include="tesser_*/anatomy/" \
    --include="tesser_*/anatomy/antsreg/" \
    --include="tesser_*/anatomy/antsreg/data/" \
    --include="tesser_*/anatomy/antsreg/data/funcunwarpspace/" \
    --include="tesser_*/anatomy/antsreg/data/funcunwarpspace/rois/" \
    --include="tesser_*/anatomy/antsreg/data/funcunwarpspace/rois/mni/" \
    --include="b_hip*.nii.gz" \
    --exclude="*" "$@"
