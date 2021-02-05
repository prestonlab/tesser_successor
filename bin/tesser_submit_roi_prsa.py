#!/usr/bin/env python
#
# Submit partial representational similarity analysis.

import os
import argparse
from tesser import util


def main(subjects, rois, study_dir, res_name, block):
    if subjects is None:
        subjects = util.subj_list()

    beh_dir = os.path.join(study_dir, 'batch', 'behav')
    options = f'--study-dir={study_dir} -b {block}'
    for roi in rois:
        for subject in subjects:
            print(f'tesser_roi_prsa.py {subject} {beh_dir} {roi} {res_name} {options}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('rois', help="name of mask to use.")
    parser.add_argument('res_name', help="name of directory to save results.")
    parser.add_argument('--block', '-b', help='block to include (walk, random)')
    parser.add_argument('--study-dir', help="path to main study data directory.")
    parser.add_argument('--subjects', '-s', help="ID of subjects to process.")
    args = parser.parse_args()

    if args.study_dir is None:
        if 'STUDYDIR' not in os.environ:
            raise RuntimeError('STUDYDIR environment variable not set.')
        env_study_dir = os.environ['STUDYDIR']
    else:
        env_study_dir = args.study_dir

    if args.subjects is not None:
        inc_subjects = args.subjects.split(',')
    else:
        inc_subjects = None

    inc_rois = args.rois.split(',')
    main(inc_subjects, inc_rois, env_study_dir, args.res_name, block=args.block)
