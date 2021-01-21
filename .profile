# Profile to be sourced before running tesser scripts.

case $USER in
    morton)
        STUDYDIR=$HOME/Dropbox/work/tesser
        conda activate tesser
        ;;

    mortonne)
        STUDYDIR=$STOCKYARD/lonestar/tesser
        unset PYTHONPATH  # remove global site-packages directory
        . $STOCKYARD/lonestar/venv/tesser/bin/activate
        export BATCHDIR=$STOCKYARD/lonestar/tesser/batch/launchscripts
        ;;

    *)
        echo "Error: unknown user $USER."
        ;;
esac
export STUDYDIR