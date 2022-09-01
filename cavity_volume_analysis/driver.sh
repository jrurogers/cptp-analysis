TRJCONV="gmx trjconv"
EDITCONF="gmx editconf"
# location of mdpocket
MDPOCKET="/usr/local/bin/mdpocket"

# directory with trajectories to analyze
TRJDIR="trajs"
# name of trajectory to analyze
XTC="100_300ns.xtc"
# corresponding tpr for trajectory (needed for alignment)
TPR="prod/start.tpr"
NAME=`basename $XTC | rev | cut -c5- | rev`

# reference structure to fit to
REFGRO="pocket_ref_prot.gro"
# index file of Calphas of residues 10-208 for fitting
NDX="pocket.ndx"
NDX="/home/juliarogers/Desktop/LipTransProt/cptpc1pc160_solv/05_analysis/trajs_04_md/pocket.ndx"


###########################################
################ FUNCTIONS ################
###########################################

function align () {
   if [ ! -f $REFGRO ]; then echo "$REFGRO does not exist! Create refernce protein gro! Exiting ..."; exit; fi
   if [ ! -f $NDX ]; then echo "$NDX does not exist! Create ndx for structured Ca based on $REFGRO! Exiting ..."; exit; fi

   cdir=`pwd`
   cd $TRJDIR
      # creat prot only traj
      echo 1 | $TRJCONV -f $XTC -s $TPR -o ${NAME}_prot.xtc
      # fit traj to reference gro
      echo 10 1 | $TRJCONV -f ${NAME}_prot.xtc -s $REFGRO -o ${NAME}_fit_cptpc1pc160_solv.xtc -n $NDX -fit rot+trans
   cd $cdir
}

function mdpocket () {
   # create pdb input
   if [ ! -f prot.pdb ]; then $EDITCONF -f $REFGRO -o prot.pdb; fi

   # run mdpocket
   xtc=$TRJDIR/${NAME}_fit_cptpc1pc160_solv.xtc
   echo "Running mdpocket on $xtc"
   echo ""
   $MDPOCKET --trajectory_file $xtc --trajectory_format xtc -f prot.pdb
}

function mdpocket_characterize () {
   local refpocket=$1
   local prefix=$2

   # create pdb input
   if [ ! -f prot.pdb ]; then $EDITCONF -f $REFGRO -o prot.pdb; fi

   # run mdpocket to characterize specific pocket
   xtc=$TRJDIR/${NAME}_fit_cptpc1pc160_solv.xtc
   echo "Running mdpocket on $xtc to characterize pocket $refpocket"
   echo "Outputs will have prefix $prefix"
   echo ""
   $MDPOCKET --trajectory_file $xtc --trajectory_format xtc -f prot.pdb --selected_pocket $refpocket -o $prefix
}

###########################################
################### MAIN ##################
###########################################

# align traj to CPTP+C1P solvated ref
align

# run mdpocket
mdpocket

# run mdpocket to analyze specific pocket
REFPOCKET="c1p_pocket_0_2.pdb"
mdpocket_characterize $REFPOCKET "c1p"


