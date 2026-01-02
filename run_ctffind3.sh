#!/bin/bash
NCPUS=4        # threads per micrograph
JOBS=6         # micrographs in parallel
CTF=/data/csun/softwares/ctf/ctffind3_mp.exe

run_one () {
  mrc="$1"
  log="${mrc%.mrc}_ctffind3.log"

  [[ -f "$log" ]] && { echo "Skipping $mrc"; return 0; }

  export NCPUS="$NCPUS"
  "$CTF" << EOF > "$log" 2>&1
$mrc
${mrc%.mrc}_ctf.mrc
2.7,300.0,-0.1,28902,2.5
512,50,5,50000,5000,500,900
EOF
  echo "Done $mrc"
}

export -f run_one
export NCPUS CTF

parallel -j "$JOBS" run_one ::: *patch_aligned.mrc
