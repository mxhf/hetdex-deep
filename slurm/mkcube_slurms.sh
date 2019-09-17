cat ifulist_COSMOSA.txt | awk '{print "COSMOSA", $1}' >> worklist
cat ifulist_COSMOSB.txt | awk '{print "COSMOSB", $1}' >> worklist
cat ifulist_COSMOSC.txt | awk '{print "COSMOSC", $1}' >> worklist
cat ifulist_COSMOSD.txt | awk '{print "COSMOSD", $1}' >> worklist
cat ifulist_GOODSN.txt | awk '{print "GOODSN", $1}' >> worklist
python mkcube_slurms.py --dir_rebin=/scratch/04287/mxhf/rebin2 -i worklist 
