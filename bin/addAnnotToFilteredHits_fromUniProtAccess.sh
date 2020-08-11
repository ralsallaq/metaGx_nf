module load parallel
module load python/2.7.15-rhel7
input=`sed 1d hits_in_1729716_EW1H-0179_S86_L001.csv`
#echo "$input" | parallel --block 100k --pipe python splitFile.py 
echo "$input" | parallel  --block 10k --pipe python addAnnotToFilteredHits_fromUniProtAccess.py
