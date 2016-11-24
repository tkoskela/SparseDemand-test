# gethost jid hostfile
#   $1 = jid is the job id of the job 
#   $2 = hostfile is the name of the hostfile to produce
#
# This script uses qstat to produce a hostfile for an sge job.
# the script calls qstat, and lists all the hosts that are assigned to jid
# the list of hosts is output to hostfile.
qstat -u uctpln0 | awk '($1 =='$1') {print $8}' > $2
awk -F@ '{print $2}' $2 > temphosts.txt
awk -F. '{print $1}' temphosts.txt > $2
rm temphosts.txt
