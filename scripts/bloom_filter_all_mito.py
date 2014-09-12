import os
import sys
import shutil
from Bio import SeqIO

"""
Submit to cluster as follow,
 -V
    means inherit all the envionment variables (so finds my Python libraries)
 -t "1-20:1"
    means run 20 sub tasks (there are 20 files to produce)
 -pe smp 4
    means each task wants 4 CPUs (really I just want to grab a machine in case of RAM shortages
    and because this will be IO heavy)


$ qsub -V -t "1-20:1" -pe smp 4 bloom_filter_all_mito.sh
Your job-array 457.1-20:1 ("bloom_filter_all_mito.sh") has been submitted


$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
    457 0.00000 bloom_filt pc40583      qw    10/01/2012 17:30:18                                    4 1-20:1



$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-wash              4 1
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-kaylee            4 2
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-simon             4 3
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-river             4 4
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-mal               4 5
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-jayne             4 6
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n09-08-144-biggus            4 7
    457 0.55500 bloom_filt pc40583      t     10/01/2012 17:30:32 all.q@n09-08-144-biggus            4 8
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n11-04-048-cortana           4 9
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-002-baldrick          4 10
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-zoe               4 11
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-book              4 12
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-008-inara             4 13
    457 0.55500 bloom_filt pc40583      r     10/01/2012 17:30:32 all.q@n08-04-002-blackadder        4 14
    457 0.00000 bloom_filt pc40583      qw    10/01/2012 17:30:18                                    4 15-20:1
"""

host = os.environ["HOSTNAME"]

if "SGE_TASK_ID" in os.environ:
    job = int(os.environ["SGE_TASK_ID"])
    print "Doing task %i, this is %s" % (job, host)
else:
    #Assuming doing all of them:
    job = None #All of them
    #job = -1 #None of them

ref = "all_mito2.fas"
python = "python2.7"
cmd_template = python + " interlace_fastq.py %s %s | " + python + " blooming_reads.py -c %s -k 20 -m 1 -f fastq | bgzip > %s"
tmp_template = "/tmp/bloom_filter_all_mito_k20_1sid_rc/%s"
new_template = "../data/bloom_filter_all_mito_k20_1sid_rc/%s"
raw_fastq_folder = "../raw_data"

if not os.path.isdir(new_template % ""):
    os.mkdir(new_template % "")

if not os.path.isdir(tmp_template % ""):
    os.mkdir(tmp_template % "")

def run(cmd):
    print cmd
    rc = os.system(cmd)
    if rc:
        sys.stderr.write("Return code %i from running command:\n%s\n" % (rc, cmd))
        sys.exit(rc)

this_job = 0
for root, dirs, files in os.walk(raw_fastq_folder):
    for f in files:
        if f.endswith("_1.fastq.gz"):
            this_job += 1
            base = f[:-11]
            old1 = os.path.join(root, base + "_1.fastq.gz")
            assert os.path.isfile(old1)
            old2 = os.path.join(root, base + "_2.fastq.gz")
            assert os.path.isfile(old2)
            new = new_template % base + ".fastq.bgz"
            tmp= tmp_template % base + ".fastq.bgz"
            idx = new + ".idx"

            if job and job != this_job:
                print "%i: %s (SKIPPING, NOT ASSIGNED)" % (this_job, new)
                continue
            if not os.path.isfile(new):
                if os.path.isfile(new[:-4] + ".gz"):
                    print "%i: %s (RECOMPRESSING...)" % (this_job, new)
                    cmd = "cat %s | gunzip | bgzip > %s" % (new[:-4] + ".gz", tmp)
                else:
                    print "%i: %s (RUNNING...)" % (this_job, new)
                    cmd = cmd_template % (old1, old2, ref, tmp)
                run(cmd)
                if not os.path.isfile(tmp):
                    print "%s (FAILED), missing:\n%s" % (new, tmp)
                    continue
                elif not os.path.getsize(tmp):
                    print "%s (FAILED), empty:\n%s" % (new, tmp)
                    os.remove(tmp)
                    continue
                shutil.copyfile(tmp, new)
                os.remove(tmp)
            if os.path.isfile(idx + "-journal"):
                #print "%i: %s (ABORTING - SQLite3 journal file found!)" % (this_job, new)
                sys.stderr.write("%i: %s (ABORTING - SQLite3 journal file found!)\n" % (this_job, new))
                continue
            if not os.path.isfile(idx):
                print "%i: %s (INDEXING...)" % (this_job, new)
            else:
                print "%i: %s (CHECKING INDEX...)" % (this_job, new)
            #Build index if not there, load it if it is there:
            d = SeqIO.index_db(idx, new, "fastq")
            print "%i: %s has %i entries" % (this_job, new, len(d))
            assert d, "Empty dict?"
