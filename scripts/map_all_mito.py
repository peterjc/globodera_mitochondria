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


$ qsub -V -t "1-20:1" -pe smp 4 map_all_mito.sh
Your job-array 463.1-20:1 ("map_all_mito.sh") has been submitted


$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
    463 0.00000 map_all_mi pc40583      qw    10/02/2012 13:31:20                                    4 1-20:1


$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n09-08-144-biggus            4 2
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n11-04-048-cortana           4 4
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-simon             4 5
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-jayne             4 6
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-wash              4 7
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-book              4 8
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-kaylee            4 9
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-inara             4 10
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-river             4 11
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-002-baldrick          4 12
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-008-mal               4 13
    463 0.55500 map_all_mi pc40583      r     10/02/2012 13:31:31 all.q@n08-04-002-blackadder        4 14
    463 0.00000 map_all_mi pc40583      qw    10/02/2012 13:31:20 


When/if need more RAM,

$ qsub -l hostname="n11-04-048-cortana|n09-08-144-biggus"  -V -t "1-20:1" -pe smp 4 map_all_mito.sh
...

$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
    464 0.55500 map_all_mi pc40583      r     10/02/2012 17:33:16 all.q@n11-04-048-cortana           4 10
    464 0.55500 map_all_mi pc40583      r     10/02/2012 17:33:16 all.q@n09-08-144-biggus            4 11
    464 0.55500 map_all_mi pc40583      r     10/02/2012 17:33:16 all.q@n09-08-144-biggus            4 12
    464 0.55500 map_all_mi pc40583      qw    10/02/2012 17:32:29                                    4 13-20:1

Note since it has 8 cores, two jobs are being run at once on biggus.

"""

host = os.environ["HOSTNAME"]

if "SGE_TASK_ID" in os.environ:
    job = int(os.environ["SGE_TASK_ID"])
    print "Doing task %i, this is %s" % (job, host)
else:
    #Assuming doing all of them:
    job = None #Use None for all of them
    #job = -1 #Use -1 for none of them

double_ref = "all_mito2.fas"
single_ref = "all_mito.fas"
folder = "bloom_filter_all_mito_k20_1sid_rc"
dedup_circular_sam = "~/repositories/picobio/blooming_reads/dedup_circular_sam.py"
re_pair_circular_sam = "~/repositories/picobio/blooming_reads/re_pair_circular_sam.py"

cmd_template = "mrfast --seq %s --seqcomp --search %s -o %s"
dedup_template = dedup_circular_sam + " -i %s -c %s | samtools view -u -S - > %s"
bam_queryname_sort_template = "samtools sort -n %s %s"

repair_template = "samtools view -h %s | " + re_pair_circular_sam + " -r %s -c %s -v %s | samtools view -u -S - > %s"
bam_coord_sort_template = "samtools sort %s %s"

tmp_template = "/tmp/"  + folder + "/%s"
new_template = folder + "/%s"

if not os.path.isdir(tmp_template % ""):
    os.mkdir(tmp_template % "")

def run(cmd):
    print cmd
    rc = os.system(cmd)
    if rc:
        sys.stderr.write("Return code %i from running command:\n%s\n" % (rc, cmd))
        sys.exit(rc)

this_job = 0
single_ref = "all_mito.fas"
for f in sorted(os.listdir(folder)):
    if f.endswith(".fastq.bgz"):
        this_job += 1
        base = f[:-10]
        old = os.path.join(folder, f)
        new = new_template % base + "_vs_all_mito2_mrfast.sam"
        dedup = new_template % base + "_vs_all_mito_mrfast_dedup.bam"
        repair = new_template % base + "_vs_all_mito_mrfast_repair.bam"
        cov = new_template % base + "_vs_all_mito_mrfast_repair.cov"

        reads = new_template % base + ".fastq.bgz"

        if job and job != this_job:
            print "%i: %s (SKIPPING, NOT ASSIGNED)" % (this_job, new)
            continue
        if not os.path.isfile(new) and not os.path.isfile(dedup):
            print "%i: %s (MAPPING...)" % (this_job, new)
            tmp = tmp_template % base + "_vs_all_mito2_mrfast.sam"
            cmd = cmd_template % (old, double_ref, tmp)
            #print "Meh... %s" % cmd
            #continue
            run(cmd)
            if not os.path.isfile(tmp):
                print "%i: %s (MAPPING FAILED), missing:\n%s" % (this_job, new, tmp)
                continue
            if not os.path.getsize(tmp):
                print "%i: %s (MAPPING FAILED), empty output:\n%s" % (this_job, new, tmp)
                continue
            shutil.copyfile(tmp, new)
            os.remove(tmp)
        if os.path.isfile(new) and not os.path.isfile(dedup):
            print "%i: %s (DEDUP...)" % (this_job, dedup)
            tmp = tmp_template % base + "_vs_all_mito_mrfast_dedup_unsorted.bam"
            if os.path.isfile(tmp):
                os.remove(tmp)
            cmd = dedup_template % (new, single_ref, tmp)
            run(cmd)
            if not os.path.isfile(tmp):
                print "%i: %s (DEDUP FAILED), missing:\n%s" % (this_job, dedup, tmp)
                continue
            if not os.path.getsize(tmp):
                print "%i: %s (DEDUP FAILED), empty output:\n%s" % (this_job, dedup, tmp)
                continue
            sorted_prefix = tmp_template % base + "_vs_all_mito_mrfast_dedup"
            sorted_bam = sorted_prefix + ".bam"
            if os.path.isfile(sorted_bam):
                os.remove(sorted_bam)
            cmd = bam_queryname_sort_template % (tmp, sorted_prefix)
            run(cmd)
            if not os.path.isfile(sorted_bam):
                print "%i: %s (SORT FAILED), missing:\n%s" % (this_job, dedup, sorted_bam)
                continue
            if not os.path.getsize(sorted_bam):
                print "%i: %s (SORT FAILED), empty output:\n%s" % (this_job, dedup, sorted_bam)
                continue
            shutil.copyfile(sorted_bam, dedup)
            os.remove(sorted_bam)
            os.remove(tmp)
        if os.path.isfile(dedup) and not os.path.isfile(repair) or not os.path.isfile(cov):
            print "%i: %s (re-pair...)" % (this_job, repair)
            tmp = tmp_template % base + "_vs_all_mito_mrfast_repair_unsorted.bam"
            cov_tmp = tmp_template % base + "_vs_all_mito_mrfast_repair.cov"
            cmd = repair_template % (dedup, reads, single_ref, cov_tmp, tmp)
            run(cmd)
            if not os.path.isfile(cov_tmp):
                print "%i: %s (RE-PAIR FAILED), missing:\n%s" % (this_job, repair, cov_tmp)
                continue
            if not os.path.getsize(cov_tmp):
                print "%i: %s (RE-PAIR FAILED), empty output:\n%s" % (this_job, repair, cov_tmp)
                continue
            if not os.path.isfile(tmp):
                print "%i: %s (RE-PAIR FAILED), missing:\n%s" % (this_job, repair, tmp)
                continue
            if not os.path.getsize(tmp):
                print "%i: %s (RE-PAIR FAILED), empty output:\n%s" % (this_job, repair, tmp)
                continue
            sorted_prefix = tmp_template % base + "_vs_all_mito_mrfast_repair"
            sorted_bam = sorted_prefix + ".bam"
            if os.path.isfile(sorted_bam):
                os.remove(sorted_bam)
            cmd = bam_coord_sort_template % (tmp, sorted_prefix)
            run(cmd)
            if not os.path.isfile(sorted_bam):
                print "%i: %s (SORT FAILED), missing:\n%s" % (this_job, dedup, sorted_bam)
                continue
            if not os.path.getsize(sorted_bam):
                print "%i: %s (SORT FAILED), empty output:\n%s" % (this_job, dedup, sorted_bam)
                continue
            shutil.copyfile(sorted_bam, repair)
            os.remove(sorted_bam)
            os.remove(tmp)
            shutil.copyfile(cov_tmp, cov)
            os.remove(cov_tmp)
        print "%i: %s (DONE)" % (this_job, dedup)
