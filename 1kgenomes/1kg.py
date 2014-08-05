#!/usr/bin/env python

"""
1kg.py

Author: Ben Langmead
Date: 3/18/2012
Contact: blangmea@jhsph.edu

Various functionality related to downloading and summarizing reads from the
1000 Genomes Project.

Stratifies results in several ways:
1. By read length and by SAM flag
2. By input file and by SAM flag

Assumes that 'sequence.index' is present in current directory.  Use
get_indexes.sh to get that file.

TODO:
- Paired-end alignment
"""

# python 1kg.py 1KG_P2_NA12878_GAII.manifest 4 /home/langmead/svn/bowtie2_trunk/bowtie2 -x $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa -p 7 --sample 0.01
# python 1kg.py 1KG_ASW_HISEQ.manifest       4 /home/langmead/svn/bowtie2_trunk/bowtie2 -x $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa -p 7 --sample 0.01

import sys
import os
import string
import Queue
import subprocess
import threading
import gzip
import time
import errno
import StringIO
import argparse

manifest       = sys.argv[1]      # manifest file with inputs
ncpus          = int(sys.argv[2]) # number of cpus to use per bowtie2 process
bt2_exe        = sys.argv[3]      # path to bowtie2 binary
bt2_args       = sys.argv[4:]     # arguments for bowtie2
outdir         = 'out'            # subdirectory where all results go
download_first = True             # download all the reads first
keep_reads     = True             # don't delete reads after aligning them
skip_existing  = True             # if sam file exists, skip

queue = Queue.Queue()
glock = threading.Lock()

class Total(object):
    def __init__(self):
        self.n = 0

class AlnSummary(object):
    ''' Summary statistics for one stratum of reads '''
    def __init__(self):
        self.quals_by_cyc = []  # Quality values by cycle
        self.avg_qual     = []  # Average per-base quality value
        self.median_qual  = []  # Median quality value
        self.nns          = []  # Number of Ns
        self.nupdate      = 0   # Number times new read was added
        self.readbuf      = StringIO.StringIO()

def median(l):
    assert len(l) > 0
    off = int(len(l)/2)
    return sorted(l)[off]

def mean(l):
    assert len(l) > 0
    return round(sum(map(lambda x: float(x) / len(l), l)))

def ensure(l, off):
    ''' Ensure list of ints l is long enough to permit access to element at
        offset off '''
    if len(l) <= off:
        diff = off + 1 - len(l)
        l.extend([0] * diff)

def ensureList(l, off):
    ''' Ensure list of lists l is long enough to permit access to element at
        offset off '''
    if len(l) <= off:
        diff = off + 1 - len(l)
        l.extend([[]] * diff)
        assert len(l) >= off+1

summs = dict()
tot = Total()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

class Worker(threading.Thread):
    ''' Run bowtie2 on a sample of the reads.  Collect data on how many reads
        align and how many fail to align.  Also, collect data on the quality
        values in the reads that align and that fail to align. '''
    def __init__(self, tid, queue, glock, summs, tot, outdir):
        threading.Thread.__init__(self)
        self.queue = queue
        self.glock = glock
        self.summs = summs  # overall summary, stratified by flag and length
        self.n = 0
        self.nglob = tot
        self.tid = tid
        self.outdir = outdir
    
    def update(self, st, ln):
        ''' Update counters for given stratum summary '''
        ln = ln.rstrip()
        ts = string.split(ln, '\t')
        rdname, flag, seq, qual = ts[0], int(ts[1]), ts[9], ts[10]
        nns = ts[9].count('N') + ts[9].count('.')
        qs = map(lambda x: ord(x)-33, list(qual))
        mea = int(mean(qs))
        med = int(median(qs))
        k = (flag, len(seq))
        if k not in st: st[k] = AlnSummary()
        ensureList(st[k].quals_by_cyc, len(qs))
        assert len(st[k].quals_by_cyc) >= len(qs)+1
        for i in xrange(0, len(qs)):
            q = qs[i]
            ensure(st[k].quals_by_cyc[i], q)
            st[k].quals_by_cyc[i][q] += 1
        ensure(st[k].avg_qual, mea)
        st[k].avg_qual[mea] += 1
        ensure(st[k].median_qual, med)
        st[k].median_qual[med] += 1
        ensure(st[k].nns, nns)
        st[k].nns[nns] += 1
        st[k].nupdate += 1
        st[k].readbuf.write('@%s\n%s\n+\n%s\n' % (rdname, seq, qual))
    
    def flush(self, st, nupdate, outd):
        ''' Write out info from the given stratum summary to the appropriate
            output directory '''
        str = ""
        mkdir_p(outd)
        ofh_basics  = open(os.path.join(outd, 'basics.tsv'),      'w')
        ofh_medians = open(os.path.join(outd, 'medians.tsv'),     'w')
        ofh_means   = open(os.path.join(outd, 'means.tsv'),       'w')
        ofh_nns     = open(os.path.join(outd, 'nns.tsv'),         'w')
        ofh_qbc     = open(os.path.join(outd, 'qual_by_cyc.tsv'), 'w')
        ofh_unal    = open(os.path.join(outd, 'unal.fastq'),      'a')
        for k in sorted(st.keys()):
            flag, rdlen = k
            if flag == 4:
                ofh_unal_len = open(os.path.join(outd, 'unal_%d.fastq' % rdlen), 'a')
                ofh_unal.write(st[k].readbuf.getvalue())
                ofh_unal_len.write(st[k].readbuf.getvalue())
                st[k].readbuf = StringIO.StringIO()
                ofh_unal_len.close()
            str += "%d:%d:%d:%0.3f%%\t" % \
                (flag, rdlen, st[k].nupdate, 100.0 * float(st[k].nupdate) / nupdate)
            ofh_basics.write("%d\t%d\t%d\n" % (flag, rdlen, st[k].nupdate))
            for m in xrange(0, len(st[k].nns)):
                if st[k].nns[m] > 0:
                    ofh_nns.write("%d\t%d\t%d\t%d\n" % (flag, rdlen, m, st[k].nns[m]))
            for m in xrange(0, len(st[k].avg_qual)):
                if st[k].avg_qual[m] > 0:
                    ofh_means.write("%d\t%d\t%d\t%d\n" % (flag, rdlen, m, st[k].avg_qual[m]))
            for m in xrange(0, len(st[k].median_qual)):
                if st[k].median_qual[m] > 0:
                    ofh_medians.write("%d\t%d\t%d\t%d\n" % (flag, rdlen, m, st[k].median_qual[m]))
            for m in xrange(0, len(st[k].quals_by_cyc)):
                for n in xrange(0, len(st[k].quals_by_cyc[m])):
                    if st[k].quals_by_cyc[m][n] > 0:
                        ofh_qbc.write("%d\t%d\t%d\t%d\t%d\n" % (flag, rdlen, m, n, st[k].quals_by_cyc[m][n]))
        for ofh in (ofh_basics, ofh_medians, ofh_means, ofh_nns, ofh_qbc, ofh_unal):
            ofh.close()
        return str
    
    def run(self):
        with self.glock:
            print "Thread %d started ..." % self.tid
        while True:
            tup = self.queue.get()
            group, name, url, url2 = tup
            with self.glock:
                print "Thread %d handling '%s' ..." % (self.tid, name)
            dir = os.path.join(outdir, group, name)
            n = 0
            
            # Set up SAM output file
            sam_fn = os.path.join(dir, "all.sam")
            if skip_existing and os.path.exists(dir) and os.path.exists(sam_fn):
                print "  skipping because sam file '%s' already exists" % sam_fn
            else:
                mkdir_p(dir)
                sam_fh = open(sam_fn, 'w')
                
                # Set up per-file dict
                file_summ = dict()
                
                cmd = ''
                if download_first:
                    bs = os.path.basename(url)
                    if not os.path.exists(bs):
                        bstmp = bs + ".tmp"
                        if os.path.exists(bstmp): os.remove(bstmp)
                        if os.path.exists(bs): os.remove(bs)
                        cmd = 'curl %s --O %s 2>/dev/null' % (url, bstmp)
                        print "Downloaded %s" % bs
                        fails = 0
                        fails_lim = 5
                        pause = 10
                        while fails < fails_lim:
                            el = os.system(cmd)
                            if el == 0: break
                            print >> sys.stderr, 'Bad return code from "%s": %d' % (cmd, el)
                            fails += 1
                            if os.path.exists(bstmp): os.remove(bstmp)
                            time.sleep(pause)
                        if fails == fails_lim:
                            raise RuntimeError('curl command failed %d times' % fails_lim)
                        os.rename(bstmp, bs)
                    cmd = 'cat ' + bs
                else: cmd = 'curl %s 2>/dev/null' % url
                if url.endswith('.gz'):  cmd += ' | gzip -dc'
                if url.endswith('.bz2'): cmd += ' | bzip2 -dc'
                cmd += ' | %s %s -U - --mm' % (bt2_exe, ' '.join(bt2_args))
                proc = subprocess.Popen(cmd, shell=True, bufsize=128*1024, stdout=subprocess.PIPE)
                
                for ln in proc.stdout:
                    sam_fh.write(ln)
                    if ln.startswith('@'): continue
                    with self.glock:
                        self.nglob.n += 1
                        self.update(self.summs, ln)
                        if self.nglob.n % 100000 == 0:
                            str = self.flush(self.summs, self.nglob.n, self.outdir)
                            print "%s\tThread%d\t%d tasks left" \
                                % (str[:-1], self.tid, self.queue.qsize())
                    n += 1
                    self.update(file_summ, ln)
                    if n % 100000 == 0:
                        self.flush(file_summ, n, dir)
                
                self.flush(file_summ, n, dir)
                sam_fh.close();
            
            with self.glock:
                str = self.flush(self.summs, self.nglob.n, self.outdir)
                print "Thread %d finished '%s' ..." % (self.tid, name)
            
            self.queue.task_done()

workers = [ Worker(x, queue, glock, summs, tot, os.path.join(outdir, group)) \
           for x in xrange(0, ncpus) ]
for w in workers:
    w.setDaemon(True)
    w.start()

fh = open(manifest, 'r')
while True:
    ln = fh.readline().rstrip()
    if len(ln) == 0: break
    ts = string.split(ln, '\t')
    group, name = None, None
    url1, url2 = None, None
    if len(ts) > 3:
        group, name, url1, url2 = ts[0], ts[1], ts[2], ts[3]  # paired-end
    else:
        group, name, url1       = ts[0], ts[1], ts[2]         # unpaired
    with glock: print "Master thread adding '%s' to queue ..." % name
    queue.put([group, name, url1, url2])
fh.close()

print "Joining queue ..."
queue.join()
print "DONE"
