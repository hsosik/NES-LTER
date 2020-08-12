# Script to excise primer and barcode regions from Illumina MiSeq runs; laready merged....need to remove front and back
#Kristen R. Hunter-Cevera, Marine Biological Laboratory

#Input data, line by line, then rewrite as it comes in?

import sys
import regex

# print sys.argv[1]
# print sys.argv[2]
fastqfile = open(sys.argv[1], 'r')
primer_match = open('%sprimer_match_%s' %(sys.argv[2], sys.argv[3]),'w')
noprimer_match = open('%snoprimer_match_%s' %(sys.argv[2],sys.argv[3]),'w')
# primer_match = open('/home/kristen/Documents/V6V8_analysis/data/MVCO_7_13Jan2020_run/merged_reads/primer_match_%s'%sys.argv[1],'w')
# noprimer_match = open('/home/kristen/Documents/V6V8_analysis/data/MVCO_7_13Jan2020_run/merged_reads/noprimer_match_%s'%sys.argv[1],'w')

#fastqfile = open('/Users/kristenhunter-cevera/MVCO_V6V8/raw_sequence_reads/KHC_1/testfile.R1','r')
# primer_match = open('/Users/kristenhunter-cevera/MVCO_V6V8/raw_sequence_reads/KHC_1/testfile_primermatch.R1','w')
# noprimer_match = open('/Users/kristenhunter-cevera/MVCO_V6V8/raw_sequence_reads/KHC_1/testfile_noprimermatch.R1','w')
# m=regex.search('.R\d{1}',sys.argv[1])
# print m.group()
# if m.group()=='.R1': #forward reads
for line in fastqfile:
    line = line.rstrip()
    if line[0]=='>': #we have a header line!
        name=line
        count=1 #new block
        # print 'found header'
        # print line
    elif count==1: # you should be at the sequence, check for primers!
        seq=line
        # print 'found seq'
        #print seq
        f = regex.search("(.....AAACT[CT]AAA[GT]GAATTGACGG){e<=3}",seq,regex.IGNORECASE) #allow up to two errors- insertions/deletions/subsitutions
        r = regex.search("G([TC]ACACACCGCCCGT){e<=2}",seq,regex.IGNORECASE) #s here is substition, no insertions or deletions allowed
        #print f
        #print r
        #r = regex.search("(ACGGG CGG TGTGT [AG]C){e<=2}",seq)
        count=2
        if f is not None and r is not None:
            fpos=f.end()
            rpos=r.start()
            if fpos < 30 and rpos > len(seq)-20: #a double check that the primers are in the right place for this sequence...
                primer_match.write(str(name) + '\n' + str(seq[fpos:rpos]) + '\n')
            else:
                noprimer_match.write(str(name) + '\n' + str(seq) + '\n')
        # elif f is None and r is not None:
        #     #run the forward primer search again, this time with greater tolerance:
        #     f = regex.search("(.....AAACT[CT]AAA[GT]GAATTGAC){s<=4,i=2,d=2}GG",seq,regex.IGNORECASE) #allow up to two errors- insertions/deletions/subsitutions
        #     if f is not None:
        #         print "found a forward primer"
        #         print f
        #         fpos=f.end()
        #         rpos=r.start()
        #         primer_match.write(str(name) + '\n' + str(seq[fpos:rpos]) + '\n')
        # elif r is None and f is not None:
        #     #run the reverse primer search again, this time with greater tolerance:
        #     r = regex.search("G([TC]ACACACCGCCCGT){e<=4}",seq,regex.IGNORECASE,) #s here is substition, no insertions or deletions allowed
        #     if r is not None:
        #         print "found a reverse primer"
        #         print r
        #         fpos=f.end()
        #         rpos=r.start()
        #         primer_match.write(str(name) + '\n' + str(seq[fpos:rpos]) + '\n')
        else:
            noprimer_match.write(str(name) + '\n' + str(seq) + '\n')

fastqfile.close()
noprimer_match.close()
primer_match.close()
