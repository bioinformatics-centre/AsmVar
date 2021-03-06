"""
=================================================
Use pysam model to calculate the AS from bam file
=================================================
Author : Shujia Huang
Date   : 2014-04-28 
"""
import sys
import re
import optparse
import os
import string
import pysam 
import matplotlib.pyplot as plt
import numpy as np
import AlterAlign as ATA


def main ( opt ) : 
    vcfInfile = opt.vcfInfile
    bamInfile = opt.bamInfile
    faInfile  = opt.refInfile
    sampleID  = opt.sample
    outPrefix = opt.outPrefix
    mapq      = opt.mapq
    newID     = opt.exc

    fa     = ATA.LoadFaSeq( faInfile )
    print >> sys.stderr, '# [INFO] Fa Loading finish ***'
    if bamInfile[-4:] == '.bam' :
        samInHandle = pysam.Samfile( bamInfile, 'rb' )
    else :
        samInHandle = pysam.Samfile( bamInfile, 'r'  )

    samOutHandle = pysam.Samfile(outPrefix + '.bam', 'wb', template=samInHandle)
    vcfOutHandle = open(outPrefix + '.vcf', 'w')

    print >> sys.stderr, '# [INFO] Now Scaning the VCF and doing alternate align ... ...'
    if vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcfInfile )
    else :
        I = open ( vcfInfile )

    frist = True
    while 1 :
        lines = I.readlines( 100000 )
        if not lines : break;

        for line in lines :
        
            line = line.strip('\n')
            col  = line.split()
            if re.search(r'^##', line) : 
                if frist and re.search(r'^##FORMAT=', line) : 
                    frist = False
                    vcfOutHandle.write('##FORMAT=<ID=AA,Number=1,Type=String,Description="Information of Alternate Align. Format: Ref_perfect,Alt_Perfect,Both_Perfect,Both_Imperfect">\n')
                vcfOutHandle.write('%s\n' % line)
                continue
            elif re.search(r'^#CHROM', line) :
                if len(col) < 10 :
                    print >> sys.stderr, '# [ERROR] The input vcf file (%s) does not contain the "Sample" fields!\n' % vcfInfile
                    sys.exit(1)
                sam2col = { sam:i+9 for i,sam in enumerate(col[9:]) }
                if sampleID not in sam2col : 
                    print >> sys.stderr, '# [ERROR] The input sample id (%s) is not match in Vcf file(%s)\n' % ( sampleID, vcfInfile)
                if len( newID ) > 0 : vcfOutHandle.write('%s\t%s\n' % ('\t'.join(col[:9]), newID) )
                else : vcfOutHandle.write('%s\t%s\n' % ('\t'.join(col[:9]), sampleID) )
                continue

            if len(col[3]) == len( col[4] ) and len(col[3]) == 1 : continue # SNP
            format = { t:i for i,t in enumerate(col[8].split(':')) } # Get Format
            if 'AA' not in format :
                print >> sys.stderr, '# [ERROR] The format of VCF file is not right which you input, it did not contian %s field' % 'AA'
                sys.exit(1)

            pos    = string.atoi(col[1])
            idx    = sam2col[sampleID]
            if col[idx] == './.' : continue
            fi     = col[idx].split(':')
            gt     = fi[0].split('/')
            if '|' in fi[0] : gt = fi[0].split('|')
            gtIdx  = 1
            if len(gt) == 2 and gt[1] != '.' : gtIdx = string.atoi( gt[1] )

            zr,za  = string.atof( fi[format['AA']].split(',')[0] ), string.atof( fi[format['AA']].split(',')[1] )
            zc,zi  = string.atof( fi[format['AA']].split(',')[2] ), string.atof( fi[format['AA']].split(',')[3] )
            if zc == 0.0 : zc = 1.0
            col[4] = col[4].split(',')[gtIdx-1]
            altseq = col[4][1:]
            probR2A= zc / (zc + zr)
            probA2R= zc / (zc + za)

            if pos > 0 :
                ATA.Align ( samInHandle, samOutHandle, fa, col[0], pos, col[2], col[3], altseq, mapq )
         

    I.close()

    samInHandle.close()
    samOutHandle.close()
    vcfOutHandle.close()
    print >> sys.stderr, '# [INFO] Closeing The two Ouput files :\n(1) %s\n(2) %s' % (outPrefix + '.bam', outPrefix + '.vcf') 

########################################################################
########################################################################
if __name__ == '__main__' :

    usage = "\nUsage : %prog [option] [-v vcfInfile] > Output"
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option("-v", "--vcf", dest="vcfInfile", metavar="VCF", help="Variants. VCF format.", default=[] )
    optp.add_option("-b", "--bam", dest="bamInfile", metavar="BAM", help="Bam Alignment file. " , default=[] )
    optp.add_option("-r", "--ref", dest="refInfile", metavar="REF", help="Reference fa format." , default=[] )
    optp.add_option("-s", "--smp", dest="sample"   , metavar="SMP", help="Sample ID."           , default=[] )
    optp.add_option("-o", "--out", dest="outPrefix", metavar="OUT", help="The prefix of output. [out]"       , default = 'out' )
    optp.add_option("-q", "--qul", dest="mapq"     , metavar="QUL", help="Threshold of Mapping Quality. [20]", default = '20'  )
    optp.add_option("-e", "--exc", dest="exc"      , metavar="EXC", help="Change Sample ID(-s) to be -e in output", default=[] )

    opt, _ = optp.parse_args()
    if len( opt.vcfInfile ) == 0 : optp.error( "Required [-v vcfInfile]\n" )
    if len( opt.bamInfile ) == 0 : optp.error( "Required [-b bamInfile]\n" )
    if len( opt.refInfile ) == 0 : optp.error( "Required [-r reference] Fa format\n" )
    if len( opt.sample    ) == 0 : optp.error( "Required [-s sample ID]\n" )

    opt.mapq = string.atoi( opt.mapq )
    print >> sys.stderr, '# [INFO] Paraeters: python' , sys.argv[0], '\n\t-v', opt.vcfInfile, \
          '\n\t-b', opt.bamInfile, '\n\t-r', opt.refInfile, '\n\t-s', opt.sample, '\n\t-o', opt.outPrefix, \
          '\n\t-q', opt.mapq
    if len( opt.exc ) > 0 : print >> sys.stderr, '\t-e', opt.exc, '\n'
    else                  : print >> sys.stderr, '\n'

    main ( opt )
    print >> sys.stderr, '*********************** ALL DONE ***********************\n'














