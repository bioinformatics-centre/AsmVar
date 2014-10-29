"""
===============================================================================
Use pysam model to calculate the AS from bam file
===============================================================================
Author : Shujia Huang
Date   : 2014-03-25 14:29:08
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

def IsSNP(refbase, alleles):

    isSnp = True
    for ale in alleles: 
        if len(ale) > 1 or len(refbase) > 1: isSnp = False

    return isSnp

def main(opt): 

    vcfInfile = opt.vcfInfile
    bamInfile = opt.bamInfile
    faInfile  = opt.refInfile
    sampleID  = opt.sample
    refId     = opt.refChrId
    outPrefix = opt.outPrefix
    mapq      = opt.mapq
    newID     = opt.exc
    fa = ATA.LoadFaSeq(faInfile, refId)

    print >> sys.stderr, '# [INFO] Fa Loading finish ***'
    if bamInfile[-4:] == '.bam':
        samInHandle = pysam.Samfile(bamInfile, 'rb')
    else :
        samInHandle = pysam.Samfile(bamInfile, 'r')

    #samOutHandle = pysam.Samfile(outPrefix + '.bam', 'wb', template=samInHandle)
    vcfOutHandle = open(outPrefix + '.vcf', 'w')

    print >> sys.stderr, '# [INFO] Now Scaning the VCF and doing alternate align ... ...'
    if vcfInfile[-3:] == '.gz':
        if refId == "ALL":
            I = os.popen('gzip -dc %s' % vcfInfile)
        else:
            I = os.popen('/home/siyang/Bin/software_pip/tabix-0.2.6/tabix -h %s %s' % (vcfInfile, refId))
    else :
        I = open(vcfInfile)

    frist = True
    while 1:
        lines = I.readlines(100000)
        if not lines : break;

        for line in lines :
        
            line = line.strip('\n')
            col  = line.split()
            if re.search(r'^##', line): 
                if frist and re.search(r'^##FORMAT=', line): 
                    frist = False
                    vcfOutHandle.write('##FORMAT=<ID=AA,Number=4,Type=Integer,Description="Information of Alternate Align. Format: Ref_perfect,Alt_Perfect,Both_Perfect,Both_Imperfect">\n')
                vcfOutHandle.write('%s\n' % line)
                continue
            elif re.search(r'^#CHROM', line):
                if len(col) < 10 :
                    print >> sys.stderr, '# [ERROR] The input vcf file (%s) does not contain the "Sample" fields!\n' % vcfInfile
                    sys.exit(1)
                sam2col = {sam:i+9 for i, sam in enumerate(col[9:])}
                if sampleID not in sam2col: 
                    print >> sys.stderr, '# [ERROR] The input sample id (%s) is not match in Vcf file(%s)\n' % (sampleID, vcfInfile)

                if len(newID) > 0: 
                    vcfOutHandle.write('%s\t%s\n' % ('\t'.join(col[:9]), newID))
                else: 
                    vcfOutHandle.write('%s\t%s\n' % ('\t'.join(col[:9]), sampleID))
                continue

            if refId != 'ALL' and refId != col[0]: continue
            # This is same in 'if col[4] != '.' and not IsSNP(col[3], [col[4]]):'
            if re.search(r'REFCALL', col[6]) or re.search(r'GAP', col[6]):
                continue

            idx = sam2col[sampleID]
            fi  = col[idx].split(':')

            gt  = fi[0].split('/')
            if '|' in fi[0]: gt = fi[0].split('|')
            gtIdx = 1 # Default is the first ALT Sequence
            if len(gt) == 2 and gt[1] != '.': 
                gtIdx = string.atoi(gt[1])

            col[4] = col[4].split(',')[gtIdx-1] # Match to the identity sample
            isAltAlign = False
            zr,za,zc,zi = 0,0,0,0
            if col[4] != '.' and not IsSNP(col[3], [col[4]]):
                # Not SNP, INTERGAP
                isAltAlign = True
            	#if col[2]  == '.': col[2] = 'V_' + col[0] + '_' + col[1]
                zr,za,zc,zi = ATA.Align(samInHandle, 
                                        #samOutHandle, Don't output bam file
                                        fa, 
                                        col[0], 
                                        string.atoi(col[1]), 
                                        col[2], 
                                        col[3], 
                                        col[4][1:], #col[4][0] is reference
                                        mapq)
            # Ignore the position which is meanless
            if not isAltAlign and col[idx] == './.': continue

            fm = {t:i for i, t in enumerate(col[8].split(':'))} # Get Format
            if col[idx] != './.' and len(fi) != len(fm): 
                raise ValueError('[ERROR] The format of "FORMAT"' + 
                                 'fields is not match sample '    + 
                                 '%r in %r' % (col[idx], fm))
            for type in ['VS', 'VT']:
                if type not in fm:
                    raise ValueError('[ERROR] The format of VCF file is ' + 
                                     'not right which you input, it did ' + 
                                     'not contian %s field in FORMAT')
            format = {}
            for k, i in fm.items(): 
                if k != 'GT' and col[idx] != './.': format[k] = fi[i]

            # Use first sample which is not './.' to set VT and VS if col[idx] == './.'
            # This is the same idea with what we do above for 'gtIdx = 1'
            if col[idx] == './.':
                if isAltAlign:
                    isam = [sam for sam in col[9:] if sam != './.' and not re.search(r'^0/0:', sam)]
                    if len(isam) == 0: 
                    # This may happen if appear duplication position and pick the 
                    # REFCALL instand of Variant call when CombineVar with GATK
                        isam = [sam for sam in col[9:] if sam != './.'][0].split(':')
                    else:
                        isam = isam[0].split(':')
                    format['VT'] = isam[fm['VT']]
                    format['VS'] = isam[fm['VS']]

            if isAltAlign:
                format['AA'] = ','.join(str(a) for a in [zr,za,zc,zi])

            gt = fi[fm['GT']].split('/')
            if '|' in fi[fm['GT']]: gt = fi[fm['GT']].split('|')
            for i, g in enumerate(gt):
                if g != '.' and string.atoi(g) > 1: gt[i] = '1'
            if '|' in fi[fm['GT']]:
                fi[fm['GT']] = '|'.join(gt)
            else:
                fi[fm['GT']] = '/'.join(gt)

            col[8] = 'GT:' + ':'.join(sorted(format.keys()))
            # Still keep the origin genotype 
            col[idx] = fi[fm['GT']] + ':' + ':'.join([format[k] for k in sorted(format.keys())])
            vcfOutHandle.write('%s\t%s\n' % ('\t'.join(col[:9]), col[idx]))

    I.close()

    samInHandle.close()
    vcfOutHandle.close()
    #samOutHandle.close()
    #print >> sys.stderr, '# [INFO] Closing the two Ouput files :\n(1) %s\n(2) %s' % (outPrefix + '.bam', outPrefix + '.vcf') 
    print >> sys.stderr, '# [INFO] Closing the two Ouput files :\n  -- %s' % (outPrefix + '.vcf') 

########################################################################
########################################################################
if __name__ == '__main__' :

    usage = "\nUsage : %prog [option] [-v vcfInfile] > Output"
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option("-v", "--vcf", dest="vcfInfile", metavar="VCF", help="Variants. VCF format.", default=[]   )
    optp.add_option("-b", "--bam", dest="bamInfile", metavar="BAM", help="Bam Alignment file.  ", default=[]   )
    optp.add_option("-c", "--chr", dest="refChrId" , metavar="CHR", help="The chr ID of Re."    , default='ALL')
    optp.add_option("-r", "--ref", dest="refInfile", metavar="REF", help="Reference fa format. ", default=[]   )
    optp.add_option("-s", "--smp", dest="sample"   , metavar="SMP", help="Sample ID."           , default=[]   )
    optp.add_option("-o", "--out", dest="outPrefix", metavar="OUT", help="The prefix of output. [out]"       , default = 'out')
    optp.add_option("-q", "--qul", dest="mapq"     , metavar="QUL", help="Threshold of Mapping Quality. [20]", default = '20' )
    optp.add_option("-e", "--exc", dest="exc"      , metavar="EXC", help="Change Sample ID(-s) to be -e in output", default=[])

    opt, _ = optp.parse_args()
    if len(opt.vcfInfile) == 0: optp.error("Required [-v vcfInfile]\n")
    if len(opt.bamInfile) == 0: optp.error("Required [-b bamInfile]\n")
    if len(opt.refInfile) == 0: optp.error("Required [-r reference] Fa format\n")
    if len(opt.sample   ) == 0: optp.error("Required [-s sample ID]\n")

    opt.mapq = string.atoi(opt.mapq)
    print >> sys.stderr, '#[INFO] Paraeters: python', sys.argv[0], '\n\t-v', opt.vcfInfile, \
          '\n\t-b', opt.bamInfile, '\n\t-r', opt.refInfile, '\n\t-s', opt.sample, '\n\t-o', opt.outPrefix, \
          '\n\t-q', opt.mapq     , '\n\t-c', opt.refChrId

    if len(opt.exc) > 0: 
        print >> sys.stderr, '\t-e', opt.exc, '\n'
    else : 
        print >> sys.stderr, '\n'

    main(opt)
    print >> sys.stderr, '*********************** ALL DONE ***********************\n'














