'''
Created on Apr 13, 2016

@author: mluksza
'''
import sys,os
from PhyloWGS import PhyloWGS

def main(argv):
    step = float(argv[1])
    phylowgs = PhyloWGS()
    if step == 1.0: #convert maf files to vcf
        maffile = argv[2] 
        phylowgsdir = argv[3]
        phylowgs.MAF2VCF(maffile, phylowgsdir)
    if step == 1.5: #convert Vinod's files to vcf
        maffile = argv[2] 
        phylowgsdir = argv[3]
        phylowgs.Vinod2VCF(maffile, phylowgsdir)
    elif step == 2.0: #correct phylowgs parser output
        mdir = sys.argv[2]
        phylowgsdir = sys.argv[3]
        samples0 = os.listdir(mdir)
        print "samples:"
        print samples0
        samples = []
        for fsample in samples0:
            if ".vcf" in fsample:
                samples.append(fsample)
        samples = list(set(map(lambda fsample: fsample.split(".vcf")[0],samples)))
        for sample in samples:
            print sample
            phylowgs.correct_ssm_files(mdir, phylowgsdir, sample)
#    elif step==3.0:
#        mdir=sys.argv[2]
#        ctype=sys.argv[3]
#        suff=""
#        if ctype=="ssm":
#            suff="_"+sys.argv[3]
#        samples0=os.listdir(mdir)
#        samples=[]
#        for ssmpath0 in samples0:
#            sample=ssmpath0
#            print "ssmpath:",ssmpath0
#            ssmpath=ssmpath0+"/"+ssmpath0+"_ssm.txt"
#            if "_ssm.txt" in ssmpath:
#                samples.append(sample)
#                print "sample:",sample
#                jsonpath=mdir+"/"+sample+suff+"/summ_"+sample+".json.gz"
#                if not os.path.exists(jsonpath):
#                    jsonpath=mdir+"/"+sample+suff+"/summ_"+sample+".json"
                                
#                print jsonpath, os.path.exists(jsonpath)
#                phylowgs.processJsonFrequencyFile(jsonpath, sample,mdir+"/"+ssmpath)
        
if __name__ == '__main__':
    main(sys.argv)

