'''
Created on May 24, 2016

@author: mluksza
'''
import os
from Log import Log

class PhyloWGS(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        

    def MAF2VCF(self,maffile,outdir):
        
        def getCounts(col2ind,tab):
            tcovcol="Tcov"
            if not tcovcol in col2ind:
                tcovcol="t_depth"
            ncovcol="Ncov"
            if not ncovcol in col2ind:
                ncovcol="n_depth"
            taccol="Tac"
            if not taccol in col2ind:
                taccol="t_alt_count"
            naccol="Nac"
            if not naccol in col2ind:
                naccol="n_alt_count"
                
            Tcov=int(tab[col2ind[tcovcol]])
            Ncov=int(tab[col2ind[ncovcol]])
            if taccol in col2ind:
                Tac=int(tab[col2ind[taccol]])
            else:
                Taf=float(tab[col2ind["Taf"]])
                Tac=int(round(Tcov*Taf))
            if naccol in col2ind:
                Nac=int(tab[col2ind[naccol]])
            else:
                Naf=float(tab[col2ind["Naf"]])
                Nac=int(round(Ncov*Naf))
            return [Tcov,Tac,Ncov,Nac]
        try:
            os.stat(outdir)
        except:
            os.mkdir(outdir)
             
        def getTaf(col2ind,tab):
            try:
                Taf=int(tab[col2ind["Taf"]])    
            except:
                Taf=1
            return Taf
            
        f=open(maffile)
        header=f.readline().strip()
        htab=header.split("\t")
        col2ind={}
        for i in range(0,len(htab)):
            col2ind[htab[i]]=i        
    
        
        poss=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR"]
        vcfheader='#'+"\t".join(poss)
        poss2=poss
        poss2[2]="POS"
        poss2[5]="POS"
#        poss2[7]="POS"
        namemap={}
        for pos in poss:
            namemap[pos]=pos
            
        if "Hugo_Symbol" in col2ind:
            namemap["INFO"]="Hugo_Symbol"
        else:
            poss2[7]="POS"
        if "Reference_Allele" in col2ind:
            namemap["REF"]="Reference_Allele"
        REFcol=namemap["REF"]
        if "Tumor_Seq_Allele2" in col2ind:
            namemap["ALT"]="Tumor_Seq_Allele2"
        ALTcol=namemap["ALT"]
        if "Chromosome" in col2ind:
            namemap["CHROM"]="Chromosome"
        chromcol=namemap["CHROM"]
        if "Start_Position" in col2ind:
            namemap["POS"]="Start_Position"
        poscol=namemap["POS"]
        if "Filter" in col2ind:
            namemap["FILTER"]="Filter"
        namemap["FILTER"]=poscol
        namemap["FORMAT"]=poscol
        namemap["NORMAL"]=poscol
        namemap["TUMOR"]=poscol
        for pos in poss:
            Log.log(str(pos)+" "+str(namemap[pos] in htab),3)
    
        sample2f={}
        mutannot2f={}
        line=f.readline()
        
        while line:
            line=line.strip()
            line=line.replace('"',"")
            tab=line.split("\t")
            [Tcov,Tac,Ncov,Nac]=getCounts(col2ind, tab)
            Taf=getTaf(col2ind, tab)
            line=f.readline()
            sample=tab[col2ind["Sample"]]
            
            if not sample in sample2f:
                sample2f[sample]=outdir+'/'+sample+".vcf"
                sof=open(sample2f[sample],'w')
                sof.write(vcfheader+"\n")
                sof.close()
                mutannot2f[sample]=outdir+'/mutations_'+sample+".txt"
                mof=open(mutannot2f[sample],'w')
                mof.close()
                
            sof=open(sample2f[sample],'a')
            mof=open(mutannot2f[sample],'a')
            nid=str(tab[col2ind[chromcol]])
            nid+="_"+str(tab[col2ind[poscol]])
            nid0=nid
            nid+="_"+str(tab[col2ind[REFcol]])
            nid+="_"+str(tab[col2ind[ALTcol]])
            nid+="_"+str(tab[col2ind["Sample"]])
            nid=nid.replace("chr","")
            mof.write(nid0+"\t"+nid+"\n")
            mof.close()
            ntab=map(lambda pos: tab[col2ind[namemap[pos]]],poss2)
            if not "chr" in ntab[0]:
                ntab[0]="chr"+ntab[0]
            ntab[2]=nid
            ntab[5]="."
            ntab[6]='PASS'
            ntab[8]="DP:AP"
            ntab[9]=str(Ncov)+":"+str(Ncov-Nac)
            ntab[10]=str(Tcov)+":"+str(Tcov-Tac)
            nline="\t".join(ntab)
            sof.write(nline+"\n")
            sof.close()
        f.close()

    def Vinod2VCF(self,maffile,outdir):
        
        def getCounts(col2ind,tab):
            try:
                Tcov=int(tab[col2ind["t_depth"]])
            except:
                Tcov=0
            try:
                Ncov=int(tab[col2ind["n_depth"]])
            except:
                Ncov=0
            try:
                Tac=int(tab[col2ind["t_alt_count"]])
            except:
                Tac=0
            try:
                Nac=int(tab[col2ind["n_alt_count"]])
            except:
                Nac=0
            return [Tcov,Tac,Ncov,Nac]
        try:
            os.stat(outdir)
        except:
            os.mkdir(outdir)
             
        f=open(maffile)
        header=f.readline().strip()
        htab=header.split("\t")
        col2ind={}
        for i in range(0,len(htab)):
            col2ind[htab[i]]=i        
    
        
        poss=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL","TUMOR"]
        vcfheader='#'+"\t".join(poss)
        poss2=poss
        poss2[2]="POS"
        poss2[5]="POS"
        #poss2[7]="POS"
        namemap={}
        for pos in poss:
            namemap[pos]=pos
        if "Hugo_Symbol" in col2ind:
            namemap["INFO"]="Hugo_Symbol"
        else:
            poss2[7]="POS"
        if "Reference_Allele" in col2ind:
            namemap["REF"]="Reference_Allele"
        REFcol=namemap["REF"]
        if "Tumor_Seq_Allele2" in col2ind:
            namemap["ALT"]="Tumor_Seq_Allele2"
        ALTcol=namemap["ALT"]
        if "Chromosome" in col2ind:
            namemap["CHROM"]="Chromosome"
        chromcol=namemap["CHROM"]
        if "Start_Position" in col2ind:
            namemap["POS"]="Start_Position"
        poscol=namemap["POS"]
        if "Filter" in col2ind:
            namemap["FILTER"]="Filter"
        namemap["FILTER"]=poscol
        namemap["FORMAT"]=poscol
        namemap["NORMAL"]=poscol
        namemap["TUMOR"]=poscol
    
        sample2f={}
        mutannot2f={}
        line=f.readline()
        
        while line:
            line=line.strip()
            line=line.replace('"',"")
            tab=line.split("\t")
            
            alt=str(tab[col2ind[ALTcol]])
            ref=str(tab[col2ind[REFcol]])
            if len(alt)!=1 or len(ref)!=1 or alt=="-" or ref=="-":
                line=f.readline()
                continue
            
            [Tcov,Tac,Ncov,Nac]=getCounts(col2ind, tab)
            line=f.readline()
            samplecol="Tumor_Sample_Barcode"
            sample=tab[col2ind[samplecol]]
            sample=str(sample).replace("_", "")
            if not sample in sample2f:
                sample2f[sample]=open(outdir+'/'+sample+".vcf",'w')
                sample2f[sample].write(vcfheader+"\n")
                sample2f[sample].close()
                mutannot2f[sample]=open(outdir+'/mutations_'+sample+".txt",'w')
                mutannot2f[sample].close()
                
            sample2f[sample]=open(outdir+'/'+sample+".vcf",'a')
            mutannot2f[sample]=open(outdir+'/mutations_'+sample+".txt",'a')
            nid=str(tab[col2ind[chromcol]])
            nid+="_"+str(tab[col2ind[poscol]])
            nid0=nid
            nid+="_"+str(tab[col2ind[REFcol]])
            nid+="_"+str(tab[col2ind[ALTcol]])
            nid+="_"+sample
            mutannot2f[sample].write(nid0+"\t"+nid+"\n")
            mutannot2f[sample].close()
            ntab=map(lambda pos: tab[col2ind[namemap[pos]]],poss2)
            ntab[0]="chr"+ntab[0]
            ntab[2]=nid
            ntab[5]="."
            ntab[6]='PASS'
            ntab[8]="DP:AP"
            ntab[9]=str(Ncov)+":"+str(Ncov-Nac)
            ntab[10]=str(Tcov)+":"+str(Tcov-Tac)
            nline="\t".join(ntab)
            sample2f[sample].write(nline+"\n")
            sample2f[sample].close()
        f.close()

    
    def correct_ssm_files(self,mdir,phylowgsdir,sample):
        mf=open(mdir+"/"+sample+".vcf",'r')
        header=mf.readline()
        header=header[1:].strip().split()
        col2pos={}
        pos=0
        for col in header:
            col2pos[col]=pos
            pos+=1
        line=mf.readline()
        mannot={}
        while line:
            line=line.strip().split()
            chrom=line[col2pos["CHROM"]]
            chrom=chrom.replace("chr","").replace("ch","")
            pos=line[col2pos["POS"]]
            ref=line[col2pos["REF"]]
            alt=line[col2pos["ALT"]]
            oldmutid="_".join([chrom,pos])
            newmutid="_".join([chrom,pos,ref,alt,sample])
            mannot[oldmutid]=newmutid
            line=mf.readline()
        mf.close()        
        print mannot
        f=open(phylowgsdir+"/"+sample+"_ssm.txt")
        lines=f.readlines()
        f.close()
        of=open(phylowgsdir+"/"+sample+"_ssm.txt",'w')
        of.write(lines[0])
        for line in lines[1:]:
            tab=line.split("\t")
            tab[1]=mannot[tab[1]]
            nline="\t".join(tab)
            of.write(nline)
        of.close()
        

        
        
