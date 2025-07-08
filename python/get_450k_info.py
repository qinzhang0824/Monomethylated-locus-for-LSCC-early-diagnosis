import os
import os
import sys
import click
import argparse
import re
import codecs
from argparse import RawTextHelpFormatter
import random

parser = argparse.ArgumentParser(description="TMB",formatter_class=RawTextHelpFormatter)
parser.add_argument('--delete',help="high",required=True)
parser.add_argument('--geosample',help="high",required=True)
parser.add_argument('--clinical',help="high",required=True)
parser.add_argument('--tcgarawbeta',help="high",required=True)
parser.add_argument('--georawbeta',help="high",required=True)
parser.add_argument('--allpickedbeta',help="high",required=True)
parser.add_argument('--allpickedsam',help="high",required=True)
parser.add_argument('--trainpickedbeta',help="high",required=True)
parser.add_argument('--trainpickedsam',help="high",required=True)
parser.add_argument('--valpickedbeta',help="high",required=True)
parser.add_argument('--valpickedsam',help="high",required=True)
parser.add_argument('--trainbloodbeta',help="high",required=True)
parser.add_argument('--dir',help="low",required=True)
#parser.add_argument('--out',help="low",required=True)
argv = parser.parse_args()

#sample_info={}
group_info={}

sample_vs_index={}

tumor=[]
normal=[]

geo_sample_list=[]

train_tumor=[]
train_blood=[]

remove=[]

#validation_tumor=[]
#validation_normal=[]
with open(argv.delete) as f:
    for line in f:
        a=line.rstrip().split('\t')
        if a[0] not in remove:
            removetaticappend(a[0])

with open(argv.geosample) as f:
    head=f.readline()
    for line in f:
        info=line.rstrip().split('\t')
        sam=info[0]
        if sam not in geo_sample_list:
            geo_sample_list.append(sam)
            group_info[sam]='Blood'
with open(argv.clinical) as f:
    head=f.readline()
    for line in f:
        item=line.rstrip().split('\t')
        sam_info=item[6].split("-")
        sample_batch=re.match(r'\d+',sam_info[3]).group()
        sampleid='.'.join(sam_info[0:3])+'.'+sample_batch
        if item[7]=='Primary Tumor':
            if sampleid not in remove:
                tumor.append(sampleid)
                group_info[sampleid]='Tumor'
        elif item[7]=='Solid Tissue Normal':
            if sampleid not in remove:
                normal.append(sampleid)
                group_info[sampleid]='Normal'

print(len(group_info),len(tumor),len(normal),len(geo_sample_list))

random.seed(1234)
train_tumor=random.sample(tumor,352)
random.seed(5678)
train_blood=random.sample(geo_sample_list,437)

#HCC_beta='/kyWGS/user_test/zhangqin/liver_methylation/TCGA_H450k/TCGA.LIHC.methylation.HM450.xls'

filter=['FileID','dir_name','filename','TCGA_Case_ID','TCGA_Sample_Type','Composite Element REF']

tcga_probe=[]
geo_probe=[]
all_probe=[]


allpicked_sample=[]
allpicked_index=[]
trainpicked_sample=[]
trainpicked_index=[]
train_sample_vs_index={}
valpicked_sample=[]
valpicked_index=[]
val_sample_vs_index={}

allpicked_sample_blood=[]
allpicked_index_blood=[]
trainpicked_sample_blood=[]
trainpicked_index_blood=[]
train_sample_vs_index_blood={}
valpicked_sample_blood=[]
valpicked_index_blood=[]
val_sample_vs_index_blood={}

'''
allpicked_beta={}
trainpicked_beta={}
valpicked_beta={}
'''

with open(argv.tcgarawbeta) as t,open(argv.georawbeta) as g,open(os.path.join(argv.dir,argv.trainpickedbeta),'w') as tb,open(os.path.join(argv.dir,argv.trainbloodbeta),'w') as bb,open(os.path.join(argv.dir,argv.trainpickedsam),'w') as tsam:
    #asam.write("Sample_Name"+"\t"+"Cancer_Type_Detailed"+"\t"+"Primary_Tumor_Site"+"\n")
    tsam.write("Sample_Name"+"\t"+"Cancer_Type_Detailed"+"\t"+"Primary_Tumor_Site"+"\n")
    #vsam.write("Sample_Name"+"\t"+"Cancer_Type_Detailed"+"\t"+"Primary_Tumor_Site"+"\n")
    for line in t:
        item=line.rstrip().split('\t')
        if item[0] in filter:
            continue
        elif item[0]=='sample':
            tumor_n=0
            #normal_n=0
            train_tumor_n=0
            #train_blood_n=0
            for s in item:
                m=s.replace('-','.')
                '''
                if m in tumor:
                    allpicked_sample.append(m)
                    allpicked_index.append(item.index(s))
                    sample_vs_index[item.index(s)]=m
                    asam.write(m+"\t"+"HNSCC"+"\t"+group_info[m]+"\n")
                    if m in tumor:
                        tumor_n+=1
                    #elif m in normal:
                        #normal_n+=1
                    else:
                        print(s)
                '''
                if m in train_tumor:
                    trainpicked_sample.append(m)
                    trainpicked_index.append(item.index(s))
                    train_sample_vs_index[item.index(s)]=m
                    tsam.write(m+"\t"+"HNSCC"+"\t"+group_info[m]+"\n")
                    train_tumor_n+=1

                '''
                else:
                    if m =="sample":
                        continue
                    else:
                        if m in tumor:
                            valpicked_sample.append(m)
                            valpicked_index.append(item.index(s))
                            val_sample_vs_index[item.index(s)]=m
                            vsam.write(m+"\t"+"HNSCC"+"\t"+group_info[m]+"\n")
                '''
            tb.write('{0}\t{1}\n'.format(item[0],'\t'.join(trainpicked_sample)))
            print(train_tumor_n)
       # elif 'NA' in item:
       #     continue
        else:
            tumor_NA=0
            picked_beta=[]
            #normal_NA=0
            #if item[0] not in allpicked_beta:
            #    allpicked_beta[item[0]]=[]
            #if item[0] not in trainpicked_beta:
            #    trainpicked_beta[item[0]]=[]
            #if item[0] not in valpicked_beta:
            #    valpicked_beta[item[0]]=[]
            #for i in allpicked_index:
            #    if item[i]=='NA' and sample_vs_index[i] in tumor:
            #        tumor_NA+=1
                #if item[i]=='NA' and sample_vs_index[i] in normal:
                    #normal_NA+=1
            #    allpicked_beta[item[0]].append(item[i])
            for i in trainpicked_index:
                if item[i]=='NA' and train_sample_vs_index[i] in train_tumor:
                    tumor_NA+=1
                picked_beta.append(item[i])
                #trainpicked_beta[item[0]].append(item[i])
            #for i in valpicked_index:
                #valpicked_beta[item[0]].append(item[i])

            if tumor_NA>=0.7*train_tumor_n:
            #if tumor_NA>=1 or normal_NA>=1:
                continue
            else:
                if item[0] not in all_probe:
                    all_probe.append(item[0])
                if item[0] not in tcga_probe:
                    tcga_probe.append(item[0])    
                tb.write('{0}\t{1}\n'.format(item[0],'\t'.join(picked_beta)))

    for line in g:
        item=line.rstrip().split('\t')
        if item[0] in filter:
            continue
        elif item[0]=='ID_REF':
            blood_n=0
            #normal_n=0
            train_blood_n=0
            #train_blood_n=0
            for m in item:
                #m=s.replace('-','.')
                '''
                if m in geo_sample_list:
                    allpicked_sample_blood.append(m)
                    allpicked_index_blood.append(item.index(s))
                    sample_vs_index_blood[item.index(s)]=m
                    asam.write(m+"\t"+"HNSCC"+"\t"+group_info[m]+"\n")
                    if m in geo_sample_list:
                        blood_n+=1
                    #elif m in normal:
                        #normal_n+=1
                    else:
                        print(s)
                '''
                if m in train_blood:
                    trainpicked_sample_blood.append(m)
                    trainpicked_index_blood.append(item.index(m))
                    train_sample_vs_index_blood[item.index(m)]=m
                    tsam.write(m+"\t"+"HNSCC"+"\t"+group_info[m]+"\n")
                    train_blood_n+=1 
                ''' 
                else:
                    valpicked_sample_blood.append(m)
                    valpicked_index_blood.append(item.index(s))
                    val_sample_vs_index_blood[item.index(s)]=m
                    vsam.write(m+"\t"+"HNSCC"+"\t"+group_info[m]+"\n")
                '''
            bb.write('{0}\t{1}\n'.format(item[0],'\t'.join(trainpicked_sample_blood)))
            #fo1.write('{0}\t{1}\t'.format(item[0],'\t'.join(allpicked_sample)))
            print(train_blood_n)
       # elif 'NA' in item:
       #     continue
        else:
            blood_NA=0
            blood_picked_beta=[]
            for i in trainpicked_index:
                if item[i]=='NA' and train_sample_vs_index_blood[i] in train_blood:
                    blood_NA+=1
                blood_picked_beta.append(item[i])
               
            if blood_NA>=0.7*train_blood_n:
                continue
            else:
                bb.write('{0}\t{1}\n'.format(item[0],'\t'.join(blood_picked_beta)))      
                if item[0] not in all_probe:
                    all_probe.append(item[0])
                if item[0] not in geo_probe:
                    geo_probe.append(item[0])
   
    
probe_out=open(os.path.join(argv.dir,'tcga_geo_overlap_probe.xls'),'w')
for each_probe in all_probe:
    if each_probe in tcga_probe and each_probe in geo_probe:
        probe_out.write(each_probe+'\n')

tsam.close()
tb.close()
bb.close()

