#!/usr/bin/python

from numpy import *
from itertools import *
import subprocess
import itertools, random
import sys, os
import csv

#print(len(sys.argv))
if len(sys.argv)!=5:
    print "Usage: ./check_duplicates datadir tempdir duplicatesfile datastem \n"
    sys.exit(1)

def calc_concordance(ids1,ids2,ped1,ped2,bim,name,freq):
    #print "Scanning .fam"
    ids1 = [id.split()[1].strip() for id in open(ids1)]
    ids2 = [id.split()[1].strip() for id in open(ids2)]
    freqs = [f.split() for f in open(freq)][1:]
    freqs = {f[1] : f[4] for f in freqs}
    phenos1 = {}
    phenos2 = {}
    #print "Scanning .ped"
    l1 =  [line.split() for line in open(ped1)]
    l2 =  [line.split() for line in open(ped2)]
    #print "Scanning .bim"

    markers = [line.split()[1] for line in open(bim)]

    #print "starting concordance calcluations"
    for l in l1:
        phenos1[l[1]] = map(lambda x : sorted(x)[0]+sorted(x)[1],zip(l[6::2],l[7::2]))
    for l in l2:
        phenos2[l[1]] = map(lambda x : sorted(x)[0]+sorted(x)[1],zip(l[6::2],l[7::2]))

    sum_incorrect_calls    = 0.0
    sum_correct_calls      = 0.0
    sum_one_not_called     = 0.0
    sum_both_not_called    = 0.0
    pairs = 0
    halved_pairs = 0

    in_mismatch = set()

    # track mismatches
    miscount= zeros(len(markers))

    # track mismatches where both are called
    misbothcall= zeros(len(markers))

    pair = [['PAIR','INCORRECT','CORRECT','ONE_NC','BOTH_NC','CONCORDANCE']]
    for id1,id2 in zip(ids1,ids2):
        if (id1 in phenos1) and (id2 in phenos2):
            mismatch   = array(map(lambda x,y : x!=y,phenos1[id1],phenos2[id2]),int32)
            miscount   = miscount + mismatch
            goodmatch  = array(map(lambda x,y : x==y,phenos1[id1],phenos2[id2]),int32)
            both_called = array(map(lambda x,y : not("0" in x or "0" in y) ,phenos1[id1],phenos2[id2]),int32)
            both_not_called = array(map(lambda x,y : ("0" in x and "0" in y) ,phenos1[id1],phenos2[id2]),int32)
            in_mismatch = in_mismatch.union(itertools.compress(markers,mismatch*both_called))
	    misbothcall = misbothcall + (mismatch * both_called) 
            # get concordance stats for each individual duplicate pair
            pair_incorrect_calls = sum(mismatch*both_called)
            pair_correct_calls = sum(goodmatch*both_called)
            pair_one_not_called = sum((1-both_called)-both_not_called)
            pair_both_not_called = sum(both_not_called)
            pair_stats = [id1 +"-"+ id2, pair_incorrect_calls, pair_correct_calls, pair_one_not_called, pair_both_not_called, float(pair_correct_calls)/(pair_correct_calls + pair_incorrect_calls)  if pair_correct_calls>0 else 0]
            pair.append(pair_stats)
            # create concordance stats for sum of markers in group (All markers)
            sum_incorrect_calls   = sum_incorrect_calls + sum(mismatch*both_called)
            sum_correct_calls     = sum_correct_calls   + sum(goodmatch*both_called)
            sum_one_not_called    = sum_one_not_called  + sum((1-both_called)-both_not_called)
            sum_both_not_called   = sum_both_not_called + sum(both_not_called)
            pairs = pairs+1
        elif (id1 in phenos1)+(id2 in phenos2)==1:
            halved_pairs = halved_pairs + 1
            #if id1 in phenos1:
            #    print id2
            #if id2 in phenos2:
            #    print id1

    # test print miscount array
    #print(miscount)
    markermis = zip(markers, miscount)
    with open(tmphome+"/"+name+"-mismatch-per-marker-all", "wb") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(markermis)

    # test print misbothcall array
    bothmis = zip(markers, misbothcall)
    with open(tmphome+"/"+name+"-misbothcall-per-marker-all", "wb") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(bothmis)

    results = [name, "All", pairs, halved_pairs, len(markers),sum_correct_calls/pairs,sum_incorrect_calls/pairs,sum_one_not_called/pairs,sum_both_not_called/pairs,float(sum_correct_calls)/(sum_correct_calls + sum_incorrect_calls)  if sum_correct_calls>0 else 0]
    print "\t".join(map(str,results))
    with open(tmphome+"/"+name+"-total-concordance","a") as f:
        f.write("\t".join(map(str,results)))
        f.write("\n")
        f.close()

    with open(tmphome+"/"+name+"-pairwise-concordance-all-markers", "wb") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(pair)

    # RARE MARKERS
    pair = [['PAIR','INCORRECT','CORRECT','ONE_NC','BOTH_NC','CONCORDANCE']]
    raremarkers = array(map(lambda x: float(x[1]) < 0.05, freqs.iteritems()), int32)

    sum_incorrect_calls    = 0.0
    sum_correct_calls      = 0.0
    sum_one_not_called     = 0.0
    sum_both_not_called    = 0.0
    halved_pairs = 0
    pairs = 0

    for id1,id2 in zip(ids1,ids2):
        if (id1 in phenos1) and (id2 in phenos2):
            mismatch   = array(map(lambda x,y : x!=y,phenos1[id1],phenos2[id2]),int32)
            goodmatch  = array(map(lambda x,y : x==y,phenos1[id1],phenos2[id2]),int32)
            both_called = array(map(lambda x,y : not("0" in x or "0" in y) ,phenos1[id1],phenos2[id2]),int32)
            both_not_called = array(map(lambda x,y : ("0" in x and "0" in y) ,phenos1[id1],phenos2[id2]),int32)
            # get concordance stats for each individual duplicate pair
            #print len(raremarkers)
            #print len(mismatch)
            #print len(both_called)
            pair_incorrect_calls  = sum(mismatch*both_called*raremarkers)
            pair_correct_calls    = sum(goodmatch*both_called*raremarkers)
            pair_one_not_called   = sum(((1-both_called)-both_not_called)*raremarkers)
            pair_both_not_called  = sum(both_not_called*raremarkers)
            pair_stats = [id1 +"-"+ id2, pair_incorrect_calls, pair_correct_calls, pair_one_not_called, pair_both_not_called, float(pair_correct_calls)/(pair_correct_calls + pair_incorrect_calls)  if pair_correct_calls>0 else 0]
            pair.append(pair_stats)
            # create concordance stats for sum of markers in group (Rare markers)
            sum_incorrect_calls   = sum_incorrect_calls + sum(mismatch*both_called*raremarkers)
            sum_correct_calls     = sum_correct_calls   + sum(goodmatch*both_called*raremarkers)
            sum_one_not_called    = sum_one_not_called  + sum(((1-both_called)-both_not_called)*raremarkers)
            sum_both_not_called   = sum_both_not_called + sum(both_not_called*raremarkers)
            pairs = pairs+1
        elif (id1 in phenos1)+(id2 in phenos2)==1:
            halved_pairs = halved_pairs + 1
    #
    #
    results = [name,"~Rare", pairs, halved_pairs,sum(raremarkers),sum_correct_calls/pairs,sum_incorrect_calls/pairs,sum_one_not_called/pairs,sum_both_not_called/pairs,float(sum_correct_calls)/(sum_correct_calls + sum_incorrect_calls) if sum_correct_calls>0 else 0]
    print "\t".join(map(str,results))

    with open(tmphome+"/"+name+"-total-concordance","a") as f:
        f.write("\t".join(map(str,results)))
        f.write("\n")
        f.close()

    with open(tmphome+"/"+name+"-pairwise-concordance-rare-markers", "wb") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(pair)

    # common markers
    pair = [['PAIR','INCORRECT','CORRECT','ONE_NC','BOTH_NC','CONCORDANCE']]
    sum_incorrect_calls    = 0.0
    sum_correct_calls      = 0.0
    sum_one_not_called     = 0.0
    sum_both_not_called    = 0.0
    halved_pairs = 0
    pairs = 0

    for id1,id2 in zip(ids1,ids2):
        if (id1 in phenos1) and (id2 in phenos2):
            mismatch   = array(map(lambda x,y : x!=y,phenos1[id1],phenos2[id2]),int32)
            goodmatch  = array(map(lambda x,y : x==y,phenos1[id1],phenos2[id2]),int32)
            both_called = array(map(lambda x,y : not("0" in x or "0" in y) ,phenos1[id1],phenos2[id2]),int32)
            both_not_called = array(map(lambda x,y : ("0" in x and "0" in y) ,phenos1[id1],phenos2[id2]),int32)
            # get concordance stats for each individual duplicate pair
            pair_incorrect_calls = sum(mismatch*both_called*(1 - raremarkers))
            pair_correct_calls = sum(goodmatch*both_called*(1 - raremarkers))
            pair_one_not_called = sum((1-both_called)-both_not_called*(1 - raremarkers))
            pair_both_not_called = sum(both_not_called*(1 - raremarkers))
            pair_stats = [id1 +"-"+ id2, pair_incorrect_calls, pair_correct_calls, pair_one_not_called, pair_both_not_called, float(pair_correct_calls)/(pair_correct_calls + pair_incorrect_calls)  if pair_correct_calls>0 else 0]
            pair.append(pair_stats)
            # create concordance stats for sum of markers in group (Common markers)
            sum_incorrect_calls   = sum_incorrect_calls + sum(mismatch*both_called*(1 - raremarkers))
            sum_correct_calls     = sum_correct_calls   + sum(goodmatch*both_called*(1 - raremarkers))
            sum_one_not_called    = sum_one_not_called  + sum(((1-both_called)-both_not_called)*(1 - raremarkers))
            sum_both_not_called   = sum_both_not_called + sum(both_not_called*(1 - raremarkers))
            pairs = pairs+1
        elif (id1 in phenos1)+(id2 in phenos2)==1:
            halved_pairs = halved_pairs + 1


    results = [name,"~Common", pairs,halved_pairs, sum(1 -  raremarkers),sum_correct_calls/pairs,sum_incorrect_calls/pairs,sum_one_not_called/pairs,sum_both_not_called/pairs, float(sum_correct_calls)/(sum_correct_calls + sum_incorrect_calls) if sum_correct_calls>0 else 0]
    print "\t".join(map(str,results))

    with open(tmphome+"/"+name+"-total-concordance","a") as f:
        f.write("\t".join(map(str,results)))
        f.write("\n")
        f.close()

    with open(tmphome+"/"+name+"-pairwise-concordance-common-markers", "wb") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(pair)
    # """
    # # testing code I used when writing out names of markers
    # mark3 = []
    # mark2 = []
    # mark1 = []
    # for i in range(len(totals)):
    #     if totals[i]==3:
    #         mark3.append(markers[i])
    #     if totals[i]==2:
    #         mark2.append(markers[i])
    #     if totals[i]==1:
    #         mark1.append(markers[i])
    #
    # with open("/home/reidarst/bad_markers3.txt","w") as f:
    #     f.write("\n".join(mark3))
    # with open("/home/reidarst/bad_markers2.txt","w") as f:
    #     f.write("\n".join(random.sample(mark2,30)))
    # with open("/home/reidarst/bad_markers1.txt","w") as f:
    #     f.write("\n".join(random.sample(mark1,30)))
    # """

#print sys.argv
datahome  = sys.argv[1]
#print("datahome:", datahome)
tmphome   = sys.argv[2]
if not os.path.exists(tmphome):
    os.makedirs(tmphome)
#print("tmphome:", tmphome)
duplicates = sys.argv[3]
#print("dup:",	 duplicates)
datastem = sys.argv[4]

ids1 = []
ids2 = []
id_to_fam = {}

for line in open(os.path.join(datahome,datastem + ".fam")):
    fields = line.split()
    id_to_fam[fields[1]] = fields[0]

for line in open(duplicates):
    id1, id2 = line.split()
    if id1 in id_to_fam and id2 in id_to_fam:
    	ids1.append(id_to_fam[id1] + " " + id1)
     	ids2.append(id_to_fam[id2] + " " + id2)

with open(tmphome+"/dup_ids1","w") as f:
    f.write("\n".join(ids1))
with open(tmphome+"/dup_ids2","w") as f:
    f.write("\n".join(ids2))

# urg, this is a bit stupid
ids1 = [i.split()[1] for i in ids1]
ids2 = [i.split()[1] for i in ids2]

def run_script(script,args):
    proc = subprocess.check_output(script+" "+" ".join(args),shell=True,stderr=subprocess.STDOUT)

#data = [f[-3] for f in os.listdir(datahome) if ".bed" in f]
#data = [os.path.splitext(f)[0] for f in os.listdir(datahome) if ".bed" in f]
#print(data)
data=[datastem]

print "Preliminary statitics for " , sys.argv[-1], " and " , sys.argv[1]
print "\t".join(["Data\t","Markers","Indiv","Dup1s","Dups2"])

for d in data:
    inds = set([line.split()[1] for line in open(datahome+"/"+d+".fam")])
    markers    = len(open(datahome+"/"+d+".bim").readlines())
    d1 = len(inds.intersection(ids1))
    d2 = len(inds.intersection(ids2))
    print "\t".join(map(str,[d,markers,len(inds),d1,d2]))

print ""
print "\t".join(["Data\t","Cat.","N","Halved","Markers","Correct","Miscall","One_NC","Both_NC","Conc"])
with open(tmphome+"/"+datastem+"-total-concordance","w") as f:
    f.write("\t".join(["Data\t","Cat.","N","Halved","Markers","Correct","Miscall","One_NC","Both_NC","Conc"]))
    f.write("\n")
    f.close()
for d in data:
    run_script("plink", [" --bfile " + datahome+ "/"+d+" --keep "+ tmphome+"/dup_ids1","--recode","--out "+tmphome+"/set1"])
    run_script("plink", [" --bfile " + datahome+ "/"+d+" --keep "+ tmphome+"/dup_ids2","--recode","--out "+tmphome+"/set2"])
    run_script("plink", [" --bfile " + datahome+ "/"+d+" --freq --out "+tmphome+"/freqs"])
    run_script("plink", [" --bfile " + datahome+ "/"+d+" --freq --nonfounders --out "+tmphome+"/allfreqs"])
    run_script("plink", [" --bfile " + datahome+ "/"+d+" --hardy --nonfounders --out "+tmphome+"/hwe"])
    calc_concordance(tmphome+"/dup_ids1",tmphome+"/dup_ids2",tmphome+"/set1.ped",tmphome+"/set2.ped", datahome + "/"+d+".bim",d,tmphome+"/freqs.frq")
