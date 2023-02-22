import re
import pickle

## Step 1: download https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/
# GRCh38
# Include CHR, POS, REF, ALT, RS, EAF in pickle files, CHR_POS: CHR_POS_REF_ALT_RS_EAF-CHR_POS_REF2_ALT2_RS2_EAF
# GRCh37
# '''
# 1_1470897_G_GA_rs975766887_0.0;1_1470897_G_A,C,T_rs746611425_4.9188391539596654e-05,0.0,0.0008362026561731431
# '21_9453888': '21_9453888_C_G,T_rs373911182_0.183585313174946,0.0'
# '21_9453945': '21_9453945_C_T_._0.0'
# '''

##Step 2： Generate pickle files containing rs number and MAF
fr=open("/work/long_lab/qli/ref/dbSNP/latest_release/VCF/hg38_ref_sequences.pickle","rb")
hg38_ref_sequences = pickle.load(fr)
fr.close()
cur_chr_index=0
cur_chr=hg38_ref_sequences[cur_chr_index]

info_dict={}
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMN10492695
with open("/work/long_lab/qli/ref/dbSNP/population_frequency/latest_release/freq.vcf","r") as fr:
    line=fr.readline()
    while line:
        if not line.startswith("#"):
            line_list=line.strip().split("\t")
            if (line_list[0] in hg38_ref_sequences) and (line_list[0] == cur_chr):
                if(len(line_list)==21):
                    key=str(cur_chr_index+1)+"_"+line_list[1]
                    REFS=line_list[3]
                    ALTS=line_list[4]
                    REFs_EAFs, ALTs_EAFs=line_list[9].split(":") #EAF
                    REFS_list,ALTS_list,REFs_EAFs_list,ALTs_EAFs_list,EAFs_list=[],[],[],[],[]
                    
                    if "," in REFS:
                        REFS_list = REFS.split(",")
                    else:
                        REFS_list.append(REFS)
                        
                    if "," in ALTS:
                        ALTS_list = ALTS.split(",")
                    else:
                        ALTS_list.append(ALTS)
                        
                    if "," in REFs_EAFs:
                        REFs_EAFs_list = REFs_EAFs.split(",")
                    else:
                        REFs_EAFs_list.append(REFs_EAFs)
                    
                    if "," in ALTs_EAFs:
                        ALTs_EAFs_list = ALTs_EAFs.split(",")
                    else:
                        ALTs_EAFs_list.append(ALTs_EAFs)
                    
                    if(len(REFS_list)==len(REFs_EAFs_list)) and (len(ALTS_list)==len(ALTs_EAFs_list)):
                        if not "0" in REFs_EAFs_list : 
                            for REF_EAF in REFs_EAFs_list:
                                if not (int(REF_EAF)==0):
                                    for ALT_EAF in ALTs_EAFs_list:
                                        EAFs_list.append(str(float(ALT_EAF)/float(REF_EAF)))
                            value = str(cur_chr_index+1)+"_"+line_list[1]+"_"+REFS+"_"+ALTS+"_"+line_list[2]+"_"+",".join(EAFs_list)
                    
                            ###Check there are previous info
                            if not key in info_dict:
                                info_dict[key]=value
                            else:
                                pre_value = info_dict.get(key)
                                new_value = pre_value+";"+value
                                info_dict[key]=new_value
            elif (line_list[0] in hg38_ref_sequences):
                ###end of one chr
                fw=open("/work/long_lab/qli/ref/dbSNP/population_frequency/latest_release/hg38."+cur_chr+".pickle","wb")
                pickle.dump(info_dict, fw)
                fw.close() 
                ###initialize dict
                if (cur_chr_index < 21):
                    cur_chr_index=cur_chr_index+1
                    cur_chr=hg38_ref_sequences[cur_chr_index]
                    print(cur_chr)
                    info_dict={}
                    if(len(line_list)==21):
                        key=str(cur_chr_index+1)+"_"+line_list[1]
                        REFS=line_list[3]
                        ALTS=line_list[4]
                        REFs_EAFs, ALTs_EAFs=line_list[9].split(":")
                        REFS_list,ALTS_list,REFs_EAFs_list,ALTs_EAFs_list,EAFs_list=[],[],[],[],[]
                        
                        if "," in REFS:
                            REFS_list = REFS.split(",")
                        else:
                            REFS_list.append(REFS)

                        if "," in ALTS:
                            ALTS_list = ALTS.split(",")
                        else:
                            ALTS_list.append(ALTS)

                        if "," in REFs_EAFs:
                            REFs_EAFs_list = REFs_EAFs.split(",")
                        else:
                            REFs_EAFs_list.append(REFs_EAFs)

                        if "," in ALTs_EAFs:
                            ALTs_EAFs_list = ALTs_EAFs.split(",")
                        else:
                            ALTs_EAFs_list.append(ALTs_EAFs)

                        if(len(REFS_list)==len(REFs_EAFs_list)) and (len(ALTS_list)==len(ALTs_EAFs_list)):
                            if not "0" in REFs_EAFs_list : 
                                for REF_EAF in REFs_EAFs_list :
                                    if not (int(REF_EAF)==0) :
                                        for ALT_EAF in ALTs_EAFs_list:
                                            EAFs_list.append(str(float(ALT_EAF)/float(REF_EAF)))
                                value = str(cur_chr_index+1)+"_"+line_list[1]+"_"+REFS+"_"+ALTS+"_"+line_list[2]+"_"+",".join(EAFs_list)

                                ###Check there are previous info
                                if not key in info_dict:
                                    info_dict[key]=value
                                else:
                                    pre_value = info_dict.get(key)
                                    new_value = pre_value+";"+value
                                    info_dict[key]=new_value
                    else:
                        break
                    
        line=fr.readline()

#save chr_index=21        
fw=open("/work/long_lab/qli/ref/dbSNP/population_frequency/latest_release/hg38."+cur_chr+".pickle","wb")
pickle.dump(info_dict, fw)
fw.close() 

##Step 3: Read your SNP_IDs and convert them to rs
fr=open("/work/long_lab/qli/ref/dbSNP/latest_release/VCF/hg38_ref_sequences.pickle","rb")
hg38_ref_sequences = pickle.load(fr)
fr.close()
cur_chr_index=0
cur_chr=hg38_ref_sequences[cur_chr_index]

fw=open("Example_rs.txt","w")
cur_chr_variant_id_2_rs_dict=pickle.load(open("hg38."+str(cur_chr)+".pickle","rb"))
with open("Example_SNP_ID.txt","r") as f:
    line=f.readline()
    while line
        chr_,pos,ref,alt,build = line.split("_")
        if "chr" in chr_.lower():
            res = re.search(r'\d{1,2}',chr_)
            if len(res>0):
                chr_num = res.group()
        else:
            chr_num =chr_
        SNP_ID=str(chr_num)+"_"+pos+"_"+ref+"_"+alt
        if chr_num != cur_chr_index:
            cur_chr = hg38_ref_sequences[int(chr_num)]
            cur_chr_variant_id_2_rs_dict = pickle.load(open("hg38."+str(cur_chr)+".pickle","rb"))
        if SNP_ID in cur_chr_variant_id_2_rs_dict:
            rs = cur_chr_variant_id_2_rs_dict.get(SNP_ID)
            fw.write(SNP_ID+"\t"+rs+"\n")
        else:
            print(SNP_ID+" does not have rs in current dbSNP")
        line=f.readline()
fw.close()
