import re

###Generate TFEnformer/*.sannot.gz from dbSNP rs files
##This *.py file can conver SNP_ID to rs based on recent dbSNP. Only applicable to point mutations

##Step 1: Download vcf file from dbSNP: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/, unzip it. The unzip vcf file is about 160G. 
##GCF_000001405.25.gz for GRCh37 and GCF_000001405.39.gz for GRCh38
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
#gunzip GCF_000001405.25.gz
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz
#gunzip GCF_000001405.39.gz


##Step 2: Choose right chromosomes dictionary based on your genome build (GRCh37.p13.json or GRCh38.p13.json)
genom_build="GRCh37" #set GRCh37 by default, please change to match your genome build
if genom_build=="GRCh37":
	target_chr_dict = json.load(open("GRCh37.p13.json","r"))
elif genom_build=="GRCh38":
	target_chr_dict = json.load(open("GRCh38.p13.json","r"))


##Step 3: Generate one json file for each chromosome to store SNP_ID:rs information 
##Firsly, let's take a look at the head 2 lines of GCF.*.gz file with "\t" as delimiter 
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
##NC_000001.10    10001   rs1570391677    T       A,C     .       .       RS=1570391677;...;...;
cur_chr=1
cur_chr_variant_id_2_rs_dict={}
point_mutation_set=("A","T","C","G")
vcf_file_br=open("GCF_000001405.25","r")
vcf_line = vcf_file_br.readline()
while vcf_line:
    if not vcf_line.startswith("#"):
        vcf_line_arr = vcf_line.split("\t")
        if vcf_line_arr[0] in target_chr_dict:
            chr_ = target_chr_dict.get(vcf_line_arr[0])
            if chr_ == str(cur_chr):
                if vcf_line_arr[3] in point_mutation_set:
                    ref = vcf_line_arr[3]
                    if "," in vcf_line_arr[4]:
                        alts = vcf_line_arr[4].split(",")
                        for alt in alts:
                            key = chr_+"_"+vcf_line_arr[1]+"_"+ref+"_"+alt
                            cur_chr_variant_id_2_rs_dict[key]=vcf_line_arr[2]
                    else:
                        alt = vcf_line_arr[4]
                        key = chr_+"_"+vcf_line_arr[1]+"_"+ref+"_"+alt
                        cur_chr_variant_id_2_rs_dict[key]=vcf_line_arr[2]
            else:
                json.dump(cur_chr_variant_id_2_rs_dict,open("GCF_000001405.25."+str(cur_chr)+".json","w"))
                cur_chr_variant_id_2_rs_dict={}
                cur_chr+=1
    vcf_line = vcf_file_br.readline()
json.dump(cur_chr_variant_id_2_rs_dict,open("GCF_000001405.25."+str(cur_chr)+".json","w"))
vcf_file_br.close()

##Step 4: Read your SNP_IDs and convert them to rs
fw=open("Example_rs.txt","w")
cur_chr="1"
cur_chr_variant_id_2_rs_dict=json.load(open("GCF_000001405.25."+str(cur_chr)+".json","r"))
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
        if str(chr_num) != cur_chr:
            cur_chr = str(chr_num)
            cur_chr_variant_id_2_rs_dict = json.load(open("GCF_000001405.25."+str(cur_chr)+".json","r"))
        if SNP_ID in cur_chr_variant_id_2_rs_dict:
            rs = cur_chr_variant_id_2_rs_dict.get(SNP_ID)
            fw.write(SNP_ID+"\t"+rs+"\n")
        else:
            print(SNP_ID+" does not have rs in current dbSNP")
        line=f.readline()
fwc.close()
