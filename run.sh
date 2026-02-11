# build index for plasmid to remove plasmid read
bwa index plasmid.fa

for i in "offins_AAVS1" "offins_VEGFA" "offins-V-NC" "offins_TRAC" 
do
    (
        NAME=$i
        min_len=90
        min_frac=0.6
        fastp --umi --umi_loc per_read --umi_len 6 --umi_prefix UMI -w 16 -i ${NAME}_raw_1.fq.gz -I ${NAME}_raw_2.fq.gz  -o ${NAME}_R1_filterd.fq.gz -O ${NAME}_R2_filterd.fq.gz
        
        python split_by_R2.py ${NAME}_R1_filterd.fq.gz ${NAME}_R2_filterd.fq.gz $NAME
        
        bwa mem plasmid.fa tmp_${NAME}_AC_R1.fq.gz | samtools sort -n -@16 -o  tmp_${NAME}_AC_plasmid_R1.bam
        python  filter_plasmid_hits.py  --bam tmp_${NAME}_AC_plasmid_R1.bam --r1 tmp_${NAME}_AC_R1.fq.gz --out1 tmp_${NAME}_AC_clean_R1.fq.gz --min_match_len $min_len --min_match_frac $min_frac --drop_policy either 

        bwa mem plasmid.fa tmp_${NAME}_BD_R1.fq.gz | samtools sort -n -@16 -o  tmp_${NAME}_BD_plasmid_R1.bam
        python  filter_plasmid_hits.py  --bam tmp_${NAME}_BD_plasmid_R1.bam --r1 tmp_${NAME}_BD_R1.fq.gz --out1 tmp_${NAME}_BD_clean_R1.fq.gz --min_match_len $min_len --min_match_frac $min_frac --drop_policy either  

        
        bwa mem -t 16 hg38_chr.fa tmp_${NAME}_AC_clean_R1.fq.gz > ${NAME}_AC_clean.sam
        bwa mem -t 16 hg38_chr.fa tmp_${NAME}_BD_clean_R1.fq.gz > ${NAME}_BD_clean.sam
        
        python hotspots_umidedup.py --a ${NAME}_AC_clean.sam --b ${NAME}_BD_clean.sam --out-prefix ${NAME}_result_500_1000_${min_len}_${min_frac} --cluster-distance 500 --pair-distance 1000 --umi-regex "(?<=:UMI_)[A-Za-z0-9_-]+"
        

    ) 
done 
