## awk -v strType=2 -f tagXSstrandedData.awk Aligned.out.sam > Aligned.XS.sam
##
## samtools view -h Aligned.out.bam | awk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > Aligned.XS.bam
##
## strType defines the strandedness of the RNA-seq library prep, and specifies the read (1 or 2) which strand agrees with the RNA strand.
## For the Illumina Tru-Seq dUTP protocol strType=2.
## If the reads are unpaired, they are assumed to be “read 1”.
BEGIN {
    OFS="\t";
    strSym[0]="+";
    strSym[1]="-";
}

{

    printf $0;
   
    if (substr($1,1,1)=="@")
    {
        printf "\n";
        next;
    };

    str=and($2,0x10)/0x10;
   
    if (and($2,0x1)==0)
    {
        mate=1;
    } else
    {
        mate=and($2,0x40)/0x40+2*and($2,0x80)/0x80;
    };

    if (mate>0 && mate <3)
    {
       if (mate!=strType) str=1-str; #revert strand if the mate is opposite
       printf "\t" "XS:A:" strSym[str];
    };    
    
    printf "\n";

}
