#include <stdio.h>
#include <stdlib.h>
#include <string.h>
struct attribute{//used to record separated attributes
    char attribute_name[50];
    char attribute_value[50];
};
struct gtfinfo{
    char seqname[15];
    char source[20];
    char feature[12];
    int start;
    int end;
    char score;
    char strand;
    char frame;
    int attribute_number;
    struct attribute attr[30];
    };
struct gene{
    char gene_id[20];
    char gene_name[20];
    char gene_type[25];
    char chromosome[15];
    int gene_start;//position on chromosome
    int gene_end;//position on chromosome
    int gene_transcriptf_start;
    int gene_transcriptf_end;
    int selected_transcript;//use to selected one transcript and merge features
    };
struct transcript{
    int gene_flag;
    char transcript_id[20];
    char transcript_type[50];
    char transcript_strand;//+ or -
    int transcript_start;//position on chromosome
    int transcript_end;//position on chromosome
    int EXON_length;//length of all exons
    int CDS_length;//length of all CDS
    int transcript_exonf_start;
    int transcript_exonf_end;
    int transcript_CDSf_start;
    int transcript_CDSf_end;
    int transcript_UTRf_start;
    int transcript_UTRf_end;
    int transcript_otherf_start;
    int transcript_otherf_end;
    int is_canonical;//boolean 1(canonical) or 0(non-canonical)
    };
struct exon{
    int transcript_flag;
    char exon_id[20];
    int exon_number;
    int exon_start;
    int exon_end;
    int exon_length;
    };
struct CDS{
    int transcript_flag;
    int exon_flag;
    int CDS_start;
    int CDS_end;
    int CDS_length;
    int CDS_frame;
    };
struct UTR{
    int transcript_flag;
    int exon_flag;
    int UTR_start;
    int UTR_end;
    int UTR_length;
    };
struct other_feature{
    int transcript_flag;
    int exon_flag;
    char feature_type[20];//start_codon,stop_codon
    int feature_start;//position on chromosome
    int feature_end;//position on chromosome
    int feature_length;
    };
struct trans_structure{
    char gene_id[20];
    char transcript_id[20];
    char strand;
    int length;
    int start_coord;//first bp of start codon
    int stop_coord;//last bp of stop codon
    };
void fgetstre(char*a,FILE*f);//scan attribute from file
int fgetattr(char*attrn,char*attrv,FILE*f);//get one attribute from file
int scan_row(FILE*gtf,struct gtfinfo*gtfinfo);//scan a row from gtf file
void read_attribute(struct attribute *attr,FILE *f,int attrn);
int read_row(FILE*gtf,struct gtfinfo*gtfinfo);
int length(char*a);//字符串长度
void swapstr(char*a,char*b);//字符串交换
void copystr(char*a,char*b);//字符串复制
void write_gene(struct gene*gene_list,struct gtfinfo*gtf,int tempg,int tempt);
void write_transcript(struct transcript*transcript_list,struct gtfinfo*gtf,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo);
void write_exon(struct exon*exon_list,struct gtfinfo*gtf,int tempt,int tempe);
void write_CDS(struct CDS*CDS_list,struct gtfinfo*gtf,int tempt,int tempe,int tempc);
void write_UTR(struct UTR*UTR_list,struct gtfinfo*gtf,int tempt,int tempe,int tempu);
void write_other(struct other_feature*other_list,struct gtfinfo*gtf,int tempt,int tempe,int tempo);
void write_close(struct gene*gene_list,struct transcript*transcript_list,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo);
int getpos(struct transcript*translist,struct exon*exonlist,struct other_feature*otherlist,struct trans_structure*transstructurelist);
int main(){
    int genenum=0,transnum=0,exonnum=0,CDSnum=0,UTRnum=0,othernum=0;//calculate number of each feature
    FILE*gtf,*out1,*out2,*out3,*out4;//gtf file
    gtf=fopen("/Users/sunkai/Desktop/gencode.v46.chr_patch_hapl_scaff.annotation.gtf","r");//open gtf file for scan
    out1=fopen("/Users/sunkai/Desktop/Gene_transcript information.tsv","w+");//open output file
    //out2=fopen("/Users/sunkai/Downloads/longest_trans_list.v44.longest_PC.tsv","w+");
    out3=fopen("/Users/sunkai/Downloads/canonical_trans_list.tsv","w+");
    out4=fopen("/Users/sunkai/Desktop/trans_list.v46.tsv","w+");
    struct gtfinfo*gtfrow;
    struct gene*genelist;
    struct transcript*translist;
    struct exon*exonlist;
    struct CDS*CDSlist;
    struct UTR*UTRlist;
    struct other_feature*otherlist;
    struct trans_structure*transstructurelist;
    gtfrow=(struct gtfinfo*)malloc(sizeof(struct gtfinfo)*1);//record a row of gtf

    //scan to count number of each feature
    printf("Start scanning the gtf file.\n");
    int sflag=scan_row(gtf,gtfrow);
    while(sflag!=EOF){
        if(strcmp(gtfrow->feature,"gene")==0)
            genenum++;
        else if(strcmp(gtfrow->feature,"transcript")==0)
            transnum++;
        else if(strcmp(gtfrow->feature,"exon")==0)
            exonnum++;
        else if(strcmp(gtfrow->feature,"CDS")==0)
            CDSnum++;    
        else if(strcmp(gtfrow->feature,"UTR")==0)
            UTRnum++;    
        else
            othernum++;
        sflag=scan_row(gtf,gtfrow);
        
    }
    fclose(gtf);
    //allocate memory for each list
    genelist=(struct gene*)malloc(sizeof(struct gene)*genenum);
    translist=(struct transcript*)malloc(sizeof(struct transcript)*transnum);
    exonlist=(struct exon*)malloc(sizeof(struct exon)*exonnum);
    CDSlist=(struct CDS*)malloc(sizeof(struct CDS)*CDSnum);
    UTRlist=(struct UTR*)malloc(sizeof(struct UTR)*UTRnum);
    otherlist=(struct other_feature*)malloc(sizeof(struct other_feature)*othernum);
    printf("Successfully scanned the gtf file!\n");
    printf("gene number:%d\ntranscript number:%d\nexon number:%d\nCDS number:%d\nUTR number:%d\nother number:%d\n",genenum,transnum,exonnum,CDSnum,UTRnum,othernum);

    //read gtf file and record information
    printf("Start reading the gtf file.\n");
    gtf=fopen("/Users/sunkai/Desktop/gencode.v46.chr_patch_hapl_scaff.annotation.gtf","r");//open gtf file
    int rflag=read_row(gtf,gtfrow);
    int gflag=-1,tflag=-1,eflag=-1,cflag=-1,uflag=-1,oflag=-1;
    while(1){
        if(strcmp(gtfrow->feature,"gene")==0){
            if(gflag>-1)
                write_close(genelist,translist,gflag,tflag,eflag,cflag,tflag,oflag);
            gflag++;
            write_gene(genelist,gtfrow,gflag,tflag);
            copystr((genelist+gflag)->chromosome,gtfrow->seqname);

        }
        else if(strcmp(gtfrow->feature,"transcript")==0){
            tflag++;
            write_transcript(translist,gtfrow,gflag,tflag,eflag,cflag,uflag,oflag);
        }
        else if(strcmp(gtfrow->feature,"exon")==0){
            eflag++;
            write_exon(exonlist,gtfrow,tflag,eflag);
        }
        else if(strcmp(gtfrow->feature,"CDS")==0){
            cflag++;
            write_CDS(CDSlist,gtfrow,tflag,eflag,cflag);
        }
        else if(strcmp(gtfrow->feature,"UTR")==0){
            uflag++;
            write_UTR(UTRlist,gtfrow,tflag,eflag,uflag);
        }
        else{
            oflag++;
            write_other(otherlist,gtfrow,tflag,eflag,oflag);
        }
        rflag=read_row(gtf,gtfrow);
        if(rflag==EOF){
            write_close(genelist,translist,gflag,tflag,eflag,cflag,tflag,oflag);
            break;
        }
    }
    
    for(int i=0;i<transnum;i++){
        translist[i].EXON_length=0;
        translist[i].CDS_length=0;
        for(int j=translist[i].transcript_exonf_start;j<=translist[i].transcript_exonf_end;j++){
            translist[i].EXON_length+=exonlist[j].exon_length;
        }
        for(int j=translist[i].transcript_CDSf_start;j<=translist[i].transcript_CDSf_end;j++){
            translist[i].CDS_length+=CDSlist[j].CDS_length;
        }
    }
    printf("successfully read the gtf file!\nassigned features:\n");
    printf("gene:%d\ntranscript:%d\nexon:%d\nCDS:%d\nUTR:%d\nother:%d\n",gflag+1,tflag+1,eflag+1,cflag+1,uflag+1,oflag+1);
    fclose(gtf);

    int pcnum=0;
    for(int i=0;i<genenum;i++){
        if(strcmp((genelist+i)->gene_type,"protein_coding")==0){
            pcnum++;
        }
    }
    transstructurelist=(struct trans_structure*)malloc(sizeof(struct trans_structure)*pcnum);
    int temp=0;
    for(int i=0;i<genenum;i++){
        if((strcmp((genelist+i)->gene_type,"protein_coding")==0)){
            for(int j=genelist[i].gene_transcriptf_start;j<=genelist[i].gene_transcriptf_end;j++){
                if(translist[j].is_canonical==1&&strcmp(translist[j].transcript_type,"protein_coding")==0){
                    copystr((transstructurelist+temp)->gene_id, (genelist+i)->gene_id);
                    copystr((transstructurelist+temp)->transcript_id,translist[j].transcript_id);
                    transstructurelist[temp].strand=translist[j].transcript_strand;
                    transstructurelist[temp].length=translist[j].EXON_length;
                    int a=getpos(translist+j,exonlist,otherlist,transstructurelist+temp);
                    temp++;
                }
            }
        }
    }


    //select longest transcript
    printf("Start selecting transcripts.\n");
    /*for(int i=0;i<genenum;i++){
        int maxl=0,maxf=genelist[i].gene_transcriptf_start;
        for(int j=genelist[i].gene_transcriptf_start;j<=genelist[i].gene_transcriptf_end;j++){
            if(translist[j].EXON_length>maxl){
                maxl=translist[j].EXON_length;
                maxf=j;
            }
        genelist[i].selected_transcript=maxf;
        }
    }
    fprintf(out1,"chromosome\tgene_id\tgene_name\tgene_type\ttranscript_name\ttranscript_type\tEXON_length\tCDS_length\texon_number\tis_canonical\n");
    for(int i=0;i<genenum;i++)
        fprintf(out1,"%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",(genelist+i)->chromosome,(genelist+i)->gene_id,(genelist+i)->gene_name,(genelist+i)->gene_type,(translist+(genelist+i)->selected_transcript)->transcript_id,(translist+(genelist+i)->selected_transcript)->transcript_type,(translist+(genelist+i)->selected_transcript)->EXON_length,(translist+(genelist+i)->selected_transcript)->CDS_length,(translist+(genelist+i)->selected_transcript)->transcript_exonf_end-(translist+(genelist+i)->selected_transcript)->transcript_exonf_start+1,(translist+(genelist+i)->selected_transcript)->is_canonical);
    for(int i=0;i<genenum;i++){
        int maxl=0,maxf=genelist[i].gene_transcriptf_start;
        for(int j=genelist[i].gene_transcriptf_start;j<=genelist[i].gene_transcriptf_end;j++){
            if(translist[j].EXON_length>maxl){
                if(strcmp(translist[j].transcript_type,"protein_coding")==0){
                    maxl=translist[j].EXON_length;
                maxf=j;
                }
            }
        genelist[i].selected_transcript=maxf;
        }
    }
    fprintf(out2,"chromosome\tgene_id\tgene_name\tgene_type\ttranscript_name\ttranscript_type\tEXON_length\tCDS_length\texon_number\tis_canonical\n");
    for(int i=0;i<genenum;i++)
        fprintf(out2,"%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",(genelist+i)->chromosome,(genelist+i)->gene_id,(genelist+i)->gene_name,(genelist+i)->gene_type,(translist+(genelist+i)->selected_transcript)->transcript_id,(translist+(genelist+i)->selected_transcript)->transcript_type,(translist+(genelist+i)->selected_transcript)->EXON_length,(translist+(genelist+i)->selected_transcript)->CDS_length,(translist+(genelist+i)->selected_transcript)->transcript_exonf_end-(translist+(genelist+i)->selected_transcript)->transcript_exonf_start+1,(translist+(genelist+i)->selected_transcript)->is_canonical);
    */
    for(int i=0;i<genenum;i++){
        int maxl=0,maxf=genelist[i].gene_transcriptf_start;
        for(int j=genelist[i].gene_transcriptf_start;j<=genelist[i].gene_transcriptf_end;j++){
            if(translist[j].is_canonical==1){
                maxf=j;
            }
            genelist[i].selected_transcript=maxf;
        }
    }
    fprintf(out3,"chromosome\tgene_id\tgene_name\tgene_type\ttranscript_name\ttranscript_type\tEXON_length\tCDS_length\texon_number\tis_canonical\n");
    for(int i=0;i<genenum;i++)
        fprintf(out3,"%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",(genelist+i)->chromosome,(genelist+i)->gene_id,(genelist+i)->gene_name,(genelist+i)->gene_type,(translist+(genelist+i)->selected_transcript)->transcript_id,(translist+(genelist+i)->selected_transcript)->transcript_type,(translist+(genelist+i)->selected_transcript)->EXON_length,(translist+(genelist+i)->selected_transcript)->CDS_length,(translist+(genelist+i)->selected_transcript)->transcript_exonf_end-(translist+(genelist+i)->selected_transcript)->transcript_exonf_start+1,(translist+(genelist+i)->selected_transcript)->is_canonical);
    printf("Successfully write the output file!\n");
    fprintf(out4,"gene_id\tgene_type\ttranscript_id\ttranscript_type\tEXON_length\tCDS_length\texon_number\tis_canonical\n");
    for(int i=0;i<transnum;i++)
        fprintf(out4,"%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",(genelist+(translist+i)->gene_flag)->gene_id,(genelist+(translist+i)->gene_flag)->gene_type,(translist+i)->transcript_id,(translist+i)->transcript_type,(translist+i)->EXON_length,(translist+i)->CDS_length,(translist+i)->transcript_exonf_end-(translist+i)->transcript_exonf_start+1,(translist+i)->is_canonical);
    fprintf(out1,"gene_id\ttranscript_id\tstrand\tlength\tstart_coord\tstop_coord\n");
    for(int i=0;i<pcnum;i++)
        fprintf(out1,"%s\t%s\t%c\t%d\t%d\t%d\n",(transstructurelist+i)->gene_id,(transstructurelist+i)->transcript_id,(transstructurelist+i)->strand,(transstructurelist+i)->length,(transstructurelist+i)->start_coord,(transstructurelist+i)->stop_coord);
    fclose(out4);
    fclose(out1);
    //fclose(out2);
    fclose(out3);
    return 0;
}
void fgetstre(char*a,FILE*f){
    int cnt=0;
    while(1){
        char c=fgetc(f);
        if(c=='\n'){
            *(a+cnt)='\0';
            break;
        }
        *(a+cnt)=c;
        cnt++;
    }
}
int fgetattr(char*attrn,char*attrv,FILE*f){
    int cnt=0;
    while(1){
        char c=fgetc(f);
        if(c==' '){
            *(attrn+cnt)='\0';
            break;
        }
        *(attrn+cnt)=c;
        cnt++;
    }
    cnt=0;
    while(1){
        char c=fgetc(f);
        if(c==';'){
            char d=fgetc(f);
            if(d=='\n'){
                *(attrv+cnt)='\0';
                return 0;
            }
            *(attrv+cnt)='\0';
            return 1;
        }
        if(c=='"'){
            continue;
        }
        *(attrv+cnt)=c;
        cnt++;
    }
}
int scan_row(FILE*gtf,struct gtfinfo*gtfinfo){
    char *temp;
    temp=(char*)malloc(sizeof(char)*1000);
    int a=fscanf(gtf,"%s\t%s\t%s\t%d\t%d\t%c\t%c\t%c\t",gtfinfo->seqname,gtfinfo->source,gtfinfo->feature,&(gtfinfo->start),&(gtfinfo->end),&(gtfinfo->score),&(gtfinfo->strand),&(gtfinfo->frame));
    if (a==EOF)
        return EOF;
    fgetstre(temp,gtf);
    return a;
}
int read_row(FILE*gtf,struct gtfinfo*gtfinfo){
    int a=0;
    a=fscanf(gtf,"%s\t%s\t%s\t%d\t%d\t%c\t%c\t%c\t",gtfinfo->seqname,gtfinfo->source,gtfinfo->feature,&(gtfinfo->start),&(gtfinfo->end),&(gtfinfo->score),&(gtfinfo->strand),&(gtfinfo->frame));
    if (a==EOF)
        return EOF;
    int count=0;
    int flag=1;
    while(flag==1){
        flag=fgetattr((gtfinfo->attr[count].attribute_name),(gtfinfo->attr[count].attribute_value),gtf);
        count++;
    }
    gtfinfo->attribute_number=count;
    return a;
}
void swapstr(char*a,char*b){
    char c[99];
    copystr(c,a);
    copystr(a,b);
    copystr(b,c);
}
void copystr(char*a,char*b){//b--->a
    int la=length(a),lb=length(b);
    if(la>lb){
        for(int i=0;i<la;i++){
            if(i<lb)
                *(a+i)=*(b+i);
            else
                *(a+i)='\0';
        }
    }
    else{
        for(int i=0;i<lb;i++){
            *(a+i)=*(b+i);
        *(a+lb)='\0';
    }
    }
}
int length(char*a){
    int cnt=0;
    while(a[cnt]!='\0'){
        if(a[cnt]=='\n'){
            a[cnt]='\0';
            break;
        }
        cnt++;
    }
    return cnt;
}
void write_gene(struct gene*gene_list,struct gtfinfo*gtf,int tempg,int tempt){
    copystr((gene_list+tempg)->chromosome,gtf->seqname);
    copystr((gene_list+tempg)->gene_id,gtf->attr[0].attribute_value);
    copystr((gene_list+tempg)->gene_name,gtf->attr[2].attribute_value);
    copystr((gene_list+tempg)->gene_type,gtf->attr[1].attribute_value);
    (gene_list+tempg)->gene_start=gtf->start;
    (gene_list+tempg)->gene_end=gtf->end;
    (gene_list+tempg)->gene_transcriptf_start=tempt+1;
}
void write_transcript(struct transcript*transcript_list,struct gtfinfo*gtf,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo){
    if(tempt>0){
        transcript_list[tempt-1].transcript_exonf_end=tempe;
        transcript_list[tempt-1].transcript_CDSf_end=tempc;
        transcript_list[tempt-1].transcript_UTRf_end=tempu;
        transcript_list[tempt-1].transcript_otherf_end=tempo;
    }
    transcript_list[tempt].gene_flag=tempg;
    copystr((transcript_list+tempt)->transcript_id,gtf->attr[1].attribute_value);
    copystr((transcript_list+tempt)->transcript_type,gtf->attr[4].attribute_value);
    transcript_list[tempt].transcript_strand=gtf->strand;
    transcript_list[tempt].transcript_start=gtf->start;
    transcript_list[tempt].transcript_end=gtf->end;
    transcript_list[tempt].transcript_exonf_start=tempe+1;
    transcript_list[tempt].transcript_CDSf_start=tempc+1;
    transcript_list[tempt].transcript_UTRf_start=tempu+1;
    transcript_list[tempt].transcript_otherf_start=tempo+1;
    transcript_list[tempt].is_canonical=0;
    for(int i=0;i<gtf->attribute_number;i++){
        if(strcmp(gtf->attr[i].attribute_name,"tag")==0){
            if(strcmp(gtf->attr[i].attribute_value,"Ensembl_canonical")==0)
                transcript_list[tempt].is_canonical=1;
        }
    }
}
void write_exon(struct exon*exon_list,struct gtfinfo*gtf,int tempt,int tempe){
    exon_list[tempe].transcript_flag=tempt;
    copystr((exon_list+tempe)->exon_id,gtf->attr[7].attribute_value);
    (exon_list+tempe)->exon_number=atoi(gtf->attr[6].attribute_value);
    (exon_list+tempe)->exon_start=gtf->start;
    (exon_list+tempe)->exon_end=gtf->end;
    (exon_list+tempe)->exon_length=gtf->end-gtf->start+1;
}
void write_CDS(struct CDS*CDS_list,struct gtfinfo*gtf,int tempt,int tempe,int tempc){
    CDS_list[tempc].transcript_flag=tempt;
    CDS_list[tempc].exon_flag=tempe;
    (CDS_list+tempc)->CDS_start=gtf->start;
    (CDS_list+tempc)->CDS_end=gtf->end;
    (CDS_list+tempc)->CDS_length=gtf->end-gtf->start+1;
    (CDS_list+tempc)->CDS_frame=atoi(&gtf->frame);
}
void write_UTR(struct UTR*UTR_list,struct gtfinfo*gtf,int tempt,int tempe,int tempu){
    UTR_list[tempu].transcript_flag=tempt;
    UTR_list[tempu].exon_flag=tempe;
    (UTR_list+tempu)->UTR_start=gtf->start;
    (UTR_list+tempu)->UTR_end=gtf->end;
    (UTR_list+tempu)->UTR_length=gtf->end-gtf->start+1;
}
void write_other(struct other_feature*other_list,struct gtfinfo*gtf,int tempt,int tempe,int tempo){
    other_list[tempo].transcript_flag=tempt;
    other_list[tempo].exon_flag=tempe;
    copystr((other_list+tempo)->feature_type,gtf->feature);
    (other_list+tempo)->feature_start=gtf->start;
    (other_list+tempo)->feature_end=gtf->end;
    (other_list+tempo)->feature_length=gtf->end-gtf->start+1;
}
void write_close(struct gene*gene_list,struct transcript*transcript_list,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo){
    gene_list[tempg].gene_transcriptf_end=tempt;
    transcript_list[tempt].transcript_exonf_end=tempe;
    transcript_list[tempt].transcript_CDSf_end=tempc;
    transcript_list[tempt].transcript_UTRf_end=tempu;
    transcript_list[tempt].transcript_otherf_end=tempo;
}
int getpos(struct transcript*translist,struct exon*exonlist,struct other_feature*otherlist,struct trans_structure*transstructurelist){
    if (transstructurelist->strand=='+'){
        int templength=0;
        int startcodon=0;
        int stopcodon=0;
        int flagt=0,flagp=0;
        for(int i=translist->transcript_otherf_start;i<=translist->transcript_otherf_end;i++){
            if(strcmp(otherlist[i].feature_type,"start_codon")==0&&flagt==0){
                startcodon=i;
                flagt=1;
            }
            if(strcmp(otherlist[i].feature_type,"stop_codon")==0){
                stopcodon=i;
            }
        }
        for(int i=translist->transcript_exonf_start;i<=translist->transcript_exonf_end;i++){
            if(otherlist[startcodon].exon_flag==i){
                transstructurelist->start_coord=templength+otherlist[startcodon].feature_start-exonlist[i].exon_start+1;
            }
            if(otherlist[stopcodon].exon_flag==i){
                transstructurelist->stop_coord=templength+otherlist[stopcodon].feature_end-exonlist[i].exon_start+1;
            }
            templength+=exonlist[i].exon_length;
        }
        if (templength==translist->EXON_length){
            return 0;
        }
        else{
            return 1;
        }
    }
    if (transstructurelist->strand=='-'){
        int templength=0;
        int startcodon=0;
        int stopcodon=0;
        int flagt=0;
        for(int i=translist->transcript_otherf_start;i<=translist->transcript_otherf_end;i++){
            if(strcmp(otherlist[i].feature_type,"start_codon")==0&&flagt==0){
                startcodon=i;
                flagt=1;
            }
            if(strcmp(otherlist[i].feature_type,"stop_codon")==0){
                stopcodon=i;
            }
        }
        for(int i=translist->transcript_exonf_start;i<=translist->transcript_exonf_end;i++){
            if(otherlist[startcodon].exon_flag==i){
                transstructurelist->start_coord=templength+exonlist[i].exon_end-otherlist[startcodon].feature_end+1;
            }
            if(otherlist[stopcodon].exon_flag==i){
                transstructurelist->stop_coord=templength+exonlist[i].exon_end-otherlist[stopcodon].feature_start+1;
            }
            templength+=exonlist[i].exon_length;
        }
        if (templength==translist->EXON_length){
            return 0;
        }
        else{
            printf("%d\n",templength);
            return 1;
        }
    }
}