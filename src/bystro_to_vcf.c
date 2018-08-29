/*

The code itself is Copyright (C) 2017, by David J. Cutler.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
This library is distributed in the hope that it will be useful,free
but WITHOUT ANY WARRANTY; without even the implied warranty offs
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define TRUE 1
#define FALSE 0
#define stdprn  ((FILE *)NULL)
#include <time.h>

#define minim(a,b) ((a<b)?a:b)
#define maxim(a,b) ((a>b)?a:b)

typedef struct sample_node
{
  int number;
  char *name;
  char *call;
} SAMPLE_NODE;

int *ivector(int nl,int nh);
char *cvector(int nl,int nh);
inline int qual_to_phred(double qual);
unsigned char *ucvector(int nl,int nh);
unsigned long *ulvector(int nl,int nh);
long *lvector(int nl,int nh);
unsigned int *uvector(int nl,int nh);
unsigned long long *ullvector(int nl,int nh);
double *dvector(int nl, int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
char **cmatrix(int nrl,int nrh,int ncl,int nch);
unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch);
unsigned long **ulmatrix(int nrl,int nrh,int ncl,int nch);
unsigned int **umatrix(int nrl,int nrh,int ncl,int nch);
void free_cvector(char *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_ucvector(unsigned char *v, int nl, int nh);
void free_ulvector(unsigned long *v, int nl, int nh);
void free_uvector(unsigned int *v, int nl, int nh);
void free_dvector(double *v, int nl,int nh);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
void free_cmatrix(char **m,int nrl,int nrh,int ncl,int nch);
void free_ucmatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch);
void free_ulmatrix(unsigned long **m,int nrl,int nrh,int ncl,int nch);
void free_umatrix(unsigned int **m,int nrl,int nrh,int ncl,int nch);
void dump_error(char *error_text);
int **imatrix(int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
int find_chrom(unsigned int *pos,int first, int last,int try,unsigned this);
void make_allele_count(SAMPLE_NODE **sn,int *count, int tot_samples);
SAMPLE_NODE *sample_alloc(char *s,int n);
int sort_sample(const void *a,const void *b);
int search_sample(const void *a,const void *b);
char *process_samples(SAMPLE_NODE **sorted_samples,char *fields,char *geno,int tot_samples);


static FILE *outfile;
static long genome_size;
static int MAX_FILE_BUFFER = 512000000;

int main(int argc,char *argv[])
{
	char sss[4196];
	char sdxname[1024];
	char **contig_names,*token;
	unsigned int *contig_starts;
	int max_contigs,no_contigs;
	char *genome_buffer;
	int i,j;

	FILE *sfile,*samplefile;
	gzFile reffile,snpfile;


	if( (argc != 4) )
	{
		printf("\nUsage: %s sdx_file seqant_file sample_list_file \n",argv[0]);
		exit(1);
	}

	strcpy(sdxname,argv[1]);
	if((sfile=fopen(sdxname,"r"))==(FILE *)NULL)
	{
	    printf("\n Can not open file %s\n",sdxname);
	    exit(1);
	}
	if(strstr(sdxname,".sdx") != NULL)
	{
	    for(i=strlen(sdxname)-1;i>0;i--)
	      if(sdxname[i] == '.')
	      {
		  sdxname[i] = '\0';
		  i = 0;
	      }
	}


	fgets(sss,256,sfile);
	max_contigs = atoi(sss);
	no_contigs = max_contigs;
	contig_starts = uvector(0,max_contigs);
	contig_names = cmatrix(0,max_contigs,0,256);

	contig_starts[0] = 0;
	for(i=0;i<max_contigs;i++)
	{
	    fgets(sss,1024,sfile);
	    token = strtok(sss,"\t \n");
	    contig_starts[i+1] = atoi(token);
	    token = strtok(NULL,"\t \n");
	    strcpy(contig_names[i],token);
	    // printf("\nFor contig %d we have offset %d",i,contig_starts[i]);
	}
	fclose(sfile);
	for(i=1;i<=max_contigs;i++)
	  contig_starts[i] += contig_starts[i-1];

	// for(i=0;i<=max_contigs;i++)
	// printf("\n Contig %d starts at position %u \n",i,contig_starts[i]);

	genome_size = (long)contig_starts[max_contigs]+(long)15*max_contigs;
	genome_buffer = (char *)malloc(sizeof(char)*(genome_size+1));
	if(!genome_buffer)
	  dump_error("\n Failed to allocate memory for the Genome Buffer \n");
        // else
	// printf("\n Genome size is %ld \n\n",genome_size);
	sprintf(sss,"%s.seq",sdxname);
	if((reffile=gzopen(sss,"r"))==(gzFile)NULL)
	{
	    printf("\n Can not open file %s for reading\n",sss);
	    exit(1);
	}
	gzbuffer(reffile,33554432);
	long g_temp = 0;
	// printf("\n About to read genome \n\n");
	if(genome_size < MAX_FILE_BUFFER)
		g_temp = gzread(reffile,(void *)genome_buffer,(unsigned int)sizeof(char)*genome_size);
	else
	{
		long count = 0;
		while(count < genome_size)
		{
		    int ttemp = minim((long)genome_size - count,MAX_FILE_BUFFER);
		    g_temp += gzread(reffile,(void *)&genome_buffer[count],ttemp);
		    count += ttemp;
		}
	}
	gzclose(reffile);

        for(i=1,j=15;i<=no_contigs;i++,j+=15)
                contig_starts[i] += j;

	printf("##fileformat=VCFv4.0\n");
	time_t ttt = time(NULL);
	struct tm tm = *localtime(&ttt);
	printf("##fileDate=%d%d%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
	printf("##reference=%s\n",argv[1]);
	printf("##phasing=none\n");
	printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	printf("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n");
	printf("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Called\">\n");
	printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	if((snpfile=gzopen(argv[2],"r"))==(gzFile)NULL)
	{
	    printf("\n Can not open file %s for reading\n",argv[2]);
	    exit(1);
	}
	char *buffer;
	int buff_len = 41959999;
	buffer = cvector(0,buff_len+1);
	int tot_samples = 0;
	if((samplefile=fopen(argv[3],"r"))==(FILE *)NULL)
	{
	    printf("\n Can not open file %s which should be the sample file\n",argv[3]);
	    exit(1);
	}
	fgets(buffer,4196,samplefile);
	SAMPLE_NODE **all_samples,**sorted_samples;
	int max_samples = 65000;

	if( (all_samples = malloc(sizeof(SAMPLE_NODE *)*max_samples)) == NULL)
	{
	    printf("\n Can't allocate memory for %d max samples. Try decreasing it \n",max_samples);
	    exit(1);
	}

	while(strlen(buffer)>1)
	{
	    token = strtok(buffer,"\t \n");
	    printf("\t%s",token);
	    all_samples[tot_samples] = sample_alloc(token,tot_samples);
	    tot_samples++;
	    if(tot_samples >= max_samples)
	    {
		printf("\n We just encountered more samples than we can handle.  Max = %d \n",max_samples);
		exit(1);
	    }
	    buffer[0] = '\0';
	    if(!feof(samplefile))
	      fgets(buffer,4196,samplefile);
	}

	if( (sorted_samples = malloc(sizeof(SAMPLE_NODE *)*tot_samples)) == NULL)
	{
	    printf("\n Can't allocate memory for %d sorted samples. Try decreasing it \n",tot_samples);
	    exit(1);
	}
	for(i=0;i<tot_samples;i++)
	  sorted_samples[i] = all_samples[i];
	// printf("\n About to sort samples \n\n");
	qsort(sorted_samples,tot_samples,sizeof(SAMPLE_NODE *),sort_sample);
	int ref_field = -1;
	int alt_field = -1;
	int missing_field = -1;
	int rs_field_gnomad_genomes = -1;
	int rs_field_gnomad_exomes = -1;
	int rs_field_dbSNP = -1;
	int hom_field = -1;
	int type_field = -1;
	int het_field = -1;
	int max_field = -1;
	int chrom_field = -1;
	int pos_field = -1;
	int field_counter = 0;
	char slabel[1024];

	gzbuffer(snpfile,33554432);
	gzgets(snpfile,buffer,buff_len);
	// printf("\n About to start reading the file header \n\n");

	token = strtok(buffer,"\n\t ");
	while(token)
	{
	    if(strcmp(token,"alt") == 0)
	    {
		alt_field = field_counter;
		max_field = field_counter;
	    }
	    else
	    if(strcmp(token,"chrom") == 0)
	    {
		chrom_field = field_counter;
		max_field = field_counter;
	    }
	    else
	      if(strcmp(token,"pos") == 0)
		{
		  pos_field = field_counter;
		  max_field = field_counter;
		}
	      else
		if(strcmp(token,"type") == 0)
		  {
		    type_field = field_counter;
		    max_field = field_counter;
		  }
	    else
	      if(strcmp(token,"heterozygotes") == 0)
	      {
		  het_field = field_counter;
		  max_field = field_counter;
	      }
	      else
		if(strcmp(token,"homozygotes") == 0)
		  {
		    hom_field = field_counter;
		    max_field = field_counter;
		  }
		else
		  if(strcmp(token,"missingGenos") == 0)
		    {
		      missing_field = field_counter;
		      max_field = field_counter;
		    }
		  else
		    if(strcmp(token,"ref") == 0)
		      {
			ref_field = field_counter;
			max_field = field_counter;
		      }
		    else
		      if(strcmp(token,"gnomad.genomes.id") == 0)
			{
			  rs_field_gnomad_genomes = field_counter;
			  max_field = field_counter;
			}
				else
		      if(strcmp(token,"gnomad.exomes.id") == 0)
			{
			  rs_field_gnomad_exomes = field_counter;
			  max_field = field_counter;
			}
				else
		      if(strcmp(token,"dbSNP.name") == 0)
			{
			  rs_field_dbSNP = field_counter;
			  max_field = field_counter;
			}


	    field_counter++;
	    token = strtok(NULL,"\n\t ");
	    // token = strtok(NULL,"\n\t ");
	}
	if(ref_field < 0)
	  dump_error("Could not find the reference field");
	if(chrom_field < 0)
	  dump_error("Could not find the chromosome field");
	if(pos_field < 0)
	  dump_error("Could not find the position field");
	if(type_field < 0)
	  dump_error("Could not find the type field");
	if(alt_field < 0)
	  dump_error("Could not find the alt field");
	if(het_field < 0)
	  dump_error("Could not find the heterozygotes field");
	if(hom_field < 0)
	  dump_error("Could not find the heterozygotes field");
	if(missing_field < 0)
	  dump_error("Could not find the missing field");
	if(rs_field_gnomad_genomes < 0)
	  dump_error("Could not find the gnomad.genomes.id field");
	if(rs_field_gnomad_exomes < 0)
	  dump_error("Could not find the gnomad.exomes.id field");
	if(rs_field_dbSNP < 0)
	  dump_error("Could not find the dbSNP.name field");

	char **fields;

	if( (fields = malloc( (max_field+1)*sizeof(char *) )) == NULL)
	  dump_error("Can't allocate space for fields");


	gzgets(snpfile,buffer,buff_len);
	int len = strlen(buffer);
	char *chrom;
	char ref;
	char ref_string[8196];
	char *alt_a_temp;
	char alt_a_final[8196];
	char alt_count_final[8196];
	char snp_name[8196];
	int pos;
	int drop_it;

	char last_chr[128];
	int last_chr_no = 0;
	sprintf(last_chr,"!!!!!!");
        char token_copy[8192];
	int *allele_count;
	allele_count = ivector(0,100);
	// printf("\n About to start looping the file \n\n");
	while(len > 5)
	{
	    int i = 0;
	    int j = 0;

	    for(i=0;i<100;i++)
		allele_count[i] = 0;

	    char *last_field;
	    i = 0;
	    last_field = buffer;
	    while(i <= max_field)
	    {
		if(buffer[j] == '\t')
		{
		    buffer[j] = '\0';
		    fields[i] = last_field;
		    last_field = &buffer[j+1];
		    i++;
		}
		j++;
	    }

	    chrom = fields[chrom_field];
	    pos = atoi(fields[pos_field]);
	    ref = fields[ref_field][0];
	    strcpy(ref_string,fields[ref_field]);
	    alt_a_temp = fields[alt_field];

	    // Apparently even gnomad will give imprecise rs#
	    // For instance: http://gnomad.broadinstitute.org/dbsnp/rs373321043
	    // For now, just use chr_pos
	    // if(strcmp(fields[rs_field_gnomad_genomes],"!") != 0) {
	    // 	strcpy(snp_name, fields[rs_field_gnomad_genomes]);
	    // } else if(strcmp(fields[rs_field_gnomad_exomes],"!") != 0) {
	    // 	strcpy(snp_name, fields[rs_field_gnomad_exomes]);
	    // } else {
	    	// dbSNP will not give accurate rs# when alt != dbSNP.alleles (minor allele)
	    	// so use chr_pos
	    	// mappens most often in multiallelic cases, and makes it difficult
	    	// to correctly weigh variants in SKAT
	    	// In multiallelic cases with no rs#, we want to distinguish alleles
	      // if(strcmp(fields[type_field], "MULTIALLELIC") == 0) {
	      	sprintf(snp_name, "%s_%d_%s", chrom, pos, fields[alt_field]);
	      // } else {
	      	// sprintf(snp_name, "%s_%d", chrom, pos);
	      // }
	    	// strcpy(snp_name,fields[rs_field_dbSNP]);
		    // token = strtok(snp_name,"\t \n;,|/");
		    // if(token[0] != '!')
		    //   strcpy(snp_name,token);
		    // else {

		    // 	// In multiallelic cases with no rs#, we want to distinguish alleles
		    //   if(strcmp(fields[type_field], "MULTIALLELIC") == 0) {
		    //   	sprintf(snp_name, "%s_%d_%s", chrom, pos, fields[alt_field]);
		    //   } else {
		    //   	sprintf(snp_name, "%s_%d", chrom, pos);
		    //   }
		    // }
		  // }

	    token = fields[type_field];
	    if(token[0] == 'D')
	    {
		strcpy(token_copy,token+7);
		token = token_copy;
	    }
	    drop_it = FALSE;

	    if(strcmp(token,"LOW")==0)
	      drop_it = TRUE;
	    else
	      if(strcmp(token,"MESS")==0)
		drop_it = TRUE;


	    if(!drop_it)
	    {

	      // printf("\n I am not going to drop this SNP \n\n");
	      // printf("\n chr = %s pos = %d ref = %c ref_string = %s alt_a_temp = %s snp_name = %s\n\n",
	      //      chrom,pos,ref,ref_string,alt_a_temp,snp_name);
		for(i=0;i<tot_samples;i++)
		  strcpy(sorted_samples[i]->call,"0/0");
		// printf("\n Filled samples with hom ref \n\n");
		strcpy(slabel,"PASS");
		sprintf(ref_string,"%c",ref);

		// printf("\n ABout to process missing samples \n\n");
		process_samples(sorted_samples,fields[missing_field],"./.",tot_samples);
		strcpy(alt_a_final,alt_a_temp);

		// printf("\n Back from missing \n\n");

		if(strcmp(token,"SNP")==0)
		{
		    process_samples(sorted_samples,fields[het_field],"0/1",tot_samples);
		    process_samples(sorted_samples,fields[hom_field],"1/1",tot_samples);
		}
		else
		  if(strcmp(token,"MULTIALLELIC")==0)
		  {
		      int this_a = 1;
		      int this_a_pos = 0;
		      int has_del = FALSE;
		      char alt_a = 'N';
		      alt_a_final[0] = '\0';
		      int this_stop = strlen(alt_a_temp);
		      char this_geno[10];
		      j = 1;
		      while(fields[het_field])
		      {
			  sprintf(this_geno,"0/%d",j++);
			  fields[het_field] = process_samples(sorted_samples,fields[het_field],this_geno,tot_samples);
		      }
		      j = 1;
		      while(fields[hom_field])
		      {
			  sprintf(this_geno,"%d/%d",j,j);
			  // printf("\n About the call with j=%d fields = %s and geno = %s\n",j,fields[hom_field],this_geno);
			  j++;
			  fields[hom_field] = process_samples(sorted_samples,fields[hom_field],this_geno,tot_samples);

		      }

		      while(this_a_pos < this_stop)
		      {
			  if(alt_a_temp[this_a_pos] == ref)
			    this_a_pos+=2;
			  else
			    if(alt_a_temp[this_a_pos] == '+')
			    {
				alt_a = 'I';
				// printf("\n In here for an insertion with previous ac = %s and ref = %s\n",alt_a_final,ref_string);
				if(!has_del)
				{
				    if(this_a == 1)
				      sprintf(alt_a_final,"%c",ref);
				    else
				      sprintf(alt_a_final,"%s,%c",alt_a_final,ref);
				}
				else
				  sprintf(alt_a_final,"%s,%s",alt_a_final,ref_string);

				this_a_pos++;
				while( (this_a_pos < this_stop) && (alt_a_temp[this_a_pos] != '/'))
				{
				    if(isalpha(alt_a_temp[this_a_pos]))
				      sprintf(alt_a_final,"%s%c",alt_a_final,alt_a_temp[this_a_pos]);
				    this_a_pos++;
				}

				this_a_pos++;
				this_a++;
				sprintf(slabel,".");
			    }
			    else
			      if(alt_a_temp[this_a_pos] == '-')
			      {
				  alt_a = 'D';

				  if(strcmp(chrom,last_chr) != 0)
				  {
				      last_chr_no = -1;
				      for(i=0;i<no_contigs;i++)
					if(strcmp(chrom,contig_names[i])==0)
					{
					    strcpy(last_chr,chrom);
					    last_chr_no = i;
					    i = no_contigs;
					}
				      if(last_chr_no < 0)
				      {
					  printf("\n Failed to find chrom = %s \n",chrom);
					  exit(1);
				      }
				  }
				  pos--;

				  long this_offset = pos + contig_starts[last_chr_no]-1;
				  char sn[128];
				  has_del = TRUE;
				  ref = genome_buffer[this_offset];
				  this_a_pos++;
				  i = 0;
				  while((this_a_pos < this_stop) && (alt_a_temp[this_a_pos] != '/'))
				    sn[i++] = alt_a_temp[this_a_pos++];
				  sn[i] = '\0';
				  int del_len = atoi(sn);
				  del_len++;

				  char gb[4196];
				  strncpy(gb,&genome_buffer[this_offset],del_len);
				  gb[del_len] = '\0';
				  sprintf(ref_string,"%s",gb);
				  sprintf(slabel,".");
				  if(this_a == 1)
				    sprintf(alt_a_final,"%c",ref);
				  else
				  {
				      strcpy(gb,alt_a_final);
				      strcpy(sn,ref_string);
				      sn[1] = gb[0];
				      sprintf(alt_a_final,"%s",sn);
				      for(i=2,j=2;i<this_a;i++,j+=2)
				      {
				      	strcpy(sn,ref_string);
				      	sn[1] = gb[j];
				      	sprintf(alt_a_final,"%s,%s",alt_a_final,sn);
                                      }
				      sprintf(alt_a_final,"%s,%c",alt_a_final,ref);
				  }
				  this_a_pos++;
				  this_a++;
				  sprintf(slabel,".");
			      }
			      else
			      {
				  alt_a = alt_a_temp[this_a_pos];
				  if(this_a == 1)
				    sprintf(alt_a_final,"%c",alt_a);
				  else
				    sprintf(alt_a_final,"%s,%c",alt_a_final,alt_a);
				  this_a++;
				  this_a_pos+=2;
			      }
		      }
		  }
		  else
		    if(strcmp(token,"INS")==0)
		    {
			int mono_allelic = TRUE;
			int al = strlen(alt_a_temp);
			for(i=1;i<=al;i++)
			  if(alt_a_temp[i] == '/')
			    mono_allelic = FALSE;
			if(!mono_allelic)
			  sprintf(alt_a_final,"%c%s",ref,&alt_a_temp[3]);
			else
			  sprintf(alt_a_final,"%c%s",ref,&alt_a_temp[1]);
			sprintf(slabel,".");
			process_samples(sorted_samples,fields[het_field],"0/1",tot_samples);
			process_samples(sorted_samples,fields[hom_field],"1/1",tot_samples);

		    }
		    else // Deletion
		    {
			if(strcmp(chrom,last_chr) != 0)
			{
			    last_chr_no = -1;
			    for(i=0;i<no_contigs;i++)
			      if(strcmp(chrom,contig_names[i])==0)
			      {
				  strcpy(last_chr,chrom);
				  last_chr_no = i;
				  i = no_contigs;
			      }
			    if(last_chr_no < 0)
			    {
				printf("\n Failed to find chrom = %s \n",chrom);
				exit(1);
			    }
			}
			pos--;
			long this_offset = pos + contig_starts[last_chr_no]-1;
			ref = genome_buffer[this_offset];
			char sn[128];
			int mono_allelic = TRUE;
			int al = strlen(alt_a_temp);
			for(i=1;i<=al;i++)
			  if(alt_a_temp[i] == '/')
			    mono_allelic = FALSE;
			if(!mono_allelic)
			  sprintf(sn,"%s",&alt_a_temp[3]);
			else
			  sprintf(sn,"%s",&alt_a_temp[1]);

			int del_len = atoi(sn);
			del_len++;
			strncpy(ref_string,&genome_buffer[this_offset],del_len);
			ref_string[del_len] = '\0';
			sprintf(slabel,".");
			sprintf(alt_a_final,"%c",ref);
			process_samples(sorted_samples,fields[het_field],"0/1",tot_samples);
			process_samples(sorted_samples,fields[hom_field],"1/1",tot_samples);
		    }


		make_allele_count(all_samples,allele_count,tot_samples);
		int AN = 0;
		j = 0;
		for(i=0;i<100;i++)
		{
		    AN += allele_count[i];
		    if(i>0)
		    {
			if(allele_count[i] > 0)
			  {
			    if(j == 0)
			      sprintf(alt_count_final,"AC=%d",allele_count[i]);
			    else
			      sprintf(alt_count_final,"%s;AC=%d",alt_count_final,allele_count[i]);
			    j++;
			  }
			else
			  i = 100;
		    }
		}
		sprintf(alt_count_final,"%s;AN=%d",alt_count_final,AN);
		// printf("\n About to print a bunch of stuff \n\n");

		printf("\n%s\t%d\t%s\t%s\t%s\t.\t%s\t%s;NS=%d\tGT",
		       chrom,pos,snp_name,ref_string,alt_a_final,slabel,alt_count_final,tot_samples);
		for(i=0;i<tot_samples;i++)
		  printf("\t%s",all_samples[i]->call);

	    }

	    buffer[0] = '\0';
	    if(!gzeof(snpfile))
	      gzgets(snpfile,buffer,buff_len);
	    len = strlen(buffer);
	}
	printf("\n");

	exit(0);

}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
char *process_samples(SAMPLE_NODE **sorted_samples,char *fields,char *geno,int tot_samples)
{
  int j,last;
  SAMPLE_NODE *this;
  char old_char;


  // printf("\n Processing %s to give them genotype %s\n\n",fields,geno);

  j = 0;
  last = strlen(fields);
  while(j < last)
  {
    // printf("\n j = %d char = %c",j,fields[j]);
      if( (fields[j] == ';') || (fields[j] == '/') )
      {
	  old_char = fields[j];
	  fields[j] = '\0';
	  // printf("\n Got here with fields = %s geno = %s\n\n",fields,geno);
	  if(fields[0] == '!')
	    return &fields[j+1];

	  this = *((SAMPLE_NODE **)bsearch(fields,sorted_samples,tot_samples,sizeof(SAMPLE_NODE *),search_sample));
	  if(!this)
	  {
	      printf("\n Failed to find the name of %s in the sample list \n",fields);
	      exit(1);
	  }
	  // printf("\n Previous call was %s sample is %s\n\n",this->call,this->name);
	  if( (this->call[2] == '0') || (this->call[0] == '.'))
	    strcpy(this->call,geno);
	  else
	  {
	      int k,l,m,n;
	      char ss[80];
	      sscanf(geno,"%d/%d",&k,&l);
	      sscanf(this->call,"%d/%d",&m,&n);
	      sprintf(ss,"%d/%d",l,n);
	      strcpy(this->call,ss);
	  }
	  if(old_char == '/')
	      return &fields[j+1];

	  fields = &fields[j+1];
	  j = 0;
	  last = strlen(fields);
      }
      j++;
  }

  // printf("\n made it out with fields = %s \n\n",fields);
  if(fields[0] == '!')
    return NULL;
  this =  *((SAMPLE_NODE **) bsearch(fields,sorted_samples,tot_samples,sizeof(SAMPLE_NODE *),search_sample));
  if(!this)
  {
      printf("\n Failed to find the name of %s in the sample list \n",fields);
      exit(1);
  }
  // printf("\n Found %s with previous call %s\n\n",this->name,this->call);
  if( (this->call[2] == '0') || (this->call[0] == '.'))
    strcpy(this->call,geno);
  else
  {
      int k,l,m,n;
      char ss[80];
      sscanf(geno,"%d/%d",&k,&l);
      sscanf(this->call,"%d/%d",&m,&n);
      sprintf(ss,"%d/%d",l,n);
      strcpy(this->call,ss);
  }

  return NULL;
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
SAMPLE_NODE *sample_alloc(char *s,int n)
{
  SAMPLE_NODE *sn;

  if((sn = malloc(sizeof(SAMPLE_NODE))) == NULL)
    {
      printf("\n Failed to allocated sample node for %s which is number %d\n",
	     s,n);
      exit(1);
    }
  sn->name = cvector(0,4196);
  strcpy(sn->name,s);
  sn->number = n;
  sn->call = cvector(0,32);
  return sn;
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
void make_allele_count(SAMPLE_NODE **sn,int *count, int tot_samples)
{
  int i,j,k;
  for(i=0;i<tot_samples;i++)
  {
      if(sn[i]->call[0] != '.')
      {
	  sscanf(sn[i]->call,"%d/%d",&j,&k);
	  // printf("\n For sample %s we have genomeype %s broken down into %d and %d",sn[i]->name,sn[i]->call,j,k);
	  count[j]++;
	  count[k]++;
      }
  }
  return;
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */
inline int qual_to_phred(double qual)
{
	if(qual > 0.999999)
		return 60;

	if(qual < 1e-9)
		return 0;

	return floor(-10 * log(1.0 - qual));
}
/*-------------------------------------------------------------------------------------------------------------------------------------- */


char *cvector(int nl,int nh)
{
	char *v;

	v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
	if (!v) dump_error("allocation failure in cvector()");
	return v-nl;
}

unsigned char *ucvector(int nl,int nh)
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned) (nh-nl+1)*sizeof(unsigned char));
	if (!v) dump_error("allocation failure in ucvector()");
	return v-nl;
}

unsigned long long *ullvector(int nl,int nh)
{
	unsigned long long *v;

	v=(unsigned long long *)malloc((unsigned) (nh-nl+1)*sizeof(unsigned long long));
	if (!v) dump_error("allocation failure in ullvector()");
	return v-nl;
}

int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) dump_error("allocation failure in ivector()");
	return v-nl;
}

unsigned int *uvector(int nl,int nh)
{
	unsigned int *v;

	v=(unsigned int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) dump_error("allocation failure in ulvector()");
	return v-nl;
}
long *lvector(int nl,int nh)
{
	long *v;

	v=(long *)malloc((nh-nl+1)*sizeof(long));
	if (!v) dump_error("allocation failure in lvector()");
	return v-nl;
}

unsigned long *ulvector(int nl,int nh)
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned) (nh-nl+1)*sizeof(long));
	if (!v) dump_error("allocation failure in ulvector()");
	return v-nl;
}

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) dump_error("allocation failure in dvector()");
	return v-nl;
}


int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) dump_error("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) dump_error("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned int **umatrix(int nrl,int nrh,int ncl,int nch)
{
        int i;
	unsigned int **m;

	m=(unsigned int **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned int*));
	if (!m) dump_error("allocation failure 1 in ulmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned int *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned int));
		if (!m[i]) dump_error("allocation failure 2 in ulmatrix()");
		m[i] -= ncl;
	}
	return m;
}

unsigned long **ulmatrix(int nrl,int nrh,int ncl,int nch)
{
        int i;
	unsigned long **m;

	m=(unsigned long **)malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned long*));
	if (!m) dump_error("allocation failure 1 in ulmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned long *)malloc((unsigned) (nch-ncl+1)*sizeof(unsigned long));
		if (!m[i]) dump_error("allocation failure 2 in ulmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_umatrix(unsigned int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}


void free_ulmatrix(unsigned long **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}


void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}


double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) dump_error("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) dump_error("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}

char **cmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	char **m;

	m=(char **) malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(char *) malloc((unsigned) (nch-ncl+1)*sizeof(char));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}


unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	unsigned char **m;

	m=(unsigned char **) malloc((unsigned) (nrh-nrl+1)*sizeof(unsigned char*));
	if (!m) dump_error("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(unsigned char *) malloc((unsigned) (nch-ncl+1)*sizeof(unsigned char));
		if (!m[i]) dump_error("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}



void free_cmatrix(char **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}


void free_ucmatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((void *) (m[i]+ncl));
	free((void *) (m+nrl));
}

void free_cvector(char *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_ucvector(unsigned char *v, int nl, int nh)
{
	free((void *) (v+nl));
}


void free_ivector(int *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_ulvector(unsigned long *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_uvector(unsigned int *v, int nl, int nh)
{
	free((void *) (v+nl));
}

void free_dvector(double *v, int nl,int nh)
{
	free((void *) (v+nl));
}

/*---------------------------------------------------------------------*/

int sort_sample(const void *a,const void *b)
{
  SAMPLE_NODE *as,*bs;

  as = *((SAMPLE_NODE **)a);
  bs = *((SAMPLE_NODE **)b);

  return strcmp(as->name,bs->name);
}
/*---------------------------------------------------------------------*/

int search_sample(const void *a,const void *b)
{
  SAMPLE_NODE *bs;
  char *as;

  as = ((char *)a);
  bs = *((SAMPLE_NODE **)b);
  // printf("\n Comparing %s to %s \n\n",as,bs->name);

  return strcmp(as,bs->name);
}

/*---------------------------------------------------------------------*/

void dump_error(char *error_text)
{

	fprintf(outfile,"Seqant_to_Vcf error...\n");
	fprintf(outfile,"%s\n",error_text);
	fprintf(outfile,"...now exiting to system...\n");
	exit(1);
}

/*---------------------------------------------------------------------*/


