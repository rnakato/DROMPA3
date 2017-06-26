/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdlib.h>
#include <string.h>
#include "outputfile.h"
#include "filehandle.h"
#include "stringp.h"
#include "alloc.h"
#include "macro.h"
#include "compress.h"

void output_bindata(char *dir, char *prefix, RefGenome *g, TYPE_WIGARRAY *array, char *gtfile, int binsize, int binnum, int chr, PWfile_Type wtype){
  char *output_prefix = alloc_str_new(dir, strlen(prefix) + 1024);
  if(wtype==TYPE_BEDGRAPH || wtype==TYPE_BIGWIG) sprintf(output_prefix, "%s/%s.%d", dir, prefix, binsize);
  else sprintf(output_prefix, "%s/%s_%s.%d", dir, prefix, g->chr[chr].name, binsize);
  char *outputfilename = alloc_str_new(output_prefix, 200);

  if(wtype==TYPE_BINARY){
    sprintf(outputfilename, "%s.bin", output_prefix);
    make_binary(array, outputfilename, binnum);
  }
  else if(wtype==TYPE_BEDGRAPH || wtype==TYPE_BIGWIG){
    sprintf(outputfilename, "%s.bedGraph", output_prefix);
    make_bedGraph(g, array, outputfilename, output_prefix, binsize, binnum, chr);
    if(chr == g->chrnum-1 && wtype==TYPE_BIGWIG) convert_bedGraph_to_bigWig(outputfilename, output_prefix, gtfile);
  }
  else if(wtype==TYPE_COMPRESSWIG || wtype==TYPE_UNCOMPRESSWIG){
    sprintf(outputfilename, "%s.wig", output_prefix);
    make_wig(g, array, outputfilename, prefix, binsize, binnum, chr, wtype);
  }

  MYFREE(output_prefix);
  MYFREE(outputfilename);
  return;
}

void make_binary(TYPE_WIGARRAY *array, char *outputfile, int binnum){
  FILE *OUT = my_fopen(outputfile, FILE_MODE_WB);
  if(fwrite(array, sizeof(TYPE_WIGARRAY) * binnum, 1, OUT) != 1){
    fprintf(stderr,"[E] fwrite error:%s\n", outputfile);
    exit(1);
  }
  fclose(OUT);
  return;
}

void make_bedGraph(RefGenome *g, TYPE_WIGARRAY *array, char *outputfile, char *prefix, int binsize, int binnum, int chr){
  int i,e;
  if(chr==1) remove_file(outputfile);
  FILE *OUT = my_fopen(outputfile, FILE_MODE_A);
  for(i=0; i<binnum; i++){
    if(i==binnum -1) e = g->chr[chr].len-1; else e = (i+1)*binsize;
    if(array[i]) fprintf(OUT, "%s %d %d %.3f\n", g->chr[chr].name, i*binsize, e, WIGARRAY2VALUE(array[i]));
  }
  fclose(OUT);

  if(chr == g->chrnum-1){
      char *command = alloc_str_new(outputfile, strlen(outputfile) + 1024);
      sprintf(command, "sort -k1,1 -k2,2n %s > %s.sort", outputfile, outputfile);
      LOG("%s\n", command);
      my_system(command);
      char *tempfile = alloc_str_new(outputfile, 1024);
      sprintf(tempfile, "%s.temp", outputfile);
      OUT = my_fopen(tempfile, FILE_MODE_A);
      fprintf(OUT, "browser position %s:%ld-%ld\n", g->chr[1].name, g->chr[1].len/3, min(g->chr[1].len/3+1000000, g->chr[1].len-1));
      fprintf(OUT, "browser hide all\n");
      fprintf(OUT, "browser pack refGene encodeRegions\n");
      fprintf(OUT, "browser full altGraph\n");
      fprintf(OUT, "track type=bedGraph name=\"%s\" description=\"Merged tag counts for every %d bp\" visibility=full\n", prefix, binsize);
      fclose(OUT);
      sprintf(command, "cat %s %s.sort > %s", tempfile, outputfile, outputfile);
      LOG("%s\n", command);
      my_system(command);
      sprintf(command, "rm %s.sort %s", outputfile, tempfile);
      LOG("%s\n", command);
      my_system(command);

      MYFREE(command);
      MYFREE(tempfile);
  }
  return;
}

void convert_bedGraph_to_bigWig(char *outputfile, char *output_prefix, char *gtfile){
  printf("convert to bigWig format...\n");
  char *command = alloc_str_new(outputfile, strlen(gtfile) + strlen(output_prefix) + 1024);
  sprintf(command, "bedGraphToBigWig %s.bedGraph %s %s.bw", output_prefix, gtfile, output_prefix);
  LOG("%s\n", command);
  my_system(command);
  remove(outputfile);
  MYFREE(command);
  return;
}

void make_wig(RefGenome *g, TYPE_WIGARRAY *array, char *outputfile, char *name, int binsize, int binnum, int chr, PWfile_Type wtype){
  int i;
  FILE *OUT = my_fopen(outputfile, FILE_MODE_WRITE);

  /* first line */
  fprintf(OUT, "track type=wiggle_0\tname=\"%s_%s\"\tdescription=\"Merged tag counts for every %d bp\"\n", name, g->chr[chr].name, binsize);

  /* second line */
  if(!strcmp(g->chr[chr].name, "2micron")) fprintf(OUT, "variableStep\tchrom=%s\tspan=%d\n", g->chr[chr].name, binsize);
  else fprintf(OUT, "variableStep\tchrom=%s\tspan=%d\n", g->chr[chr].name, binsize);

  /* data line */
  for(i=0; i<binnum; i++){
    if(array[i]) fprintf(OUT, "%d\t%.3f\n", i*binsize +1, WIGARRAY2VALUE(array[i]));
  }
  fclose(OUT);

  /* compression */
  if(wtype==TYPE_COMPRESSWIG) compress2gz(outputfile);

  return;
}
