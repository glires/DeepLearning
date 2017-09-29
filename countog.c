/*                                                                           */
/* countog.c - prepare training and test data by counting oligonucleotides   */
/*                                                                           */
/* COMPILE                                                                   */
/*   $ gcc -W -Wall -O -ansi -pedantic -Werror -lm -o countog countog.c      */
/*                                                                           */
/* SYNOPSIS                                                                  */
/*   $ countog [-c number_of_oligos] [-d] [-g genome_size] [-l label] \      */
/*       [-o size_of_oligo] [-q min_q_score] [-r] [-s size_of_shift]  \      */
/*       input_FASTA_or_FASTQ                                                */
/*                                                                           */
/* USAGE                                                                     */
/*   $ countog -d -l mouse -o 6 -t 100 mm10.fa                               */
/*                                                                           */
/* DESCRIPTION                                                               */
/*    This program reads a FASTA or FASTQ sequence file and counts           */
/*    numbers of each specified-length oligonucleotides.                     */
/*    Normalized values, which range from 0 to 1, are printed onto the       */
/*    standard output.                                                       */
/*                                                                           */
/* OPTIONS                                                                   */
/*   -c  Number of counting oligos for one-line data (default 100000)        */
/*   -d  Print the header line                                               */
/*   -g  Maximum genome size (default: 4294967296)                           */
/*   -l  Add a label for training data                                       */
/*   -o  Size of oligonucleotide in nt                                       */
/*   -q  Minimum quality score (default: 16)                                 */
/*   -r  Merge complementary oligonucleotides                                */
/*   -s  Size of shift in bp for the next round                              */
/*   -t  Number of one-line data (default: 20000)                            */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   2017-05-29  Started to be coded                                         */
/*   2017-09-15  Named as countog, count oligonucleotides                    */
/*   2017-09-16  Fixed a bug which is related to reading FASTQ               */
/*   2017-09-18  Output normalized counts                                    */
/*   2017-09-19  Support option -l, label for training data                  */
/*   2017-09-29  Released at GitHub                                          */
/*                                                                           */
/* MEMORANDOM                                                                */
/*   Next error code: Error 9, Error 13, Error 14, ...                       */
/*                                                                           */

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

#define OLIGO 8
#define SIZE_GENOME 4294967296L	/* 2^32, more than 4 billion (bases) */
#define SIZE_LINE_CHARS 1024
#define SIZE_COUNTING 100000
#define SIZE_DATA 20000
#define SIZE_SHIFT 20000
#define NUCLEOTIDES 4
#define DEFAULT_MIN_QSCORE 16
#define CODE_TO_SCORE (int)(-33)

extern char *optarg;
extern int optind;

int *counter, *complementary;
char *genome, *genomep;	/* genome and position */

long int size_genome    = SIZE_GENOME,
         gsize          = 0,	/* exclude inserted ns */
         gnsize         = 0;	/* include inserted ns */
int      size_oligo     = 1,
         size_counting  = SIZE_COUNTING,
         oligo          = OLIGO,
         size_data      = SIZE_DATA,
         size_shift     = SIZE_SHIFT,
         minimum_qscore = DEFAULT_MIN_QSCORE;
short int reduce = 0,	/* for complementary oligos */
          header = 0,	/* print the header line */
          label  = 0;	/* label for training data */


int getopt(int, char * const [], const char *);


int reset_counter(void)
{
  int i;
  for (i = 0; i < size_oligo; i++) { counter[i] = 0; }
  return i;
}


int print_header(void)
{
  int i, j, digit, fwd;

  if (header == 0) { return (int)header; }
  if (label !=0 ) { fprintf(stdout, "DATA\t"); }

  for (j = 0; j < size_oligo; j++)
  {
    if (j > 0) { fputc('\t', stdout); }
    fwd = j;
    for (i = 0; i < oligo; i++)
    {
      digit = fwd % NUCLEOTIDES;
      fwd = (int)(fwd / NUCLEOTIDES);
      if      (digit == 0) { fputc('T', stdout); }
      else if (digit == 1) { fputc('C', stdout); }
      else if (digit == 2) { fputc('A', stdout); }
      else if (digit == 3) { fputc('G', stdout); }
      else
      {
        fprintf(stderr, "Error 8: nucleotide %d\n", digit);
        return EXIT_FAILURE;
      }
    }
  }
  fputc('\n', stdout);
  return oligo;
}


int count_octamer(void)
{	/* not necessarily restrict oligomer to octamer */
  int i, n, octa = oligo, idx = 0, idc = 0;

  for (i = 0; i < octa; i++)
  {
    switch (*(genomep + i))
    {
      case 't':  n = 0; break;
      case 'c':  n = 1; break;
      case 'a':  n = 2; break;
      case 'g':  n = 3; break;
      case '\0': return -1;
      default:   return i;	/* i < octa */
    }
    idx += n * (int)pow((double)NUCLEOTIDES, (double)i);
  }
  counter[idx]++;
  if (complementary[idx] == -1)	/* build the complementary table */
  {
    for (i = 0; i < octa; i++)
    {
      switch (*(genomep + octa -1 - i))
      {
        case 't':  n = 2; break;
        case 'c':  n = 3; break;
        case 'a':  n = 0; break;
        case 'g':  n = 1; break;
        default:   fprintf(stderr, "Error 3: complementary (%s)\n", genomep);
                   return EXIT_FAILURE;
      }
      idc += n * (int)pow((double)NUCLEOTIDES, (double)i);
    }
    complementary[idx] = idc;
  }
  return i;	/* i == octa */
}


int increment_counter(int upto)
{
  int i = 0, rv, counter_shift = 1;

  while (i < upto)
  {
    rv = count_octamer();
    if (rv == -1) { genomep = genome + size_shift * counter_shift++; }
    else if (rv == oligo) { i++; }
    if (gsize >= size_shift * (counter_shift + 1)) { counter_shift = 0; }
    genomep++;
  }
  return i;
}


int get_complementary_oligo(int forward)
{
  int fwd = forward, rev = 0, i, n;
  int *digits;

  digits = (int *)malloc(sizeof(int) * oligo);
  if (digits == NULL)
  {
    fprintf(stderr, "Error 5: malloc for digits\n");
    return EXIT_FAILURE;
  }
  for (i = 0; i < oligo; i++)
  { digits[i] = fwd % NUCLEOTIDES; fwd = (int)(fwd / NUCLEOTIDES); }
  for (i = oligo - 1; i > -1; i--)
  {
    if      (digits[i] == 0) { n = 2; }
    else if (digits[i] == 1) { n = 3; }
    else if (digits[i] == 2) { n = 0; }
    else if (digits[i] == 3) { n = 1; }
    else
    {
      fprintf(stderr, "Error 4: nucleotide %d\n", digits[i]);
      free(digits);
      return EXIT_FAILURE;
    }
    rev += n * (int)pow((double)NUCLEOTIDES, (double)(oligo - i - 1));
  }
  free(digits);
  return rev;
}


int output_normalized_counts(char *tlabel)
{
  int i, j = 0, max = 0, *total;	/* j is a counter for output values */

  reset_counter();
  increment_counter(size_counting);
  if (label != 0) { fprintf(stdout, "%s\t", tlabel); }

  if (reduce == 0)
  {
    for (i = 0; i < size_oligo; i++)	/* search for maximum counts */
    { if (counter[i] > max) { max = counter[i]; } }
    for (i = 0; i < size_oligo; i++)
    {
      if (j++ != 0) { fputc('\t', stdout); }
      fprintf(stdout, "%.4f", (float)counter[i] / max);
    }
  }
  else	/* merge complementary oligos */
  {
    total = (int *)malloc(sizeof(int) * size_oligo);
    for (i = 0; i < size_oligo; i++)
    {
      if (counter[i] != -1)	/* ignore complementary already processed */
      {
        if (complementary[i] == -1)
          complementary[i] = get_complementary_oligo(i);
        assert(counter[complementary[i]] != -1);
        total[j++] = counter[i] + counter[complementary[i]];
        /* ignore complementary in the next round */
        counter[i] = counter[complementary[i]] = -1;      /* to ignore */
      }
    }	/* j is the number of the total elements */
    for (i = 0; i < j; i++)	/* search for one of the maximum counts */
    { if (total[i] > max) { max = total[i]; } }
    for (i = 0; i < j; i++)
    {
      if (i != 0) { fputc('\t', stdout); }
      fprintf(stdout, "%.4f", (float)total[i] / max);
    }
    free(total);
  }
  fputc('\n', stdout);
  return j;
}


int main(int argc, char* argv[])
{
  FILE *check;
  char line[SIZE_LINE_CHARS], qscore[SIZE_LINE_CHARS], tlabel[SIZE_LINE_CHARS];
  int num_chars;
  int basepairs;
  int i, opt, fastq = -1;

  assert(sizeof(int) >= 4);	/* int should be no less than 32 bit */

  while ((opt = getopt(argc, argv, "c:dg:l:o:q:rs:t:")) != -1)
  {
    switch (opt)
    {
      case 'c': size_counting = atoi(optarg);
                break;
      case 'd': header = 1;	/* print the header line */
                break;
      case 'g': size_genome = atol(optarg);
                break;
      case 'l': strcpy(tlabel, optarg); label = 1;
                break;
      case 'o': oligo = atoi(optarg);
                break;
      case 'q': minimum_qscore = atoi(optarg);
                break;
      case 'r': reduce = 1;	/* merge complementary oligos */
                break;
      case 's': size_shift = atoi(optarg);
                break;
      case 't': size_data = atoi(optarg);
                break;
      default:  fprintf(stderr, "Warning: unknown option -%c\n", opt);
    }
  }

  if (optind + 1 != argc)
  {
    fprintf(stderr, "Error 1: specify an input FASTA file name\n");
    return EXIT_FAILURE;
  }
  else
  {
    check = fopen(argv[optind], "r");
    if (check != NULL)
    { fclose(check); check =freopen(argv[optind], "r", stdin); }
  }

  genome = (char *)malloc(size_genome);
	/* char genome[size_genome]; does not work. */
  if (genome == NULL)
  {
    fprintf(stderr, "Error 2: malloc(size_genome)\n");
    return EXIT_FAILURE;
  }
  genomep = genome;

  for (i = 0; i < oligo; i++) { size_oligo *= NUCLEOTIDES; }
  	/* T, C, A, and G */
  counter = (int *)malloc(sizeof(int) * size_oligo);
  complementary = (int *)malloc(sizeof(int) * size_oligo);
  for (i = 0; i < size_oligo; i++) { complementary[i] = -1; }

  if (fgets(line, SIZE_LINE_CHARS, check) == NULL)
	/* read the first line */
  { fprintf(stderr, "Error 6: fgets()\n"); return EXIT_FAILURE; }
  else
  {
    if (line[0] == '>')      fastq = 0;
    else if (line[0] == '@') fastq = 1;
    else
    {
      fprintf(stderr, "Error 7: neither FASTA nor FASTQ\n");
      return EXIT_FAILURE; 
    }
  }

  do	/* build the reference sequence using the genomep pointer */
  {
    if (fastq == 1)	/* FASTQ */
    {
      assert(line[0] == '@');
      *genomep++ = 'n';	/* insert n to split the two scaffolds */
      gnsize++;
      if (fgets(line, SIZE_LINE_CHARS, check) == NULL)
      { fprintf(stderr, "Error 10: fgets()\n"); return EXIT_FAILURE; }
      if (fgets(qscore, SIZE_LINE_CHARS, check) == NULL)
      { fprintf(stderr, "Error 11: fgets()\n"); return EXIT_FAILURE; }
      if (fgets(qscore, SIZE_LINE_CHARS, check) == NULL)
      { fprintf(stderr, "Error 12: fgets()\n"); return EXIT_FAILURE; }
      basepairs = (int)strlen(line) - 1;
      assert(basepairs == (int)strlen(qscore) - 1);
      if (size_genome - 1 < gnsize + (long int)basepairs) break;
      for (i = 0; i < basepairs; i++)
      {
        if ((int)qscore[i] + CODE_TO_SCORE < minimum_qscore)
        { *genomep++ = 'n'; }
        else if (line[i] == 'T') { *genomep++ = 't'; }
        else if (line[i] == 'C') { *genomep++ = 'c'; }
        else if (line[i] == 'A') { *genomep++ = 'a'; }
        else if (line[i] == 'G') { *genomep++ = 'g'; }
        else                     { *genomep++ = line[i]; }
      }
      gsize  += (long int)basepairs;
      gnsize += (long int)basepairs;
    }
    else	/* FASTA */
    {
      if (line[0] == '>') { *genomep++ = 'n'; gnsize++; continue; }
	/* insert n to split the two scaffolds */
      basepairs = num_chars = (int)strlen(line);
      if (size_genome - 1 < gnsize + (long int)basepairs) break;
      for (i = 0; i < num_chars; i++)
      {
        if (isalpha(line[i]) == 0) { basepairs--; }
        else if (line[i] == 'T') { *genomep++ = 't'; }
        else if (line[i] == 'C') { *genomep++ = 'c'; }
        else if (line[i] == 'A') { *genomep++ = 'a'; }
        else if (line[i] == 'G') { *genomep++ = 'g'; }
        else                     { *genomep++ = line[i]; }
      }
      gsize  += (long int)basepairs;
      gnsize += (long int)basepairs;
    }
  } while (fgets(line, SIZE_LINE_CHARS, check) != NULL);
  *genomep = '\0';

  genomep = genome;	/* reset */
  if (gsize < (long int)size_shift) { size_shift = 1; }
  print_header();
  for (i = 0; i < size_data; i++) output_normalized_counts(tlabel);

  fclose(check);
  free(complementary);
  free(counter);
  free(genome);
  return EXIT_SUCCESS;
}
