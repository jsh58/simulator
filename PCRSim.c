/*
  John Gaspar
  Dec. 2015

  Finding PCR primer matches in a genome.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "PCRSim.h"
#include "jmg_utils.h"

// global variables
static char* line;

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./PCRSim {%s <file> %s <file>", GENFILE, PRIMFILE);
  fprintf(stderr, " %s <file>} [optional parameters]\n", OUTFILE);
  fprintf(stderr, "Required parameters:\n");
  fprintf(stderr, "  %s  <file>       Fasta file of reference genome\n", GENFILE);
  fprintf(stderr, "  %s  <file>       Input file listing primer sequences, one set (forward and\n", PRIMFILE);
  fprintf(stderr, "                     reverse) per line, comma- or tab-delimited. For example:\n");
  fprintf(stderr, "                       341F-926R,CCTACGGGAGGCAGCAG,AAACTCAAAKGAATTGACGG\n");
  fprintf(stderr, "                       968F-1401R,AACGCGAAGAACCTTAC,CGTTCCCGGGCCTTGTACACACCG\n");
  fprintf(stderr, "                     The sequences can contain IUPAC ambiguous DNA codes.\n");
  fprintf(stderr, "                     NOTE: Both primers should be given with respect to the plus\n");
  fprintf(stderr, "                       strand, i.e. the sequence given for the reverse primer is\n");
  fprintf(stderr, "                       the reverse-complement of actual reverse primer\n");
  fprintf(stderr, "  %s  <file>       Output file for primer-genome matches\n", OUTFILE);
  fprintf(stderr, "Optional parameters:\n");
  fprintf(stderr, "  %s  <int>        Minimum amplicon length (def. %d)\n", MINLEN, DEFMIN);
  fprintf(stderr, "  %s  <int>        Maximum amplicon length (def. %d)\n", MAXLEN, DEFMAX);
  fprintf(stderr, "  %s  <float>      Minimum primer-genome match score (in (0-1]; def. %.2f)\n", MINSCORE, DEFSCORE);

  fprintf(stderr, "  %s  <file>       Log file for stitching results\n", LOGFILE);
  fprintf(stderr, "  %s               Option to check for dovetailing of the reads\n", DOVEOPT);
  fprintf(stderr, "  %s <file>       Log file for dovetailed reads only\n", DOVEFILE);
  fprintf(stderr, "  %s               Option to produce shortest stitched read, given\n", MAXOPT);
  fprintf(stderr, "                     multiple overlapping possibilities (by default,\n");
  fprintf(stderr, "                     the longest stitched read is produced)\n");
  fprintf(stderr, "  %s              Option to print counts of stitching results to stdout\n", VERBOSE);
  fprintf(stderr, "  %s              Option to print counts of stitching results to stdout\n", VERBOSE);
  exit(-1);
}


/* void freeMemory()
 * Frees allocated memory.
 */
void freeMemory(Primer* head) {
  Primer* temp;
  for (Primer* p = head; p != NULL; ) {
    free(p->name);
    for (int i = 0; i < 4; i++)
      free(p->seq[i]);
    temp = p;
    p = p->next;
    free(temp);
  }
}


/* int getSeq()
 * Read sequence and quality scores from a fastq file.
 */
int getSeq(FILE* in, char* line, char* seq,
    char* qual, int nSeq, int nQual) {
  for (int i = 0; i < 3; i++)
    getLine(line, MAX_SIZE, in);
  int len = strlen(seq);
  if (len != strlen(qual))
    exit(error("", ERRQUAL));
  return len;
}

/* float compare()
 * Compare two sequences. Return the percent mismatch.
 */
float compare(char* seq1, char* seq2, int length,
    float mismatch, int overlap) {
  int mis = 0;       // number of mismatches
  int len = length;  // length of overlap, not counting Ns
  float allow = len * mismatch;
  for (int i = 0; i < length; i++) {
    // do not count Ns
    if (seq1[i] == 'N' || seq2[i] == 'N') {
      if (--len < overlap || mis > len * mismatch)
        return 0;
      allow = len * mismatch;
    } else if (seq1[i] != seq2[i] && ++mis > allow)
      return 0;
  }
  return (float) mis / len;
}

/* int findPos()
 * Find optimal overlapping position.
 */
int findPos (char* seq1, char* seq2, char* qual1,
    char* qual2, int len1, int len2, int overlap,
    int dovetail, float mismatch, int maxLen,
    float* best) {
  int pos = len1 - overlap + 1;  // position of match
  for (int i = len1 - overlap; i > -1; i--) {
    if (len1 - i > len2 && !dovetail)
      break;
    float res = compare(seq1 + i, seq2,
      len1-i < len2 ? len1-i : len2, mismatch, overlap);
    if (res < *best || (res == *best && !maxLen)) {
      *best = res;
      pos = i;
    }
    if (res == 0.0f && maxLen)
      return pos;  // shortcut for exact match
  }

  // check for dovetailing
  if (dovetail) {
    for (int i = 1; i < len2 - overlap + 1; i++) {
      float res = compare(seq1, seq2 + i,
        len2-i < len1 ? len2-i : len1, mismatch, overlap);
      if (res < *best || (res == *best && !maxLen)) {
        *best = res;
        pos = -i;
      }
      if (res == 0.0f && maxLen)
        return pos;  // shortcut for exact match
    }
  }

  return pos;
}

/* void createSeq()
 * Create stitched sequence (into seq1, qual1).
 */
void createSeq(char* seq1, char* seq2, char* qual1, char* qual2,
    int len1, int len2, int pos) {
  int len = len2 + pos;  // length of stitched sequence
  for (int i = 0; i < len; i++) {
    if (i - pos < 0)
      continue;
    // disagreements favor higher quality score or
    //   equal quality score that is closer to 5' end
    else if (i >= len1 ||
        (seq1[i] != seq2[i-pos] && (qual1[i] < qual2[i-pos] ||
        (qual1[i] == qual2[i-pos] && i >= len2 - i + pos)))) {
      seq1[i] = seq2[i-pos];
      qual1[i] = qual2[i-pos];
    } else if (qual1[i] < qual2[i-pos])
      qual1[i] = qual2[i-pos];
  }
  seq1[len] = '\0';
  qual1[len] = '\0';
}

/* void printRes()
 * Print stitched read.
 */
void printRes(FILE* out, FILE* log, FILE* dove,
    char* header, char* seq1, char* seq2, char* qual1,
    char* qual2, int len1, int len2, int pos, float best) {
  // log result
  if (log != NULL) {
    fprintf(log, "%s\t%d\t%d\t", header,
      pos < 0 ? (len2+pos < len1 ? len2+pos : len1) :
      (len1-pos < len2 ? len1-pos : len2),
      len2 + pos);
    best ? fprintf(log, "%.3f", best) : fprintf(log, "0");
    fprintf(log, "\n");
  }

  // log 3' overhangs of dovetailed sequence(s)
  if (dove != NULL && (len1 > len2 + pos || pos < 0)) {
    fprintf(dove, "%s\t%s\t", header, len1 > len2 + pos ?
      seq1 + len2 + pos : "-");
    if (pos < 0)
      for (int i = -1; i - pos > -1; i--)
        fprintf(dove, "%c", comp(seq2[i - pos]));
    else
      fprintf(dove, "-");
    fprintf(dove, "\n");
  }

  // print stitched sequence
  createSeq(seq1, seq2, qual1, qual2, len1, len2, pos);
  fprintf(out, "@%s\n%s\n+\n%s\n", header, seq1, qual1);
}

/* void printFail()
 * Print stitch failure reads.
 */
void printFail(FILE* un1, FILE* un2, FILE* log,
    char* header, char* head1, char* head2, char* seq1,
    char* seq2, char* qual1, char* qual2, int len) {
  if (log != NULL)
    fprintf(log, "%s\tn/a\n", header);
  if (un1 != NULL && un2 != NULL) {
    fprintf(un1, "@%s\n%s\n+\n%s\n", head1, seq1, qual1);
    // put rev sequence back
    fprintf(un2, "@%s\n", head2);
    for (int i = len - 1; i > -1; i--)
      fprintf(un2, "%c", comp(seq2[i]));
    fprintf(un2, "\n+\n");
    for (int i = len - 1; i > -1; i--)
      fprintf(un2, "%c", qual2[i]);
    fprintf(un2, "\n");
  }
}

/* int getChunk()
 * Loads the next chunk from the genome.
 */
int getChunk(char* chunk, int pos, FILE* gen) {
  int i = 0;

  // copy piece from 3' end
  if (chunk[0] != '\0')
    while (chunk[pos] != '\0')
      chunk[i++] = chunk[pos++];

  // copy next chunk of genome
  while (i < CHUNK_SIZE) {
    int ch = getc(gen);
    if (ch == EOF) {
      // end of file
      chunk[i] = '\0';
      return 0;
    } else if (ch == '>') {
      // next chromosome
      chunk[i] = '\0';
      getLine(line + 1, MAX_SIZE, gen);
      return 0;
    } else if (ch != '\n')
      chunk[i++] = toupper(ch);
  }
  chunk[i] = '\0';
  return 1;
}


/* void findMatch()
 * Find a primer match to the genome chunk.
 */
void findMatch(Primer* p, char* chunk, int pos,
    Match* dummy, int minLen, int maxLen, float minScore) {
  // check for first primer match (p->seq[0] or p->seq[3])
  for 
  // p->fmax, p->rmax are maximum matches

  // check for second primer match

}


/* int readFile()
 * Parses the input file. Produces the output file(s).
 */
int readFile(FILE* out, FILE* gen, Primer* head,
    int minLen, int maxLen, float minScore) {

  char* chunk = (char*) memalloc(1 + CHUNK_SIZE);
  char* chrom = (char*) memalloc(1 + MAX_PRIM);
  Match* dummy = (Match*) memalloc(sizeof(Match));

  getLine(line, MAX_SIZE, gen);
  while (line[0] == '>') {

    // save chromosome
    int i;
    for (i = 0; line[i + 1] != '\n' &&
        line[i + 1] != ' ' && i < MAX_PRIM; i++)
      chrom[i] = line[i + 1];
    chrom[i] = '\0';

    chunk[0] = '\0';     // reset chunk
    dummy = NULL;        // reset list of matches for next chunk

    int count = 0;
    int pos = 0;
    int next = 1;
    while (next) {

      // load next chunk of genome
      next = getChunk(chunk, CHUNK_SIZE - MAX_PRIM, gen);

      for (Primer* p = head; p != NULL; p = p->next)
        findMatch(p, chunk, pos, dummy, minLen, maxLen, minScore);

//printf("chunk is %s\n%s\n", chunk, chrom);
//free(chrom);
//while (!getchar()) ;

//printf("chunk is\n%s\n", chunk);
//while (!getchar()) ;
      pos += CHUNK_SIZE - MAX_PRIM;
      count++;
    }

printf("chromosome %s, chunks of %dbp: %d\n", chrom, CHUNK_SIZE, count);
printf("last chunk is %s (%dbp)\n", chunk, (int) strlen(chunk));
while (!getchar()) ;

  }

  // free memory
  free(chunk);
  free(chrom);

  return 0;
}

/* int calcMax()
 * Calculates the maximum primer-genome score.
 */
int calcMax(char* prim) {
  int match = 0;
  int len = strlen(prim);
  for (int i = len - 1; i > -1; i--) {
    int val = 21 - len + i;
    if (val > 19)
      val *= 5;
    else if (val > 10)
      val *= 3;
    else if (val > 0)
      val *= 2;
    else
      val = 1;
    match += val;
  }
  return match;
}

/* Primer* loadSeqs(FILE*)
 * Loads the primers from the given file.
 */
Primer* loadSeqs(FILE* prim) {

  Primer* head = NULL, *prev = NULL;
  while (fgets(line, MAX_SIZE, prim) != NULL) {

    if (line[0] == '#')
      continue;

    // load name and sequence
    char* name = strtok(line, CSV);
    char* fwd = strtok(NULL, CSV);
    char* rev = strtok(NULL, DEL);
    if (name == NULL || fwd == NULL || rev == NULL) {
      error("", ERRPRIM);
      continue;
    }

    // check for duplicate
    for (Primer* pc = head; pc != NULL; pc = pc->next)
      if (!strcmp(pc->name, name))
        exit(error(name, ERRPREP));

    // create primer
    Primer* p = (Primer*) memalloc(sizeof(Primer));
    p->name = (char*) memalloc(1 + strlen(name));
    p->seq[0] = (char*) memalloc(1 + strlen(fwd));  // fwd primer
    p->seq[1] = (char*) memalloc(1 + strlen(fwd));  // rev-comp of fwd primer
    p->seq[2] = (char*) memalloc(1 + strlen(rev));  // rev primer
    p->seq[3] = (char*) memalloc(1 + strlen(rev));  // rev-comp of rev primer
    strcpy(p->name, name);
    strcpy(p->seq[0], fwd);
    strcpy(p->seq[2], rev);

    // save sequence rc's
    revComp(p->seq[0], p->seq[1]);
    revComp(p->seq[2], p->seq[3]);

    // calculate max. primer match scores
    p->fmax = calcMax(p->seq[0]);
    p->rmax = calcMax(p->seq[2]);

    p->next = NULL;
    if (head == NULL)
      head = p;
    else
      prev->next = p;
    prev = p;
  }

  closeFile(prim);
  return head;
}



/* void openFiles()
 * Opens the files to run the program.
 */
void openFiles(char* outFile, FILE** out,
    char* primFile, FILE** prim, char* genFile, FILE** gen,
    char* logFile, FILE** log, char* doveFile, FILE** dove,
    int dovetail) {
  // open required files
  *out = openFile(outFile, WRITE);
  *prim = openFile(primFile, READ);
  *gen = openFile(genFile, READ);

  // open optional files
  if (logFile != NULL) {
    *log = openFile(logFile, WRITE);
    fprintf(*log, "Read\tOverlapLen\tStitchedLen\tMismatch\n");
  }
  if (dovetail && doveFile != NULL) {
    *dove = openFile(doveFile, WRITE);
    fprintf(*dove, "Read\tDovetailFwd\tDovetailRev\n");
  }
}

/* void getParams()
 * Parses the command line.
 */
void getParams(int argc, char** argv) {

  char* outFile = NULL, *primFile = NULL, *genFile = NULL,
    *logFile = NULL,
    *doveFile = NULL;
  int minLen = DEFMIN, maxLen = DEFMAX;
  float minScore = DEFSCORE;
  int verbose = 0;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (!strcmp(argv[i], MAXOPT))
      maxLen = 0;
    else if (!strcmp(argv[i], VERBOSE))
      verbose = 1;
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], PRIMFILE))
        primFile = argv[++i];
      else if (!strcmp(argv[i], GENFILE))
        genFile = argv[++i];
      else if (!strcmp(argv[i], MINLEN))
        minLen = getInt(argv[++i]);
      else if (!strcmp(argv[i], MAXLEN))
        maxLen = getInt(argv[++i]);

      else if (!strcmp(argv[i], LOGFILE))
        logFile = argv[++i];
      else if (!strcmp(argv[i], DOVEFILE))
        doveFile = argv[++i];
      else if (!strcmp(argv[i], MINSCORE))
        minScore = getFloat(argv[++i]);
      else
        exit(error(argv[i], ERRPARAM));
    } else
      usage();
  }

  // check for parameter errors
  if (outFile == NULL || primFile == NULL || genFile == NULL)
    usage();
  if (minLen > maxLen)
    exit(error(LENERR, SPECERR));
  if (minScore <= 0 || minScore > 1)
    exit(error(SCOREERR, SPECERR));

  // open files
  FILE* out = NULL, *prim = NULL, *gen = NULL,
    *log = NULL, *dove = NULL;
int dovetail = 0;
  openFiles(outFile, &out, primFile, &prim, genFile, &gen,
    logFile, &log,
    doveFile, &dove, dovetail);
  Primer* head = loadSeqs(prim);
/*
  for (Primer* p = primo; p != NULL; p = p->next) {
    printf("%s\t%s\t%s\n\t%s\t%s\n", p->name, p->fwd, p->rev, p->frc, p->rrc);
    while (!getchar()) ;
  }

printf("Primer %s, score %d\nPrimer %s, score %d\n", p->fwd, p->fmax, p->rev, p->rmax);
while (!getchar()) ;
*/


  // read file
  int stitch = 0, fail = 0;  // counting variables
  int count = readFile(out, gen, head, minLen, maxLen, minScore);

  if (verbose) {
    printf("Reads analyzed: %d\n", count);
    printf("  Successfully stitched: %d\n", stitch);
    printf("  Stitch failures: %d\n", fail);
  }

  // close files
  closeFile(out);
  closeFile(gen);
  if (log != NULL);
    closeFile(log);
  if (dovetail && doveFile != NULL)
    closeFile(dove);

  freeMemory(head);
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  line = (char*) memalloc(MAX_SIZE);
  getParams(argc, argv);
  free(line);
  return 0;
}
