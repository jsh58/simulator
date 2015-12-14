/*
  John Gaspar
  Dec. 2015

  Header file for PCRSim.c.
*/

#define MAX_SIZE    1024    // maximum length for input line
#define CHUNK_SIZE  50     // maximum chunk of genome to analyze
//#define CHUNK_SIZE  65536   // maximum chunk of genome to analyze
#define MAX_PRIM    10      // maximum primer length
#define CSV         ",\t"   // delimiter for primer file
#define DEL         ",\t\n"

// command-line parameters
#define HELP        "-h"
#define PRIMFILE    "-p"
#define GENFILE     "-f"    // fasta genome
#define OUTFILE     "-o"
#define MINLEN      "-m"
#define MAXLEN      "-M"
#define MINSCORE    "-s"

#define LOGFILE     "-l"
#define DOVEOPT     "-d"
#define DOVEFILE    "-dl"
#define MAXOPT      "-n"
#define VERBOSE     "-ve"

// default parameter values
#define DEFMIN      60     // minimum amplicon length
#define DEFMAX      300    // maximum amplicon length
#define DEFSCORE    0.75f  // primer-genome match score

// third parameter to copyStr()
#define FWD         0
#define RC          1
#define REV         2

// custom error messages
#define LENERR      "Min. amplicon length cannot be larger than max."
#define SCOREERR    "Min. score must be in (0,1]"

// structs
typedef struct match {
  float fmatch;
  float rmatch;
  int fpos;
  int rpos;
  int chrom;
  struct match* next;
} Match;

typedef struct primer {
  char* name;
  char* fwd;
  char* rev;
  char* frc;
  char* rrc;
  int fmax;    // max. match score for fwd primer
  int rmax;    // max. match score for rev primer
  Match* first;
  struct primer* next;
} Primer;
