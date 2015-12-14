/*
  John M. Gaspar

  Header file for jmg_utils.c.
*/

// functions
int error(char*, int);            // prints an error message
FILE* openFile(char*, char*);     // wrapper for fopen()
void closeFile(FILE*);            // wrapper for fclose()
void* memalloc(int);              // wrapper for malloc()
float getFloat(char*);            // wrapper for strtof()
int getInt(char*);                // wrapper for strtol()
void getLine(char*, int, FILE*);  // wrapper for fgets()
char comp(char);                  // produces the complementary nucleotide
void revComp(char*, char*);       // reverse-complements a sequence
int ambig(char, char);            // checks for matches of ambiguous nucleotides


#define READ       "r"
#define WRITE      "w"

// error messages
#define ERROPEN     0
#define MERROPEN    ": cannot open file for reading"
#define ERROPENW    1
#define MERROPENW   ": cannot open file for writing"
#define ERRCLOSE    2
#define MERRCLOSE   "Cannot close file"
#define ERRMEM      3
#define MERRMEM     "Cannot allocate memory"
#define ERRINT      4
#define MERRINT     ": cannot convert to int"
#define ERRFLOAT    5
#define MERRFLOAT   ": cannot convert to float"
#define ERRLINE     6
#define MERRLINE    "Cannot read line from file"


#define ERRUNK      7
#define MERRUNK     "Unknown nucleotide"



#define ERRPRIM     6
#define MERRPRIM    "cannot load primer sequence"
#define ERRPREP     7
#define MERRPREP    ": cannot repeat primer name"


#define ERRQUAL     6
#define MERRQUAL    "Sequence/quality scores do not match"
#define ERRHEAD     7
#define MERRHEAD    ": not matched in input files"
#define ERRPARAM    10
#define MERRPARAM   ": unknown command-line parameter"
#define ERROVER     11
#define MERROVER    "Overlap must be greater than 0"
#define ERRMISM     12
#define MERRMISM    "Mismatch must be in [0,1)"
#define SPECERR     -1
#define DEFERR      "Unknown error"
