/*
  John M. Gaspar

  Various utilities and wrapper functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jmg_utils.h"

/* int error()
 * Prints an error message.
 */
int error(char* msg, int err) {
  char* msg2 = "";
  if (err == ERROPEN) msg2 = MERROPEN;
  else if (err == ERROPENW) msg2 = MERROPENW;
  else if (err == ERRCLOSE) msg2 = MERRCLOSE;
  else if (err == ERRMEM) msg2 = MERRMEM;
  else if (err == ERRINT) msg2 = MERRINT;
  else if (err == ERRFLOAT) msg2 = MERRFLOAT;
  else if (err == ERRLINE) msg2 = MERRLINE;
  else if (err == ERRUNK) msg2 = MERRUNK;

  else if (err == ERRQUAL) msg2 = MERRQUAL;
  else if (err == ERRPRIM) msg2 = MERRPRIM;
  else if (err == ERRPREP) msg2 = MERRPREP;

  else if (err == ERRHEAD) msg2 = MERRHEAD;
  else if (err == ERRPARAM) msg2 = MERRPARAM;
  else if (err == ERROVER) msg2 = MERROVER;
  else if (err == ERRMISM) msg2 = MERRMISM;
  else if (err != SPECERR) msg2 = DEFERR;

  fprintf(stderr, "Error! %s%s\n", msg, msg2);
  return -1;
}

/* FILE* openFile()
 * Opens a file for reading/writing.
 */
FILE* openFile(char* fileName, char* mode) {
  FILE* file = fopen(fileName, mode);
  if (file == NULL)
    exit(error(fileName,
      (mode[0] == READ[0] ? ERROPEN : ERROPENW)));
  return file;
}

/* void closeFile()
 * Closes a file.
 */
void closeFile(FILE* file) {
  if (fclose(file))
    exit(error("", ERRCLOSE));
}

/* void* memalloc()
 * Allocates a heap block.
 */
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char** endptr = NULL;
  int ans = (int) strtol(in, endptr, 10);
  if (endptr != '\0')
    exit(error(in, ERRINT));
  return ans;
}

/* float getFloat(char*)
 * Converts the given char* to a float.
 */
float getFloat(char* in) {
  char** endptr = NULL;
  float ans = strtof(in, endptr);
  if (endptr != '\0')
    exit(error(in, ERRFLOAT));
  return ans;
}

/* void getLine()
 * Reads the next line from the given FILE*.
 */
void getLine(char* line, int size, FILE* in) {
  if (fgets(line, size, in) == NULL)
    exit(error("", ERRLINE));
}


/* char comp(char)
 * Returns the complement of the given base.
 */
char comp(char in) {
  char out = '\0';
  if (in == 'A') out = 'T';
  else if (in == 'T') out = 'A';
  else if (in == 'C') out = 'G';
  else if (in == 'G') out = 'C';
  else if (in == 'Y') out = 'R';
  else if (in == 'R') out = 'Y';
  else if (in == 'W') out = 'W';
  else if (in == 'S') out = 'S';
  else if (in == 'K') out = 'M';
  else if (in == 'M') out = 'K';
  else if (in == 'B') out = 'V';
  else if (in == 'V') out = 'B';
  else if (in == 'D') out = 'H';
  else if (in == 'H') out = 'D';
  else if (in == 'N') out = 'N';
  else exit(error("", ERRUNK));
  return out;
}

/* void revComp()
 * Reverse-complements the given sequence.
 */
void revComp(char* in, char* out) {
  int j = 0;
  for (int i = strlen(in) - 1; i > -1; i--)
    out[j++] = comp(in[i]);
  out[j] = '\0';
}

/* int ambig(char, char)
 * Checks for matches of ambiguous DNA bases.
 */
int ambig(char x, char y) {
  if (x == 'N' ||
      (x == 'W' && (y == 'A' || y == 'T')) ||
      (x == 'S' && (y == 'C' || y == 'G')) ||
      (x == 'M' && (y == 'A' || y == 'C')) ||
      (x == 'K' && (y == 'G' || y == 'T')) ||
      (x == 'R' && (y == 'A' || y == 'G')) ||
      (x == 'Y' && (y == 'C' || y == 'T')) ||
      (x == 'B' && (y == 'C' || y == 'G' || y == 'T')) ||
      (x == 'D' && (y == 'A' || y == 'G' || y == 'T')) ||
      (x == 'H' && (y == 'A' || y == 'C' || y == 'T')) ||
      (x == 'V' && (y == 'A' || y == 'C' || y == 'G')))
    return 0;
  return 1;
}

