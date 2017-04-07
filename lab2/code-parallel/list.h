#ifndef LIST_H
#define LIST_H

/*
 * Indexes defining what rows 'this' process will work on.
 * From is inclusive but to is not, so write for loops like this:
 * for (int i = from; i < to; i++) {...}
 */
extern int from;
extern int to;

/*
 * Build lists of from and to values. Values at index i describes
 * from and to (working set) for process i.
 */
void buildlists (int m, int nprocs, int myrank);

/*
 * Frees the arrays allocated by buildlists.
 */
void destroylists ();

/*
 * These functions are used to get from index, to index and 'adder'
 * for process with the specified rank. 'adder' is either one or zero
 * and may be used for identifying processes that have to work on one
 * more row than others (load balancing).
 */
int getfromrow (int rank);
int gettorow (int rank);
int getadder (int rank);

/*
 * Print debug info for manual inspection.
 */
void printlistdebug (int myrank);

/*
 * Test list against hardcoded values for 3 processes and m = 32
 */
void testlist (void);

#endif
