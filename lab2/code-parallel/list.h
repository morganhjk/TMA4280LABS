#ifndef LIST_H
#define LIST_H

extern int from;
extern int to;

void definearea (int m);
void buildlists (int m, int nprocs);
void destroylists ();
int getfromrow (int rank);
int gettorow (int rank);

#endif
