#pragma once

/** The number of distinguishable base pairs */
#define NBPAIRS 7

const int BP_pair[5][5]=
/* @  A  C  G  U*/
{{ 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5},
 { 0, 0, 0, 1, 0},
 { 0, 0, 2, 0, 3},
 { 0, 6, 0, 4, 0}};

/* rtype[pair[i][j]]:=pair[j][i] */
const int rtype[7] = {0, 2, 1, 4, 3, 6, 5};


// constants
/** Infinity as used in minimization routines */
#define INF 10000000
