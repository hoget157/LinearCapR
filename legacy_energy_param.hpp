/*
extracted from CapR
*/

#pragma once

#ifndef ENERGY_PARAM_NAMESPACE_BEGIN
#define ENERGY_PARAM_NAMESPACE_BEGIN
#define ENERGY_PARAM_NAMESPACE_END
#define ENERGY_PARAM_NAMESPACE_INTERNAL_GUARD
#endif

ENERGY_PARAM_NAMESPACE_BEGIN

#include "miscs.hpp"

// constants

/** The gas constant */
#define GASCONST 1.98717  /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0  273.15
/** The minimum loop length */
#define TURN 3
/** The maximum loop length */
#define MAXLOOP 30

#define MULTI_MAX_UNPAIRED 30

const double temperature = 37;
const double kT = (temperature + K0) * GASCONST / 10.;

const double lxc37 = 107.856;
const int ML_intern37 = 40;
const int ML_closing37 = 340;
const int ML_BASE37 = 0;
const int MAX_NINIO = 300;
const int ninio37 = 50;
const int TerminalAU37 = 50;

// stack
const int stack37[NBPAIRS+1][NBPAIRS+1] =
/*          CG     GC     GU     UG     AU     UA     NN*/
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF,  -240,  -330,  -210,  -140,  -210,  -210,  -140}
,{   INF,  -330,  -340,  -250,  -150,  -220,  -240,  -150}
,{   INF,  -210,  -250,   130,   -50,  -140,  -130,   130}
,{   INF,  -140,  -150,   -50,    30,   -60,  -100,    30}
,{   INF,  -210,  -220,  -140,   -60,  -110,   -90,   -60}
,{   INF,  -210,  -240,  -130,  -100,   -90,  -130,   -90}
,{   INF,  -140,  -150,   130,    30,   -60,   -90,   130}};

// initiation
const int hairpin37[31] = {
  INF, INF, INF, 570, 560, 560, 540, 590, 560, 640, 650,
       660, 670, 678, 686, 694, 701, 707, 713, 719, 725,
       730, 735, 740, 744, 749, 753, 757, 761, 765, 769};
const int bulge37[31] = {
  INF, 380, 280, 320, 360, 400, 440, 459, 470, 480, 490,
       500, 510, 519, 527, 534, 541, 548, 554, 560, 565,
	   571, 576, 580, 585, 589, 594, 598, 602, 605, 609};
const int internal_loop37[31] = {
  INF, INF, 410, 510, 170, 180, 200, 220, 230, 240, 250,
       260, 270, 278, 286, 294, 301, 307, 313, 319, 325,
       330, 335, 340, 345, 349, 353, 357, 361, 365, 369};

const int mismatchI37[NBPAIRS][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0, -110,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0, -110,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -70}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0, -110,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0, -110,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -70}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
};

inline const int (*const mismatch1nI37)[5][5] = mismatchI37;
inline const int (*const mismatch23I37)[5][5] = mismatchI37;

 const int mismatchH37[NBPAIRS][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { -90, -150, -150, -140, -180}, /* A@  AA  AC  AG  AU */
   { -90, -100,  -90, -290,  -80}, /* C@  CA  CC  CG  CU */
   { -90, -220, -200, -160, -110}, /* G@  GA  GC  GG  GU */
   { -90, -170, -140, -180, -200}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { -70, -110, -150, -130, -210}, /* A@  AA  AC  AG  AU */
   { -70, -110,  -70, -240,  -50}, /* C@  CA  CC  CG  CU */
   { -70, -240, -290, -140, -120}, /* G@  GA  GC  GG  GU */
   { -70, -190, -100, -220, -150}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   20,  -50,  -30,  -30}, /* A@  AA  AC  AG  AU */
   {   0,  -10,  -20, -150,  -20}, /* C@  CA  CC  CG  CU */
   {   0,  -90, -110,  -30,    0}, /* G@  GA  GC  GG  GU */
   {   0,  -30,  -30,  -40, -110}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,  -50,  -30,  -60,  -50}, /* A@  AA  AC  AG  AU */
   {   0,  -20,  -10, -170,    0}, /* C@  CA  CC  CG  CU */
   {   0,  -80, -120,  -30,  -70}, /* G@  GA  GC  GG  GU */
   {   0,  -60,  -10,  -60,  -80}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,  -30,  -50,  -30,  -30}, /* A@  AA  AC  AG  AU */
   {   0,  -10,  -20, -150,  -20}, /* C@  CA  CC  CG  CU */
   {   0, -110, -120,  -20,   20}, /* G@  GA  GC  GG  GU */
   {   0,  -30,  -30,  -60, -110}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,  -50,  -30,  -60,  -50}, /* A@  AA  AC  AG  AU */
   {   0,  -20,  -10, -120,   -0}, /* C@  CA  CC  CG  CU */
   {   0, -140, -120,  -70,  -20}, /* G@  GA  GC  GG  GU */
   {   0,  -30,  -10,  -50,  -80}}/* U@  UA  UC  UG  UU */
};


 /* dangle5 */
const int dangle5_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -10,   -50,   -30,   -20,   -10},
/* GC */ {    -0,   -20,   -30,    -0,    -0},
/* GU */ {   -20,   -30,   -30,   -40,   -20},
/* UG */ {   -10,   -30,   -10,   -20,   -20},
/* AU */ {   -20,   -30,   -30,   -40,   -20},
/* UA */ {   -10,   -30,   -10,   -20,   -20},
/* NN */ {    -0,   -20,   -10,    -0,    -0}
};

/* dangle3 */
const int dangle3_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -40,  -110,   -40,  -130,   -60},
/* GC */ {   -80,  -170,   -80,  -170,  -120},
/* GU */ {   -10,   -70,   -10,   -70,   -10},
/* UG */ {   -50,   -80,   -50,   -80,   -60},
/* AU */ {   -10,   -70,   -10,   -70,   -10},
/* UA */ {   -50,   -80,   -50,   -80,   -60},
/* NN */ {   -10,   -70,   -10,   -70,   -10}
};

#include "legacy_intloops.hpp"

ENERGY_PARAM_NAMESPACE_END

#ifdef ENERGY_PARAM_NAMESPACE_INTERNAL_GUARD
#undef ENERGY_PARAM_NAMESPACE_BEGIN
#undef ENERGY_PARAM_NAMESPACE_END
#undef ENERGY_PARAM_NAMESPACE_INTERNAL_GUARD
#endif
