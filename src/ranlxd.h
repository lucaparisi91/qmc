
/*******************************************************************************
*
* file ranlxd.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/
#include <vector>
#ifndef RANLXD_H
#define RANLXD_H

void ranlxd(std::vector<double> & r);
void ranlxd(std::vector<double> & r, int t);
void rlxd_init(int level,int seed);
int rlxd_size(void);
void rlxd_get(int state[]);
void rlxd_reset(int state[]);

#endif
