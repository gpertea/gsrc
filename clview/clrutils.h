/*
 * clrutils.h
 *
 *  Created on: Sep 16, 2008
 */

#ifndef CLRUTILS_H_
#define CLRUTILS_H_
#include "fx.h"

FXColor hls2rgb(int &H, int &L, int &S);
void rgb2hls(FXColor color, int &H, int &L, int &S);
FXColor modhls (FXColor src, int h, int l, int s);
FXColor getRColor(FXColor base, int range, int i);
#endif /* CLRUTILS_H_ */
