/*
 * clrutils.cpp
 *
 *  Created on: Sep 16, 2008
 *      Author: gpertea
 */

/*================== useful rgb <-> hls utility functions I should put
 in a separate module
*/
#include "clrutils.h"
#include "GBase.h"

inline double h_value(double n1, double n2, double hue) {
     if (hue<0)   hue += 360;
     if (hue>360) hue -= 360;
     if (hue<60) return n1+(n2-n1)*hue/60.;
            else if (hue<180) return n2;
                   else if (hue<240)
                           return n1+(n2-n1)*(240-hue)/60.;
                         else
                           return n1;
}

FXColor hls2rgb(int &H, int &L, int &S) {
   double h,l,s;
   double r,g,b,m1,m2;
   h = ((double)H)/255. * 360;
   l = ((double)L)/255.;
   s = ((double)S)/255.;
   if (l<.5) m2 = l*(1+s);
      else m2 = l+s-l*s;
   m1 = 2*l-m2;
   if (s==0) {
        r=l;
        g=l;
        b=l;
       }
    else {
       r=h_value(m1,m2,h+120);
       g=h_value(m1,m2,h);
       b=h_value(m1,m2,h-120);
       }
   return FXRGB((int)(r*255), (int)(g*255), (int)(b*255));
}

void rgb2hls(FXColor color, int &H, int &L, int &S) {
  double r,g,b;
  double h,l,s,min,max,delta;
  r = ((double)FXREDVAL(color)) / 255. ;
  g = ((double)FXGREENVAL(color)) / 255. ;
  b = ((double)FXBLUEVAL(color)) / 255. ;
  min = FXMIN3(r,g,b);
  max = FXMAX3(r,g,b);
  l = (max+min)/2;
  if (max==min) {
          s = 0;
          h = 0;
          }
     else {
       if (l<=0.5) s = (max-min)/(max+min);
            else   s = (max-min)/(2-max-min);
       delta = max-min;
       if (r==max) h=(g-b)/delta;
           else if (g==max) h=2+(b-r)/delta;
                   else h=4+(r-g)/delta;
       h /= 6.;
       if (h<0) h ++;
       }
  H = iround(h*255.0);
  L = iround(l*255.0);
  S = iround(s*255.0);
}

//adjust h,l,s values for a color
FXColor modhls (FXColor src, int h, int l, int s) {
 int H,L,S;
 rgb2hls(src, H, L, S);
 H+=h;L+=l;S+=s;
 if (H>255) H %= 255;
 //H=FXCLAMP(0,H,255);
 L=FXCLAMP(0,L,255);
 S=FXCLAMP(0,S,255);
 return hls2rgb(H,L,S);
}

FXColor getRColor(FXColor base, int range, int i) {
  int hstep=255/(range>>2);
  if (hstep>70) hstep=70;
  if (hstep<10) hstep=35;
  return modhls(base,hstep*(i+1), (i/(range>>2))*10, (i/(range>>2)*12));
}

