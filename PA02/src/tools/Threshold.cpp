#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "image.h"

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);

int main(int argc, char *argv[])
{
 int i, j; 
 int M, N, Q;
 bool type;
 int val;
 int thresh;

 // read image header
 readImageHeader(argv[1], N, M, Q, type);

 // allocate memory for the image array

 ImageType image(N, M, Q);

 // read image
 readImage(argv[1], image);

 cout << "Enter threshold: ";
 cin >> thresh;

 // threshold image 

 for(i=0; i<N; i++)
   for(j=0; j<M; j++) {
     image.getPixelVal(i, j, val);
     if(val < thresh) 
       image.setPixelVal(i, j, 255);
     else
       image.setPixelVal(i, j, 0);
    }

 // write image
 writeImage(argv[2], image);

 return (1);
}
