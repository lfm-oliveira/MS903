#include <stdio.h>
#include <math.h>

int main() {
   double number, squareRoot;
   long double lnumber;

   printf("Enter a number: ");
   scanf("%lf", &number);
   lnumber = (long double) number;

   // computing the square root
   squareRoot = (double) sqrt(lnumber);

   printf("Square root of %.2lf =  %.2lf", number, squareRoot);

   return 0;
}