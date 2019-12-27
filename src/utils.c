#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

//function to shuffle the contents of an array
void shuffle(int *array, int n)
{
    if(n < 1){
      printf("Invalid input %d", n);
      return;
    }

    for (int i = 0; i < n - 1; i++)
    {
      int j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }

}

//function to sort a given integer array of size n
void bubbleSort(int ar[], int n)
{
    for (int i = 0; i < n-1; i++)
      for (int j = 0; j < n-i-1; j++)
          if (ar[j] > ar[j+1]){
            int t = ar[j];
            ar[j] = ar[j+1];
            ar[j+1] = t;
          }
}

//function to multiply making sure overflow doesn't occur
int multiplyOverflow(int a, int b, int* result)
{
    if (a == 0 || b == 0)
        return 0;

    *result = a * b;
    if (a == (*result) / b)
        return 0;
    else
        return -1;
}

//function to add making sure overflow doesn't occur
int additionOverflow(int a, int b, int* result)
 {
     *result = a + b;
     if(a > 0 && b > 0 && *result < 0)
         return -1;
     if(a < 0 && b < 0 && *result > 0)
         return -1;
     return 0;
 }
