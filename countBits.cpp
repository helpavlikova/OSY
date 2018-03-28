#include <stdio.h>
#include <bitset>
#include <iostream>

/* Function to get no of set bits in binary
   representation of positive integer n */
unsigned int countSetBits(unsigned int n)
{
  unsigned int count = 0;
  while (n)
  {
    count += n & 1;
    std::cout << std::bitset<4>(n) << std::endl;
    n >>= 1;
  }
  return count;
}

/* Program to test function countSetBits */
int main()
{
    int i = 9;
    printf("%d\n", countSetBits(i));
    return 0;
}
