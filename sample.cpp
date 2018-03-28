#ifndef __PROGTEST__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <algorithm>
#include <pthread.h>
#include <semaphore.h>
#include <cstdint>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <memory>
#include <condition_variable>
#include <atomic>
using namespace std;


class CFITCoin;
class CCVUTCoin;
class CCustomer;

typedef struct shared_ptr<CFITCoin>                        AFITCoin;
typedef struct shared_ptr<CCVUTCoin>                       ACVUTCoin;
typedef struct shared_ptr<CCustomer>                       ACustomer;

//=================================================================================================

class CFITCoin
{
  public:
                             CFITCoin                      ( const vector<uint32_t> & vectors,
                                                             int               distMax )
                             : m_Vectors ( vectors ),
                               m_DistMax ( distMax ),
                               m_Count ( 0 )
    {
    }
    virtual                  ~CFITCoin                     ( void ) = default;
    vector<uint32_t>         m_Vectors;
    int                      m_DistMax;
    uint64_t                 m_Count;
};

//=================================================================================================

class CCVUTCoin       
{
  public:
                             CCVUTCoin                     ( const vector<uint8_t> & data,
                                                             int               distMin,
                                                             int               distMax )
                             : m_Data ( data ),
                               m_DistMin ( distMin ),
                               m_DistMax ( distMax ),
                               m_Count ( 0 )
    {
    }
    virtual                  ~CCVUTCoin                    ( void ) = default;
    vector<uint8_t>          m_Data;
    int                      m_DistMin;
    int                      m_DistMax;
    uint64_t                 m_Count;
};

//=================================================================================================

class CCustomer
{
  public:
    virtual                  ~CCustomer                    ( void ) = default;
    virtual AFITCoin         FITCoinGen                    ( void ) = 0;
    virtual ACVUTCoin        CVUTCoinGen                   ( void ) = 0;
  
    virtual void             FITCoinAccept                 ( AFITCoin          x ) = 0;
    virtual void             CVUTCoinAccept                ( ACVUTCoin         x ) = 0;
};

//=================================================================================================
#endif /* __PROGTEST__ */


class CRig
{
  public:
    static void Solve (ACVUTCoin x);
    static void Solve (AFITCoin x);
    CRig (void);
    ~CRig (void);
    void Start (int thrCnt);
    void Stop (void);
    void AddCustomer (ACustomer c);
  private:
    static void bin(unsigned n);
    static void printBinary(int i, int j, int diff, int mask);
    static unsigned int countSetBits(int n);
    static int binomialCoeff(int n, int k);
    static int v;
    static int f;
    static void prepareVectors(AFITCoin &x);
    vector<uint32_t> shortVectors;
};

//via geeksforgeeks
int CRig::binomialCoeff(int n, int k)
{
  // Base Cases
  if (k==0 || k==n)
    return 1;

  // Recur
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k);
}

//via geeksforgeeks
void CRig::bin(unsigned n)
{
    unsigned i;
    for (i = 1 << 31; i > 0; i = i / 2)
        (n & i)? printf("1"): printf("0");
}

void CRig::printBinary(int i, int j, int diff, int mask){
    bin(i);
    printf(" (%d) XOR\n",i);
    bin(j);
    printf(" is:\n");
    bin(mask);
    printf(" = %d\n", diff);
    printf("-----------------------\n");
}

//via geeksforgeeks
unsigned int CRig::countSetBits(int n)
{
    unsigned int count = 0;
    while (n)
    {
      n &= (n-1) ;
      count++;
    }
    return count;
}

void CRig::prepareVectors(AFITCoin &x){
    vector<int> indexes;
    int mask;
    uint32_t testedBit;

    for (uint32_t i = 0; i < 32; i++) { //through all bits
        mask = 1 << i;
        printf("___________________________________\n");
        printf("i = %d Mask is: \n", i);
        bin(mask);
        printf("\n");

        for(uint32_t j = 0; j != x->m_Vectors.size(); j++) { //for each vector
            testedBit = x->m_Vectors[j] & mask;
            bin(testedBit);
            printf("\n");
        }
    }
}

void CRig::Solve (ACVUTCoin x) {

}

void CRig::Solve (AFITCoin x) {
    uint32_t mask;
    int diff, maxDiff;

    prepareVectors(x);

    for (uint32_t i = 0; i < 16; i++){ //through all numbers 0-2^32
        maxDiff = 0;
        for(uint32_t j = 0; j != x->m_Vectors.size(); j++) { //through all vectors given
            mask = i ^ x->m_Vectors[j]; // xor operation to find different bits between number i and vector j
            diff = countSetBits(mask);

            if(diff > maxDiff)
                maxDiff = diff;
        }

        if (maxDiff <= x->m_DistMax){
            x->m_Count++;
        }
    }
}


CRig::CRig (void) {
}

CRig::~CRig (void) {

}

void CRig::Start (int thrCnt) {

}

void CRig::Stop (void) {

}

void CRig::AddCustomer (ACustomer c) {

}




#ifndef __PROGTEST__
#include "test.cpp"
#endif /* __PROGTEST__ */
