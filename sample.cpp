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
    static int varCount;
    static int testCounter;
    static uint64_t fixCount;
    static int bitsLen;
    static vector<uint32_t> shortVectors;
    static void prepareVectors(AFITCoin &x);
    static void printVectors(vector<int>& vectors);
    static void printVectors(vector<uint32_t>& vectors);
    static int swapBits(unsigned int n, unsigned int p1, unsigned int p2);
    static uint64_t combinationSum(uint64_t k);
};


int CRig::varCount;
uint64_t CRig::fixCount;
int CRig::testCounter = 0;
int CRig::bitsLen = 32;
vector<uint32_t> CRig::shortVectors;

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

void CRig::printBinary(int i, int j, int diff, int mask) {
    bin(i);
    printf(" (%d) XOR\n",i);
    bin(j);
    printf(" is:\n");
    bin(mask);
    printf(" = %d\n", diff);
    printf("-----------------------\n");
}

uint64_t CRig::combinationSum(uint64_t k) { //kombinatoricka magie z prednasky
    uint64_t result = 0;
    k = min(k, fixCount); //mensi z cisel k ci f

    for (uint64_t i = 0; i <= k; i++) { //suma az do k vcetne
        result += binomialCoeff(fixCount, i); //f nad i
        printf("combinatorial magic: sum from 0 to %zu C(%zu, %zu) = %zu\n", k, fixCount, i, result);
    }
    printf("combinatorial magic: sum from 0 to %zu C(%zu, i) = %zu\n", k, fixCount, result);
    return result;
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

//via geeksforgeeks
int CRig::swapBits(unsigned int n, unsigned int p1, unsigned int p2)
{
    /* Move p1'th to rightmost side */
    unsigned int bit1 =  (n >> p1) & 1;

    /* Move p2'th to rightmost side */
    unsigned int bit2 =  (n >> p2) & 1;

    /* XOR the two bits */
    unsigned int x = (bit1 ^ bit2);

    /* Put the xor bit back to their original positions */
    x = (x << p1) | (x << p2);

    /* XOR 'x' with the original number so that the
       two sets are swapped */
    unsigned int result = n ^ x;

    return result;
}

void CRig::printVectors(vector<int>& vectors) {
    for (unsigned int i = 0; i < vectors.size(); i++) {
        printf(" %d ", vectors[i]);
    }
    printf("\n");
}

void CRig::printVectors(vector<uint32_t>& vectors) {
    for (unsigned int i = 0; i < vectors.size(); i++) {
        bin(vectors[i]);
        printf("\n");
    }
    printf("\n");
}

void CRig::prepareVectors(AFITCoin &x) {
    vector<int> indexes;
    int mask;
    uint32_t testedBit1;
    uint32_t testedBit2;

    printf("original vectors:\n");
    printVectors(x->m_Vectors);

    //find out indexes of varying bits
    for (int i = 0; i < bitsLen; i++) { //through all bits
        mask = 1 << i;
        for(uint32_t j = 0; j != x->m_Vectors.size() - 1; j++) { //for each vector
            testedBit1 = x->m_Vectors[j] & mask;
            testedBit2 = x->m_Vectors[j+1] & mask;
            if ((testedBit1 ^ testedBit2) >> i == 1){
                indexes.push_back(i);\
                break;
            }
        }
    }

    varCount = indexes.size();
    printf("varCount: %d\n", varCount);
    fixCount = bitsLen - varCount;
    printf("fixCount: %zu\n", fixCount);

    //swap varying with fixed bits
    for (uint32_t i = 0; i != x->m_Vectors.size(); i++ ) {
         uint32_t swapper = x->m_Vectors[i];
         for(uint32_t j = 0; j < indexes.size(); j++) {
            swapper = swapBits(swapper,j, indexes[j]);
         }

         //zero out the fixed bits
         swapper <<= fixCount;
         swapper >>= fixCount;

         shortVectors.push_back(swapper);
    }


    printf("shortened vectors:\n");
    printVectors(shortVectors);
    indexes.clear();
}

void CRig::Solve (AFITCoin x) {
    uint32_t mask;
    int diff = 0;
    int maxDiff;

    testCounter++;
    printf("test %d: \n", testCounter);



    prepareVectors(x);

    if (x->m_Vectors.size() == 1) {
        x->m_Count = combinationSum(x->m_DistMax);
        shortVectors.clear();
        return;
    }

    for (uint32_t i = 0; i < pow(2, varCount); i++) { //through all number from 0 to 2^v

        maxDiff = 0;

        for(uint32_t j = 0; j != shortVectors.size(); j++) { //through all vectors given

            mask = i ^ shortVectors[j]; // xor operation to find different bits between number i and vector j
            diff = countSetBits(mask);

            if(diff > maxDiff) {
                maxDiff = diff;
            }
        }

        if (maxDiff == x->m_DistMax) {
            x->m_Count++;
        }
        else if (maxDiff < x->m_DistMax) {
            x->m_Count += combinationSum(x->m_DistMax - maxDiff);
        }
    }

    shortVectors.clear();

}

void CRig::Solve (ACVUTCoin x) {

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
