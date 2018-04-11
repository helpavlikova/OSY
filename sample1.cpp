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

class CCoin {
public:
    CCoin(bool coinType, AFITCoin fitCoinRef, ACVUTCoin cvutCoinRef, int idx, int ID):
        isFit(coinType), fitCoin(fitCoinRef), cvutCoin(cvutCoinRef), customerIdx(idx), ID(ID) { }
    bool isFit; // if 1, then FitCoin, if 0 then CvutCoin
    AFITCoin fitCoin;
    ACVUTCoin cvutCoin;
    int customerIdx;
    int ID;
    void printCoin();
private:

};

void CCoin::printCoin() {
    if (isFit) {
        printf("FITCoin %d: distMax: %d (%zu)\n", ID, fitCoin->m_DistMax,  fitCoin->m_Count);
    } else {
        printf ("CVUTCoin %d: <distMin, distMax>: <%d,%d> (%zu)\n", ID, cvutCoin->m_DistMin, cvutCoin->m_DistMax, cvutCoin->m_Count);
    }
}

class CustomerWrapper {
public:
    CustomerWrapper(ACustomer newCustomer): customerRef(newCustomer) {
        sem_init ( &semCusEmpty, 0, 100 );
        sem_init ( &semCusFull, 0, 0 );
    }
    ACustomer customerRef;
    deque<CCoin> solvedCoins;
    sem_t semCusEmpty; //personal semaphore for each customer
    sem_t semCusFull;
private:
};

class CRig
{
  public:
    static void Solve (ACVUTCoin x);
    static void Solve (AFITCoin x);
    CRig (void);
    ~CRig (void) {}
    void Start (int thrCnt);
    void Stop (void);
    void AddCustomer (ACustomer c);
  private:
    static void bin(unsigned n, int bitsLen);
    static void printBinary(int i, int j, int diff, int mask);
    static void printVectors(vector<int>& vectors);
    static void printVectors(vector<bool>& vectors);
    static void printVectors(vector<uint32_t>& vectors);
    static void printVectors(vector<uint8_t>& vectors);
    void printCustomers();
    void printBuffer();
    static unsigned int countSetBits(int n);
    static uint64_t binomialCoeff(uint64_t n, uint64_t k);
    static void prepareData(ACVUTCoin &x, vector<bool> & boolVector);
    static void findPrefixSuffix(ACVUTCoin &x, vector<bool> & boolVector);
    static int editDistance(vector<bool>&  word1, vector<bool>& word2);
    static int swapBits(unsigned int n, unsigned int p1, unsigned int p2);
    static uint64_t combinationSum(uint64_t k, uint64_t fixCount);
    static void findVarBits(AFITCoin &x, uint64_t &fixCount, int& varCount,vector<int>& indexes);
    static void removeFixBits(AFITCoin &x, uint64_t &fixCount,vector<int>& indexes, vector<uint32_t>& shortVectors);
    static void compareWithNumbers(AFITCoin &x,int& varCount,uint64_t &fixCount, vector<uint32_t>& shortVectors);
    static int min(int x, int y, int z);
    static uint64_t min(uint64_t x, uint64_t y);
    void AddFitCoins(ACustomer c, int custIdx);
    void AddCvutCoins(ACustomer c,  int custIdx);
    void AcceptCoin(ACustomer c, int idx);
    void SolveCoin();
    static int FITtestCounter;
    static int CVUTtestCounter;
    static int bitsLen;
    static int customerIndex;
    deque<CCoin> coinBuffer;
    vector<thread> workThreads;
    vector<thread> customerThreads;
    vector<CustomerWrapper> customers;
    bool endFlag;
    int customerCounter;
    int workCount;
    mutex mtxCoinBuffer;
    mutex mtxSolved;
    mutex mtx;
    sem_t semCoinEmpty;
    sem_t semCoinFull;

};

//declaration of static variables
int CRig::FITtestCounter = 0;
int CRig::CVUTtestCounter = 0;
int CRig::bitsLen = 32;
int CRig::customerIndex = 0;

//--------------------------------------- Printing methods -----------------------------------------------

void CRig::printCustomers(){
    for (unsigned i = 0; i < customers.size(); i++) {
        printf ("customer %d:\n",i);
        for (unsigned j = 0; j < customers[i]. solvedCoins . size (); j ++) {
            customers[i]. solvedCoins[j].printCoin();
        }
    }
}

//via geeksforgeeks
void CRig::bin(unsigned n, int len) {
    unsigned i;
    for (i = 1 << (len - 1); i > 0; i = i / 2)
        (n & i)? printf("1"): printf("0");
}

void CRig::printBinary(int i, int j, int diff, int mask) {
    bin(i, bitsLen);
    printf(" (%d) XOR\n",i);
    bin(j, bitsLen);
    printf(" is:\n");
    bin(mask, bitsLen);
    printf(" = %d\n", diff);
    printf("-----------------------\n");
}

void CRig::printVectors(vector<int>& vectors) {
    for (unsigned int i = 0; i < vectors.size(); i++) {
        printf(" %d ", vectors[i]);
    }
    printf("\n");
}

void CRig::printVectors(vector<uint32_t>& vectors) {
    for (unsigned int i = 0; i < vectors.size(); i++) {
        bin(vectors[i], bitsLen);
        printf("\n");
    }
    printf("\n");
}

void CRig::printVectors(vector<uint8_t>& vectors) {
    for (unsigned int i = 0; i < vectors.size(); i++) {
        bin(vectors[i], 8);
        printf("\n");
    }
    printf("\n");
}

void CRig::printVectors(vector<bool>& vectors) {
    for (unsigned int i = 0; i < vectors.size(); i++) {
        printf("%d",  static_cast<int>(vectors[i]));
    }
    printf("\n");
}

void CRig::printBuffer() {
    int i = 0;
    for ( auto it = this->coinBuffer . begin (); it != this->coinBuffer . end (); it++ ) {
        printf ("[%d] %d ",i, it -> customerIdx );
        it -> printCoin();
        i++;
    }
}

//--------------------------------------- FITCoin methods -----------------------------------------------

uint64_t CRig::min(uint64_t x, uint64_t y) {
    if (x < y)
        return x;
    else
        return y;
}

//via geeksforgeeks
uint64_t CRig::binomialCoeff(uint64_t n, uint64_t k) {
    uint64_t res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k ) {
        k = n - k;
    }

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (uint64_t i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}

uint64_t CRig::combinationSum(uint64_t k, uint64_t fixCount) { //kombinatoricka magie z prednasky
    uint64_t result = 0;
    k = min(k, fixCount); //mensi z cisel k ci f

    for (uint64_t i = 0; i <= k; i++) { //suma az do k vcetne
        result += binomialCoeff(fixCount, i); //f nad i
    }
    return result;
}

//via geeksforgeeks
unsigned int CRig::countSetBits(int n) {
    unsigned int count = 0;
    while (n) {
      n &= (n-1) ;
      count++;
    }
    return count;
}

//via geeksforgeeks
int CRig::swapBits(unsigned int n, unsigned int p1, unsigned int p2) {
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

void CRig::findVarBits(AFITCoin &x, uint64_t &fixCount, int& varCount,vector<int>& indexes) {
    int mask;
    uint32_t testedBit1, testedBit2;

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
    fixCount = bitsLen - varCount;
}

void CRig::removeFixBits(AFITCoin &x, uint64_t &fixCount,vector<int>& indexes, vector<uint32_t>& shortVectors) {
    //swap varying with fixed bits
    for (uint32_t i = 0; i != x->m_Vectors.size(); i++ ) {
         uint32_t swapper = x->m_Vectors[i];
         for(uint32_t j = 0; j < indexes.size(); j++) {
            swapper = swapBits(swapper,j, indexes[j]);
         }

         //zero out the fixed bits
         swapper <<= fixCount;
         swapper >>= fixCount;

         shortVectors.push_back(swapper); //adding the number without fixed bits to the new vector
    }

}

void CRig::compareWithNumbers(AFITCoin &x, int& varCount,uint64_t &fixCount,vector<uint32_t>& shortVectors ) {
    uint32_t mask;
    int diff = 0;
    int maxDiff;

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
            x->m_Count += combinationSum(x->m_DistMax - maxDiff, fixCount);
        }
    }

    shortVectors.clear();
}

void CRig::Solve (AFITCoin x) {

    uint64_t fixCount;
    int varCount;
    vector<uint32_t> shortVectors;
    vector<int> indexes;

    //prepare vectors
    findVarBits(x, fixCount, varCount, indexes);
    removeFixBits(x, fixCount, indexes, shortVectors);

    if (x->m_Vectors.size() == 1) {
        x->m_Count = combinationSum(x->m_DistMax, fixCount); //for a single number given, we just compute the result
        shortVectors.clear();
    } else {
        compareWithNumbers(x, varCount,fixCount, shortVectors);
    }

    FITtestCounter++;
  //  printf("FITtest %d: result: %zu\n", FITtestCounter, x->m_Count);

}

//--------------------------------------- CVUTCoin methods -----------------------------------------------

int CRig::min(int x, int y, int z) {
    return min(min(x, y), z);
}

//via geeksforgeeks
int CRig::editDistance(vector<bool>& str1, vector<bool>& str2) {

    int m,n;
    m = str1.size();
    n = str2.size();

    // Create a table to store results of subproblems
    vector< vector<int> > dp(m + 1, vector<int>(n + 1)); //a vector version of int dp[m+1][n+1];

    // Fill d[][] in bottom up manner
    for (int i=0; i<=m; i++)
    {
        for (int j=0; j<=n; j++)
        {
            // If first string is empty, only option is to insert all characters of second string
            if (i==0)
                dp[i][j] = j;  // Min. operations = j

            // If second string is empty, only option is to remove all characters of second string
            else if (j==0)
                dp[i][j] = i; // Min. operations = i

            // If last characters are same, ignore last char and recur for remaining string
            else if (str1[i-1] == str2[j-1])
                dp[i][j] = dp[i-1][j-1];

            // If last character are different, consider all possibilities and find minimum
            else
                dp[i][j] = 1 + min(dp[i][j-1],  // Insert
                                   dp[i-1][j],  // Remove
                                   dp[i-1][j-1]); // Replace
        }
    }

    return dp[m][n];
}

void CRig::prepareData(ACVUTCoin &x, vector<bool> & boolVector) { // move all the bits into a single array of bools
    bool newBit;
    uint8_t currentNum;

     for (uint64_t i = 0; i != x->m_Data.size(); i++) { //for each byte in input data
        currentNum = x->m_Data[i];
        for (int j = 0; j < 8; j++) { //for each bit in byte given, put it in the vector
            newBit = currentNum & 1;
            boolVector.push_back(newBit);
            currentNum >>= 1;
        }
     }
}

void CRig::findPrefixSuffix(ACVUTCoin &x, vector<bool> & boolVector) {
    int dist;

    for (uint64_t i = 1; i <= boolVector.size(); i++) { //find all prefixes, i must not be 0
        vector<bool> prefix(boolVector.begin(), boolVector.begin() + i);
        for (uint64_t j = 0; j < boolVector.size(); j++) { //find all suffixes, j must not be max
            vector<bool> suffix(boolVector.begin()+ j, boolVector.end());

            dist = editDistance(prefix, suffix);
            if(dist >= x->m_DistMin && dist <= x->m_DistMax) {
                x->m_Count++;
            }

        }

    }


}

void CRig::Solve (ACVUTCoin x) {
    vector<bool> boolVector;

    prepareData(x, boolVector);
    findPrefixSuffix(x, boolVector);

    CVUTtestCounter++;
   // printf("CVUTtest %d: result: %zu\n", CVUTtestCounter, x->m_Count);
}

CRig::CRig (void) {
    FITtestCounter = 0;
    CVUTtestCounter = 0;
    bitsLen = 32;
    customerIndex = 0;
    endFlag = false;
    customerCounter = 0;
    workCount = 0;
    sem_init ( &semCoinEmpty, 0, 100 );
    sem_init ( &semCoinFull, 0, 0 );
}

//--------------------------------------- Parallel solution metods ----------------------------------------


void CRig::AddFitCoins(ACustomer c, int custIdx) {
    printf("----Inside addfitCoins method %d\n", custIdx);

    int i = 1;
    for ( AFITCoin x = c -> FITCoinGen (); x ; x = c -> FITCoinGen () ) {
        CCoin newFitCoin(true, x, nullptr, custIdx, custIdx * 100 + i);
     //   printf("[%d]----Adding a fitCoin %d\n",custIdx, custIdx * 100 + i);

        sem_wait(&semCoinEmpty);

        mtxCoinBuffer.lock();
            coinBuffer.push_back(newFitCoin);
            workCount++;
        mtxCoinBuffer.unlock();
        i++;

        sem_post(&semCoinFull); //zvyseni semaforu, signal konzumentovi SolveCoin
    }

    printf("----added all fit coins from customer %d\n", custIdx);
}

void CRig::AddCvutCoins(ACustomer c, int custIdx) {
    printf("-----Inside addcvutCoins method %d\n", custIdx);

    int i = 1;
    for ( ACVUTCoin x = c -> CVUTCoinGen (); x ; x = c -> CVUTCoinGen () ) {
        CCoin newCvutCoin(false, nullptr, x, custIdx, custIdx * 1000 + i);
      //  printf("[%d]----Adding a cvutCoin %d\n",custIdx, custIdx * 1000 + i);

        sem_wait(&semCoinEmpty);

        mtxCoinBuffer.lock();
            coinBuffer.push_back(newCvutCoin);
            workCount++;
        mtxCoinBuffer.unlock();
        i++;

        sem_post(&semCoinFull); //zvyseni semaforu, signal konzumentovi SolveCoin
    }


    printf("-----added all cvut coins from customer %d\n", custIdx);
}

void CRig::AcceptCoin(ACustomer c, int idx) {
    printf("inside AcceptCoin %d\n", idx);

    while (1) {
        if (customers[idx].solvedCoins. size () > 0) { //there are some solved problems in the buffer to be accepted
            //vybrat
            sem_wait(& (customers[idx].semCusFull));

            mtxSolved.lock();
                CCoin coin = customers[idx].solvedCoins.front();
                customers[idx].solvedCoins.pop_front();
            mtxSolved.unlock();

            //odevzdat
            if(coin.isFit) {
                printf("[%d] accepting fitCoin %d (result %zu)\n", idx, coin.ID, coin.fitCoin->m_Count);
              //  coin.printCoin();
                c -> FITCoinAccept (coin.fitCoin);
            }
            else {
                printf("[%d] accepting cvutCoin %d (result %zu) \n", idx, coin.ID, coin.cvutCoin->m_Count);
              //  coin.printCoin();
                c -> CVUTCoinAccept (coin.cvutCoin);
            }

            mtxSolved.lock();
                workCount--;
            mtxSolved.unlock();

            sem_post(& (customers[idx].semCusEmpty));

        }
        //konec
        if(endFlag && (workCount == 0) ) { //buffer je prazdny a zaroven fce Stop() nastavila endflag
            break;
        }
    }

}

void CRig::AddCustomer (ACustomer c) {

   // printf("---AddCustomer\n");
    CustomerWrapper customer(c);


    customerThreads . push_back ( thread (&CRig::AddFitCoins, this, c, customerIndex) ); //thread to add fitcoins
    customerThreads . push_back (  thread (&CRig::AddCvutCoins, this, c, customerIndex) ); //thread to add cvutcoins
    customerThreads . push_back (  thread (&CRig::AcceptCoin, this,  c, customerIndex) ); //thread to accept coins of both types


    customers.push_back(customer);

    customerIndex++; //index to tell which customer the coin belongs to
   // printf("customer index = %d\n", customerIndex);

}

void CRig::SolveCoin() {
    printf("-----SolveCoin\n");
    while(1) {

        //semafor
        sem_wait(&semCoinFull); //global semaphore of coinBuffer

            //vybrat
            mtxCoinBuffer.lock();
                CCoin coin = coinBuffer.front();
                coinBuffer.pop_front();
            mtxCoinBuffer.unlock();

            //solve
            if(coin.isFit)
                Solve(coin.fitCoin);
            else
                Solve(coin.cvutCoin);

           sem_wait(& (customers[coin.customerIdx].semCusEmpty)); //personal semaphore of customer

           //ulozit
            mtxSolved.lock();
                customers[coin.customerIdx].solvedCoins.push_back(coin);
            mtxSolved.unlock();

        sem_post(&semCoinEmpty); //global semaphore of coinBuffer
        sem_post(& (customers[coin.customerIdx].semCusFull)); //personal semaphore of customer

        //konec
        if((coinBuffer . size () == 0) && endFlag) { //buffer je prazdny a zaroven fce Stop() nastavila endflag
            printf("End of SolveCoin\n");
            break;
        }
    }


}

void CRig::Start (int thrCnt) {
    //create threads
    for ( int i = 0; i < thrCnt; i ++ )
      workThreads . push_back ( thread ( &CRig::SolveCoin, this ) );
}

void CRig::Stop (void) {
    printf("------STOP---Calling stop function\n");
    endFlag = true; // indikator konce pro funkci SolveCoin


            //wait for customer threads
            for ( auto & th : customerThreads ) {
              th . join ();
            }


            //wait for working threads
            for ( auto & t : workThreads ) {
                t . join ();
            }




    if (coinBuffer . size () > 0)
        printCustomers();
}



#ifndef __PROGTEST__
#include "test.cpp"
#endif /* __PROGTEST__ */
