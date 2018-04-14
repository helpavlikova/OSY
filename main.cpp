#ifndef __PROGTEST__

#include <cstdio>
#include <cstdlib>
#include <bitset>
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
class CFITCoin{

public:
    CFITCoin ( const vector<uint32_t> & vectors, int distMax ):m_Vectors ( vectors ), m_DistMax ( distMax ), m_Count ( 0 ) {
    }
    virtual                  ~CFITCoin                     ( void ) = default;
    vector<uint32_t>         m_Vectors;
    int                      m_DistMax;
    uint64_t                 m_Count;
};
//=================================================================================================
class CCVUTCoin{

public:
    CCVUTCoin ( const vector<uint8_t> & data, int distMin, int distMax )
            : m_Data ( data ),
              m_DistMin ( distMin ),
              m_DistMax ( distMax ),
              m_Count ( 0 ) {
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

#define  SEM_INIT_VARIABLE 80



class Problem;

class MyCustomer{
public:
    MyCustomer();
    ~MyCustomer();
    deque <Problem> solvedQueue;
    int workingProblems = 0;
    int id  = 0;
    sem_t sem_free_end;
    sem_t sem_full_end;


};

MyCustomer::MyCustomer() {
    sem_init(&sem_free_end, 0, SEM_INIT_VARIABLE);
    sem_init(&sem_full_end, 0, 0);

}

MyCustomer::~MyCustomer() {
    sem_destroy(&sem_free_end);
    sem_destroy(&sem_full_end);

}


class Problem{
public:
    AFITCoin afitCoin= NULL;
    ACVUTCoin acvutCoin = NULL;
    bool lastProblem = false;
    int CustomerId = 0;


};




class CRig{
public:
    static void Solve ( ACVUTCoin x );
    static void Solve ( AFITCoin x );

    vector <bool> Differences ( AFITCoin x);
    int hammingWeight (vector<bool> x);
    vector <uint32_t> prepareVector( vector<uint32_t > x, vector <bool> y);
    uint64_t Solve ( vector<uint32_t> preparedVectors, int hammingweight, int distMax);
    uint64_t Pascal ( uint32_t m, uint32_t n);
    void FillArray();
    int levenshteinDistance( const string &s1, const string &s2);
    string convertToBinaryString(vector<uint8_t>);
    void SolveThread();
    void SolveFit(uint32_t , ACustomer);
    void SolveCvut(uint32_t , ACustomer);
    void PassResult(uint32_t  ,ACustomer);


    CRig                          ( void );
    ~CRig                         ( void );
    void Start( int thrCnt );
    void Stop ( void );
    void AddCustomer ( ACustomer c );

private:

    uint64_t array[33][33];
    vector <thread> workers;
    vector <MyCustomer> customers;
    deque <Problem> problemQueue;
    int numOfWorkers = 0;
    uint32_t numOfCustomers = 0;
    uint32_t numOfCustomersId = 0;
    uint32_t numOfCustomersProblems = 0;
    vector <thread> customerThreads;



    sem_t sem_full, sem_free;
    sem_t sem_free_end, sem_full_end;
    mutex mtxCnt, mtxCustThr, mtxCust,  mtxQueue, mtxEndQ, mtxObj;

};








void CRig::Solve(AFITCoin x) {
//    cout << " START " << endl;
    CRig fitCoin;

    vector <bool> m_isDifferent;
    m_isDifferent = fitCoin.Differences(x);
    int hammingweight = fitCoin.hammingWeight(m_isDifferent);
/*
    cout << " M IS DIFFERENT" << endl;
    for (uint i = 0; i < m_isDifferent.size();++i){
        cout << m_isDifferent[i];
    }
    cout << endl;
*/
    vector <uint32_t > preparedVectors = fitCoin.prepareVector(x->m_Vectors, m_isDifferent);
//    cout << "prepared vectors " << endl;
//    for (uint i = 0; i < preparedVectors.size(); ++i) {
//        cout << bitset<32>(preparedVectors[i]) << endl;
//    }
//    cout << endl;

//    cout << endl << "hamming weight  = " << hammingweight << endl;

    uint64_t result = fitCoin.Solve(preparedVectors,hammingweight, x->m_DistMax);
    x->m_Count = result;



//    cout << " END "<< endl;

}

void CRig::Solve(ACVUTCoin x) {
    CRig cvutCoin;
    uint64_t cnt=0;


    string s = cvutCoin.convertToBinaryString(x->m_Data);

//    cout << s << endl;

    for(uint prefix = 0; prefix < s.length();prefix++){
        for(uint suffix = 0; suffix < s.length(); suffix++){
            string s1, s2;
            s1 = s.substr(0, prefix+1);
            s2 = s.substr(s.length()-1-suffix,suffix+1);
            int result = cvutCoin.levenshteinDistance(s1, s2);
            if(result<=x->m_DistMax && result >= x->m_DistMin)
                cnt++;
//            cout << "s1 " << s1 << " s2 " << s2 << " result " << result << endl;

        }
    }
//    1101110001001000100110000010001101000101
    //cout << "cnt " << cnt << endl;

    x->m_Count=cnt;

}

CRig::CRig(void) {

    FillArray();
    sem_init(&sem_free, 0, SEM_INIT_VARIABLE);
    sem_init(&sem_full, 0, 0);


}

CRig::~CRig(void) {
    sem_destroy(&sem_full);
    sem_destroy(&sem_free);

}



vector<bool> CRig::Differences(AFITCoin x) {
    vector <bool> differences;

//    cout << "DistMax = " << x->m_DistMax<< endl;

//    for (uint i = 0; i < x->m_Vectors.size(); ++i) {
//        cout <<bitset<32> (x->m_Vectors[i]) <<endl;
//    }
//    cout << endl;


//    cout << endl << "bit shift " << endl;

    vector <uint32_t> shifted;
    for(int i =0; i <32 ;i++){
        for( uint n = 0; n < x->m_Vectors.size();n++){
            uint32_t num = x->m_Vectors[n];
            uint32_t shiftedNum = num << (31-i);
            uint32_t shiftedNumBack = shiftedNum >> 31;
            shifted.push_back(shiftedNumBack );

        }
        uint8_t bit = shifted[0];
        for(uint j = 0; j< shifted.size();j++){
            if(bit != shifted[j]) {
                differences.push_back(true);
                break;
            }
            if(j == shifted.size()-1){
                differences.push_back(false);
            }
        }
        shifted.clear();
    }
//    cout << endl;

    for(uint i = 0; i<differences.size() ;i++){
//        if(differences[i] == true)
//            cout << 1;
//        else
//            cout << 0;
    }
//    cout << endl;
/*
    cout << "shifted array" << endl;
    for(int i = 0 ; i < shifted.size(); ++i){
        cout <<(shifted[i]) << endl;
    }
    */

    return differences;
}

int CRig::hammingWeight(vector <bool> v) {
    int x=0;
    for(uint i = 0; i < v.size();++i){
        if (v[i] == true){
            x+=1;
            if(i == v.size()-1) break;
            x<<=1;

        } else {
            if(i == v.size()-1) break;
            x<<=1;
        }
    }

    // cout << "X in binary "<< bitset<32>(x) << endl;

    int k = x ? 1 : 0;
    while (x &= (x-1))++k;
    return k;
}

vector<uint32_t> CRig::prepareVector(vector<uint32_t> x, vector <bool> y) {
    uint32_t number = 0;
    //int cnt=0;
    vector<uint32_t> prepared;
    for( uint vector = 0; vector < x.size();vector++){
        for (uint difference = 0; difference < y.size(); difference++){
            if(y[difference] == true){
                //cout << " x[vector]" << bitset<32> ( x[vector]  )<< endl;
                //cout << "difference" << difference << endl;
                uint32_t shift1 = x[vector] << (31 - difference);
                //cout << "shift1" << bitset<32>(shift1) << endl;
                uint32_t shift2 = shift1 >> 31;
                //cout << "shift2" << bitset<32>(shift2) << endl;
                number |= shift2;

                //cout << "cnt "<<endl <<++cnt << endl;
                //cout <<"number" << bitset<32>(number) << endl;
                number <<= 1;


            }
        }
        number >>=1;
        prepared.push_back(number);
        number = 0;
    }
    return prepared;
}

uint64_t CRig::Solve(vector<uint32_t> preparedVectors, int hammingweight, int distMax) {
//    cout << "solve start" << endl;
    vector <int> maxValues;
    int distance = 0;
    int max = 0;
    uint32_t maxaval = pow (2,hammingweight);
    for (uint32_t i = 0; i < maxaval;i++){
        for(uint32_t j = 0; j < preparedVectors.size(); j++){
            uint32_t val = i ^ preparedVectors[j];
            while (val!=0){
                val = val & (val -1);
                distance++;
            }
            if(distance>max)
                max = distance;

            distance = 0;
        }
        //cout << "max is " << max << endl;
        maxValues.push_back(max);
        max = 0;
    }
//    cout << maxValues.size() << endl;
//    cout << "maxvalues" << endl;
//    for (uint  i = 0; i < maxValues.size(); i++) {
//        cout << maxValues[i]<< " ";
//    }

//    cout << "solve mid" <<endl;


    uint64_t sum = 0;


//    cout << "data for sum " << endl;


    for(uint32_t i = 0; i < maxValues.size();i++){
        for(int  j = 0; j <= distMax-maxValues[i];j++){
//            cout << "sum from 0 to " << distMax << " - " << maxValues[i] << " of m = 32 - " << hammingweight << " over n = " << j << endl;
//            uint32_t  result = Pascal(32 - hammingweight, j);
            if( (32-hammingweight) < j) continue;
//            uint32_t result = Pascal(32-hammingweight,j);
            uint32_t  result = array[32-hammingweight][j];
            sum+=result;
//            cout << "result " << result << endl;
        }
    }

//    cout << " SUM IS " << sum << endl;




    return sum;
}
// CODE FOR THIS FUCNTION IS FROM https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/ 6.4.2018
uint64_t CRig::Pascal(uint32_t n, uint32_t k) {


    uint64_t res = 1;

    if ( k > n - k )
        k = n - k;

    for (uint i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;


}

void CRig::FillArray() {
    uint64_t i,j,x;
    for(i=0;i<33;i++)
    {
        x=1;
        for(j=0;j<=i;j++)
        {
            array[i][j] = x;
            x = x * (i - j) / (j + 1);
//            cout << Pascal(i,j) << " ";
//            cout <<array[i][j] << " ";
//            cout << "[" <<i <<"][" << j<< "]"<<array[i][j] << " ";
        }
//        cout << endl;
    }


}
//-------------------THIS ALGORITHM AND CODE IS FROM https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++ 7.4.2018
int CRig::levenshteinDistance(const string &s1, const string &s2) {
    int s1len = s1.size();
    int s2len = s2.size();

    auto column_start = (decltype(s1len))1;

    auto column = new decltype(s1len)[s1len + 1];
    iota(column + column_start - 1, column + s1len + 1, column_start - 1);

    for (auto x = column_start; x <= s2len; x++) {
        column[0] = x;
        auto last_diagonal = x - column_start;
        for (auto y = column_start; y <= s1len; y++) {
            auto old_diagonal = column[y];
            auto possibilities = {
                    column[y] + 1,
                    column[y - 1] + 1,
                    last_diagonal + (s1[y - 1] == s2[x - 1]? 0 : 1)
            };
            column[y] = min(possibilities);
            last_diagonal = old_diagonal;
        }
    }
    auto result = column[s1len];
    delete[] column;
    return result;


}


string CRig::convertToBinaryString(vector<uint8_t > x) {
    string result = "";


    for(uint i = x.size();i-- > 0;){
        for(uint j = 0; j < 8 ; j++){
            uint8_t shift1 = x[i] >> (7-j);
            uint8_t shift2 = shift1 << 7;
            uint8_t shift3 = shift2 >> 7;
            result += to_string(shift3);
        }
    }


//    cout << endl;
//    for(uint i = 0; i < 8 ; i++){
//        cout << " number" <<bitset<8>(number) << endl;
//        uint8_t shift1 = number >> (7-i);
//        cout << " shift 1" <<bitset<8>(shift1) << endl;
//        uint8_t shift2 = shift1 << 7;
//        uint8_t shift3 = shift2 >> 7;
//        cout << " shift 2" <<bitset<8>(shift2) << endl;
//        result += to_string(shift3);
//        cout << " result " << result;
//    }
//    cout << endl;
    return result;
}



void CRig::Start(int thrCnt) {
    unique_lock<mutex> lock(mtxCnt);
    numOfWorkers = thrCnt;

    for(int i = 0 ; i < thrCnt; i++){
        workers.push_back(thread(&CRig::SolveThread, this));
    }
    lock.unlock();


}

void CRig::Stop(void) {


    for(uint i = 0; i < workers.size(); i++) {
        workers[i]. join();
    }

    for(uint i = 0; i < customerThreads.size();i++){
        customerThreads[i].join();
    }


}


void CRig::AddCustomer(ACustomer c) {




    unique_lock<mutex> lock3 (mtxCust);
    numOfCustomers++;
    lock3.unlock();

    MyCustomer customer;


    customers.push_back(customer);

    unique_lock<mutex> lock (mtxCnt);
    numOfCustomersProblems+=2;
    numOfCustomersId++;
    unique_lock<mutex> lock2(mtxCustThr);
    customerThreads.push_back(thread (&CRig::SolveFit, this, numOfCustomers-1, c));
    customerThreads.push_back(thread (&CRig::PassResult,this, numOfCustomers-1, c));
    customerThreads.push_back(thread (&CRig::SolveCvut,this, numOfCustomers-1, c));
    lock2.unlock();
    lock.unlock();



}



void CRig::PassResult(uint32_t customerVectorIndex, ACustomer c) {
    while(1){

        sem_wait(&customers[customerVectorIndex].sem_full_end);
        unique_lock<mutex> lock ( mtxEndQ);





//
//        if(customers[customerVectorIndex].solvedQueue.front().acvutCoin == NULL && customers[customerVectorIndex].solvedQueue.front().afitCoin == NULL){
//            break;
//        }


        if(customers[customerVectorIndex].solvedQueue.front().acvutCoin){
            cout << " successfully passed CVUT COIN " << this_thread::get_id() << endl;
            unique_lock<mutex> read (mtxObj);
            c->CVUTCoinAccept(customers[customerVectorIndex].solvedQueue.front().acvutCoin);
            read.unlock();
            if(customers[customerVectorIndex].solvedQueue.front().lastProblem){
                customers[customerVectorIndex].solvedQueue.pop_front();
                break;
            }
            customers[customerVectorIndex].solvedQueue.pop_front();

        }
        if(customers[customerVectorIndex].solvedQueue.front().afitCoin){
            unique_lock <mutex> read(mtxObj);
            cout << " successfully passed FIT COIN " << this_thread::get_id() << endl;
            c->FITCoinAccept(customers[customerVectorIndex].solvedQueue.front().afitCoin);
            read.unlock();
            if(customers[customerVectorIndex].solvedQueue.front().lastProblem) {
                customers[customerVectorIndex].solvedQueue.pop_front();
                break;
            }
            customers[customerVectorIndex].solvedQueue.pop_front();

        }
        lock.unlock();
        sem_post(&customers[customerVectorIndex].sem_free_end);




    }





}

void CRig::SolveThread() {
    while(1){

        Problem problem;

        sem_wait(&sem_full);

        unique_lock<mutex> lock2 (mtxQueue);
        problem = problemQueue.front();
        problem.CustomerId = problemQueue.front().CustomerId;

        problemQueue.pop_front();
        lock2.unlock();

        sem_post(&sem_free);




//        cout << "suze" << problem.myCustomer->solvedQueue.size() << endl;

        if(problem.acvutCoin == nullptr && problem.afitCoin == nullptr){
//            cout << " THREAD ENDED " << endl;
            break;
        }

        if(problem.acvutCoin){
            Solve(problem.acvutCoin);
        } else {
            Solve(problem.afitCoin);
        }
        unique_lock <mutex> lock (mtxCnt);
        customers[problem.CustomerId].workingProblems--;


        if(customers[problem.CustomerId].workingProblems == 0){
            problem.lastProblem = true;
        }
        lock.unlock();


        sem_wait(&customers[problem.CustomerId].sem_free_end);
        unique_lock<mutex> lock4(mtxEndQ);
        customers[problem.CustomerId].solvedQueue.push_back(problem);
        lock4.unlock();

        sem_post(&customers[problem.CustomerId].sem_full_end);

    }

}

void CRig::SolveFit(uint32_t customerVectorIndex, ACustomer c) {
    Problem problem;
    problem.CustomerId = customerVectorIndex;
    while(1){

        unique_lock<mutex> read (mtxObj);
        AFITCoin afitCoin;
        afitCoin = c->FITCoinGen();
        read.unlock();
        problem.afitCoin = afitCoin;

        if(afitCoin == NULL){
            {
                unique_lock <mutex> lock (mtxCnt);
                numOfCustomersProblems--;

                if(numOfCustomersProblems ==0 ){
                    for(int i = 0 ; i < numOfWorkers;i++){
                        sem_wait(&sem_free);
                        unique_lock<mutex> lock2(mtxQueue);
                        problemQueue.push_back(problem);
                        lock2.unlock();
                        sem_post(&sem_full);
                    }
                }
                lock.unlock();
            }
            break;
        }

        unique_lock <mutex> lock (mtxCnt);
        customers[customerVectorIndex].workingProblems++;
        lock.unlock();

        sem_wait(&sem_free);
        unique_lock<mutex> lock2(mtxQueue);
        problemQueue.push_back(problem);
        lock2.unlock();
        sem_post(&sem_full);

    }



}

void CRig::SolveCvut(uint32_t customerVectorIndex, ACustomer c) {
    Problem problem;
    problem.CustomerId = customerVectorIndex;
    while(1){
        unique_lock<mutex> read (mtxObj);
        ACVUTCoin acvutCoin;
        acvutCoin = c->CVUTCoinGen();
        read.unlock();



        problem.acvutCoin = acvutCoin;


        if(acvutCoin == NULL){
            unique_lock <mutex> lock (mtxCnt);
            numOfCustomersProblems--;

            if(numOfCustomersProblems ==0 ){
                for(int i = 0 ; i < numOfWorkers;i++){
                    sem_wait(&sem_free);
                    unique_lock<mutex> lock2(mtxQueue);
                    problemQueue.push_back(problem);
                    lock2.unlock();
                    sem_post(&sem_full);
                }
            }
            lock.unlock();
            break;
        }

        unique_lock <mutex> lock (mtxCnt);
        customers[customerVectorIndex].workingProblems++;
        lock.unlock();

        sem_wait(&sem_free);
        unique_lock<mutex> lock2(mtxQueue);
        problemQueue.push_back(problem);
        lock2.unlock();
        sem_post(&sem_full);

    }


}

#ifndef __PROGTEST__
#include "test.inc"
#endif /* __PROGTEST__ */
