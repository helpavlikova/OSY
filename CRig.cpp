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
    static void bin(unsigned n, int bitsLen);
    static void printBinary(int i, int j, int diff, int mask);
    static unsigned int countSetBits(int n);
    static uint64_t binomialCoeff(uint64_t n, uint64_t k);
    static void prepareVectors(AFITCoin &x);
    static void prepareData(ACVUTCoin &x);
    static void findPrefixSuffix(ACVUTCoin &x);
    static void printVectors(vector<int>& vectors);
    static void printVectors(vector<bool>& vectors);
    static void printVectors(vector<uint32_t>& vectors);
    static void printVectors(vector<uint8_t>& vectors);
    static int editDistance(vector<bool>&  word1, vector<bool>& word2);
    static int swapBits(unsigned int n, unsigned int p1, unsigned int p2);
    static uint64_t combinationSum(uint64_t k);
    static void findVarBits(AFITCoin &x);
    static void removeFixBits(AFITCoin &x);
    static void compareWithNumbers(AFITCoin &x);
    static int varCount;
    static int testCounter;
    static uint64_t fixCount;
    static int bitsLen;
    static vector<uint32_t> shortVectors;
    static vector<int> indexes;
    static vector<bool> boolVector;
};

//declaration of static variables
int CRig::varCount;
uint64_t CRig::fixCount;
int CRig::testCounter = 0;
int CRig::bitsLen = 32;
vector<uint32_t> CRig::shortVectors;
vector<int> CRig::indexes;
vector<bool> CRig::boolVector;


//--------------------------------------- printing methods -----------------------------------------------

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
//--------------------------------------- FITCoin methods -----------------------------------------------

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

uint64_t CRig::combinationSum(uint64_t k) { //kombinatoricka magie z prednasky
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

void CRig::findVarBits(AFITCoin &x) {
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

void CRig::removeFixBits(AFITCoin &x) {
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

    indexes.clear();
}

void CRig::prepareVectors(AFITCoin &x) {
    findVarBits(x);
    removeFixBits(x);
}

void CRig::compareWithNumbers(AFITCoin &x) {
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
            x->m_Count += combinationSum(x->m_DistMax - maxDiff);
        }
    }

    shortVectors.clear();
}

void CRig::Solve (AFITCoin x) {
    testCounter++;
    printf("test %d: \n", testCounter);

    prepareVectors(x);

    if (x->m_Vectors.size() == 1) {
        x->m_Count = combinationSum(x->m_DistMax); //for a single number given, we just compute the result
        shortVectors.clear();
    } else {
        compareWithNumbers(x);
    }

}

//--------------------------------------- CVUTCoin methods -----------------------------------------------
int min(int x, int y, int z)
{
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

void CRig::prepareData(ACVUTCoin &x) { // move all the bits into a single array of bools
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

     //printVectors(boolVector);
}

void CRig::findPrefixSuffix(ACVUTCoin &x) {
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

    testCounter++;
    printf("test %d: \n", testCounter);

    prepareData(x);
    findPrefixSuffix(x);
    boolVector.clear();

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
