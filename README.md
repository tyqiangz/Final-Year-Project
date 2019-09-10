# Final-Year-Project
# Quasi-Cyclic Moderate Density Parity-Check (QC-MDPC) codes

A research on optimising binary QC-MDPC codes for use in the McEliece cryptosystem for my Bachelor's Thesis at the National University of Singapore, Department of Mathematics.
<br>
- [x] Algorithms for matrix operations on circulant matrices
- [x] Algorithm for counting 4 cycles in Tanner graph of a QC-MDPC code
- [x] Optimising the algorithm for counting 4 cycles
- [x] Bit-Flipping algorithm
- [x] Sum-Product algorithm
- [ ] Algorithm to assess error-correcting capability of a QC-MDPC code 

For the mathematical background of QC-MDPC codes, refer to sections 2 & 3 of [2].  For the mathematical background of circulant matrices, refer to section 2 of the document [1]. For the pseudocode of the bit-flipping algorithm and Sum-Product algorithm, refer to the Appendix, Algorithm 1 and 2 of the document [1]. For the mathematical background of the bit-flipping algorithm and Sum-Product algorithm, refer to section 5 of the document [1].  

References (available in this repo):  
1. A Code-Based Key Agreement Scheme using QC-MDPC Codes  
2. MDPC-McEliece - New McEliece Variants from Moderate Density Parity-Check Codes

Note: I am only working on binary QC-MDPC matrices, hence by default all input vectors and matrices are to be in binary wherever neccesary.  

---  
**API Specification:**  

```genFirstRow(r, wi)```  
**Input**: Integers r, wi  
**Output**: Numpy array of length r and Hamming weight wi

```genCirculant(firstRow)```  
**Input**: First row of a circulant matrix  
**Output**: Circulant matrix based on the inputed first row

```genTransposePoly(firstRow)```  
**Input**: First row of a circulant matrix, H  
**Output**: The first row of the transpose of H

```genSumPoly(firstRowA, firstRowB)```  
**Input**: First row of circulant matrices A and B  
**Output**: First row of circulant matrix A+B

```genProdPoly(firstRowA, firstRowB)```  
**Input**: First row of circulant matrices A and B  
**Output**: First row of circulant matrix AB

```genInvPoly(firstRow)```  
**Input**: First row of circulant matrix H  
**Output**: First row of circulant matrix inverse of H

```convertBinary(v)```  
**Input**: integer array v  
**Output**: integer array of v modulo 2  

```convertNumpyToSympy(f)```  
**Input**: Numpy array containing coefficients of desired polynomial f(x) in ascending power of x  
**Output**: Sympy polynomial f

```convertSympyToNumpy(f)```  
**Input**: Sympy polynomial f  
**Output**: Numpy array containing coefficients of f(x) in ascending power of x  

```genQCMDPC(n, r, w)```  
**Input**: Integers n (length of QC-MDPC code), r (length of each ciruclant block), w (sum of weight of all circulant blocks)  
**Output**: (n, r, w)-QC-MDPC matrix

```genGenQCMDPC(H)```  
**Input**: Parity-check matrix H
**Output**: The generator matrix G of H

```count4Cycles(H, n, r, w)```  
**Input**: (n, r, w)-QC-MDPC matrix  
**Output**: Number of 4 cycles in the Tanner graph of the input (n, r, w)-QC-MDPC matrix

```drawTanner(H)```  
**Input**: Parity-check matrix H, not necessarily QC-MDPC  
**Output**: Tanner graph of H

```bitFlipping(H, c, N)```  
**Input**: H (Parity-check matrix), c (word to be decoded), N (cutoff for the number of sum-product iterations)  
**Output**: if decoding is successful, return the decoded word, else return 0

```sumProduct(H, y, N, p)```  
**Input**: H (Parity-check matrix), y (word to be decoded), N (cutoff for the number of sum-product iterations), p (probability of a bit being digit 0)  
**Output**: if decoding is successful, return the decoded word, else return 0

```genRandomVector(k, t)```  
**Input**: k (dimension of vector), t (Hamming weight of vector)  
**Output**: A binary array of size k with Hamming weight t.

```encryptMcEliece(G, m, t)```  
**Input**: G (Generator Matrix of a QC-MDPC matrix), m (Plaintext), t (Hamming weight of the error array to be added to m)  
**Note**: If G is an a by b  matrix then m should be of length (b - a).  
**Output**: Encrypted message

```decryptMcEliece(H, y, method, N, p)```  
**Input**: H (QC-MDPC matrix), y (ciphertext), method (either 'BF' or 'SP', representing Bit-Flipping and Sum-Product resp.), N (max no. of decoding iterations), p (probability of bit error, only for method = 'SP')  
**Output**: decryptedText (decrypted text of length n - r)  
**Note**: If H is a r by n matrix then y should be of length n. decryptedText is the integer 0 if the decoding algorithm specified by ```method``` took more than N iterations. decryptedText is the integer -1 if the Sum-Product algorithm incurred computational error.

```decryptSuccess(plaintext, decryptedText)```  
**Input**: plaintext, decryptedText  
**Output**: true if plaintext == decryptedText element-wise and have same length, false otherwise
