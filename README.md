# Final-Year-Project
# Quasi-Cyclic Moderate Density Parity-Check (QC-MDPC) codes

A research on optimising the QC-MDPC code for use in the McEliece cryptosystem for my Bachelor's Thesis in the NUS Math department

- [x] Functions for matrix operations on QC-MDPC matrices
- [ ] Bit-Flipping algorithm
- [ ] Sum-Product algorithm

<br>

**Functions built:**  

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

```convertNumpyToSympy(f)```  
**Input**: Numpy array containing coefficients of desired polynomial f(x) in ascending power of x  
**Output**: Sympy polynomial f

```convertSympyToNumpy(f)```  
**Input**: Sympy polynomial f  
**Output**: Numpy array containing coefficients of f(x) in ascending power of x  

```genParityCheck(n, r, w)```  
**Input**: Integers n (length of QC-MDPC code), r (length of each ciruclant block), w (sum of weight of all circulant blocks)  
**Output**: (n, r, w)-QC-MDPC matrix

```count4Cycles(H, n, r, w)```  
**Input**: (n, r, w)-QC-MDPC matrix  
**Output**: Number of 4 cycles in the Tanner graph of the input (n, r, w)-QC-MDPC matrix

```drawTanner(H)```  
**Input**: Parity-check matrix H, not necessarily QC-MDPC  
**Output**: Tanner graph of H



