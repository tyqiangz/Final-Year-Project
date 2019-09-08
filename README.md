# Final-Year-Project
# Quasi-Cyclic Moderate Density Parity-Check (QC-MDPC) codes

### Functions built:

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
