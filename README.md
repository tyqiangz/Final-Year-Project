# Final-Year-Project
# Quasi-Cyclic Moderate Density Parity-Check (QC-MDPC) codes

### Functions built:

```genFirstRow(r, wi)```
Input: Integers r, wi  
Output: Numpy array of length r and Hamming weight wi

```genCirculant(firstRow)```  
Input: First row of a circulant matrix  
Output: Circulant matrix based on the inputed first row

```genTransposePoly(firstRow)```  
Input: First row of a circulant matrix, H  
Output: The first row of the transpose of H

```genSumPoly(firstRowA, firstRowB)```  
Input: First row of circulant matrices A and B
Output: First row of circulant matrix A+B

