#include <Rcpp.h>
using namespace Rcpp;

#include <string.h>
#include <vector>
using namespace std;
// #include <stdio.h>

// This code is loosely based on functions from the pcaPP package.
// See https://github.com/cran/pcaPP/blob/master/src/cov.kendall.cpp

void insertionSort(double* y, double* data, int* ids, int* swapCounts, size_t len) {
  size_t maxJ, i, tieCount;
  
  maxJ = len - 1;
  tieCount = 0; 
  
  if(len < 2) {
    return;
  }
  
  for(i = len - 2; i < len; --i) {
    size_t j = i;
    
    double tempVal = y[i];
    int tempId = ids[i];
    for(; j < maxJ && y[j + 1] < tempVal; ++j) {
      y[j] = y[j + 1];
      ids[j] = ids[j + 1];
      if(data[ids[j]] != data[tempId]){
        swapCounts[ids[j]] += 1;
      }
    }
    
    y[j] = tempVal;
    ids[j] = tempId;
    for(size_t k=i; k< j; ++k){
      if(data[ids[k]] == data[tempId]){
        tieCount++;
      }
    }
    swapCounts[tempId] += (j - i);
    swapCounts[tempId] -= tieCount; 
    tieCount = 0; 
  }
}


static void merge(double* y, double* x, double* data, 
                  int* ids, int* swapCounts,
                  size_t middle, size_t len) {
  
  vector<int> copyIds;
  double* left;
  double* right;
  size_t bufIndex, leftLen, rightLen, tieCount;
  
  
  for(size_t k=0; k<len; ++k) copyIds.push_back(ids[k]);
  left = y;
  right = y + middle;
  leftLen = middle;
  rightLen = len - middle;
  bufIndex = 0;
  tieCount = 0; 
  
  while(leftLen && rightLen) {
    
    if(right[0] < left[0]) {
      tieCount = 0; 
      for(size_t i = middle - leftLen; i < middle; i++) {
        if(data[copyIds[i]] == data[copyIds[len-rightLen]]) {
          tieCount++;
        }
      }
      x[bufIndex] = right[0];
      ids[bufIndex] = copyIds[len - rightLen];
      swapCounts[copyIds[len - rightLen]] += leftLen;
      swapCounts[copyIds[len - rightLen]] -= tieCount;
      rightLen--;
      right++;
    } else {
      tieCount = 0;
      x[bufIndex] = left[0];
      ids[bufIndex] = copyIds[middle - leftLen];
      for(size_t i = middle; i < len- rightLen; i++) {
        if(data[copyIds[i]] == data[copyIds[middle-leftLen]]) {
          tieCount++;
        }
      }
      swapCounts[copyIds[middle - leftLen]] += len - middle - rightLen;
      swapCounts[copyIds[middle - leftLen]] -= tieCount;
      leftLen--;
      left++;
    }
    bufIndex++;
  }
  
  if(leftLen) {
    for(size_t k=0; k<leftLen; ++k){
      tieCount = 0;
      for(size_t i = middle; i < len; i++) {
        if(data[copyIds[i]] == data[copyIds[middle-leftLen + k]]) {
          tieCount++;
        }
      }
      ids[bufIndex + k] = copyIds[middle - leftLen + k];
      swapCounts[copyIds[middle - leftLen + k]] += len - middle;
      swapCounts[copyIds[middle - leftLen + k]] -= tieCount;
    }
    memcpy(x + bufIndex, left, leftLen * sizeof(double));
  } else if(rightLen) {
    memcpy(x + bufIndex, right, rightLen * sizeof(double));
  }
  
  return;
}

void mergeSort(double* y, double* x, double* data, 
               int* ids, int* swapCounts,
               size_t len) {
  
  if(len < 2) {
    return;
  }
  if(len < 10) {
    insertionSort(y, data, ids, swapCounts, len);
    return;
  }
  
  size_t half;
  half = len / 2;
  mergeSort(y, x, data, ids, swapCounts, half);
  mergeSort(y + half, x + half, data, ids + half, swapCounts, len - half);
  merge(y, x, data, ids, swapCounts, half, len);
  
  memcpy(y, x, len * sizeof(double));
  return;
}



// [[Rcpp::export]]
NumericVector countSwaps(NumericVector y, NumericVector data) {
  
  int n = y.size();
  
  // Allocate memory for arrays
  double* y_ptr = new double[n];
  double* x_ptr = new double[n];
  double* data_ptr = new double[n];
  int* ids_ptr = new int[n];
  int* swapCounts_ptr = new int[n];
  
  // Initialize arrays
  for (int i = 0; i < n; i++) {
    y_ptr[i] = y[i];
    data_ptr[i] = data[i];
    x_ptr[i] = 0;
    ids_ptr[i] = i;
    swapCounts_ptr[i] = 0;
  }
  
  // Call mergeSort
  mergeSort(y_ptr, x_ptr, data_ptr, ids_ptr, swapCounts_ptr, n);
  
  // Copy sorted array to output vector
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = swapCounts_ptr[i];
  }
  
  // Free memory
  delete[] y_ptr;
  delete[] x_ptr;
  delete[] ids_ptr;
  delete[] data_ptr;
  delete[] swapCounts_ptr;
  
  return out;
}
