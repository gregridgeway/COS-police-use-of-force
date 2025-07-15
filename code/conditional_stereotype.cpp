// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;


// numerically stable log(exp(a) + exp(b))
double logSumExp(double a, double b) 
{
  double maxVal = std::max(a, b);
  return maxVal + std::log(std::exp(a - maxVal) + std::exp(b - maxVal));
}


// B.R. Heap (1963). Permutations by interchanges. The Computer Journal 6(3):293-298.
//  modified to avoid repeats if multiple officers use the same kind of force
double denomHeaps(std::vector<int>& viYTemp,
                  const std::vector<double>& vdLambda,
                  const std::vector<double>& vdS,
                  int index,
                  int n)
{
  // Base case: if index reaches the end, return 1.0
  if (index >= n) {
    return 1.0;
  }
  
  double dSum = 0.0;
  
  for(int i=index; i < n; i++) 
  {
    bool fCheck = true;
    
    // Check for duplicates to avoid repeats
    for(int j=index; j < i; j++) 
    {
      if(viYTemp[j] == viYTemp[i]) 
      {
        fCheck = false;
        break;
      }
    }
    
    if(fCheck) 
    {
      std::swap(viYTemp[index], viYTemp[i]);
      
      // Compute the factor for the current index
      double factor = std::exp(vdLambda[index] * vdS[viYTemp[index]]);
      
      // Recursive call to process the next index
      double recursiveResult = denomHeaps(viYTemp, vdLambda, vdS, index + 1, n);
      
      dSum += factor * recursiveResult;
      
      // Swap back to restore ivYTemp
      std::swap(viYTemp[index], viYTemp[i]);
    }
  }
  
  return dSum;
}


// starting point for denomHeaps...
//   Efficient to allocate viYTemp once rather than O(m!) times in recursion
double denomHeapsWrapper(const std::vector<int>& viY,
                         const std::vector<double>& vdLambda,
                         const std::vector<double>& vdS)
{
  std::vector<int> viYTemp = viY;
  return denomHeaps(viYTemp, vdLambda, vdS, 0, viYTemp.size());
}

// starting point for denomHeaps...
//   Efficient to allocate viYTemp once rather than O(m!) times in recursion
double logDenomHeapsWrapper(const std::vector<int>& viY,
                            const std::vector<double>& vdLambda,
                            const std::vector<double>& vdS)
{
  std::vector<int> viYTemp = viY;
  return std::log(denomHeaps(viYTemp, vdLambda, vdS, 0, viYTemp.size()));
}


// R interface to denomHeaps()
//   useful for testing
// [[Rcpp::export]]
double denomHeaps_R(IntegerVector ivY,
                    NumericVector nvLambda,
                    NumericVector nvS)
{
  std::vector<int>    viY      = Rcpp::as<std::vector<int>>(ivY);
  std::vector<double> vdLambda = Rcpp::as<std::vector<double>>(nvLambda);
  std::vector<double> vdS      = Rcpp::as<std::vector<double>>(nvS);
  
  // return denomHeapsWrapper(viY, vdLambda, vdS);
  return logDenomHeapsWrapper(viY, vdLambda, vdS);
}


// Discrete Fourier transform to compute denominator
//   not used, but useful for verifying calculations in testing
// Lin, Z., Wang, Y. & Hong, Y. The computing of the Poisson multinomial distribution and applications in ecological inference and machine learning. Comput Stat 38, 1851–1877 (2023). https://doi.org/10.1007/s00180-022-01299-0
//   Equation 7, 8, 9
//   time O(n^(m-1)), memory O(n*m)
//   Faster than no-repeat Heaps when n>8, still works when n=60
//   This implementation evaluates at one value of y
//      R dpmd() computes entire pdf using FFT
//      This implementation uses DFT
//      172x faster than dpmd() for n=20, m=4
//      215x faster than dpmd() for n=10, m=4
//      146x faster than dpmd() for n=5, m=4
//      40x  faster than dpmd() for n=2, m=4
// [[Rcpp::export]]
double logDenomDFT(NumericVector nvLambda,
                   NumericVector nvS,
                   IntegerVector ivY)
{
  int m0 = nvS.length();      // # force types
  int n0 = nvLambda.length(); // # officers
  int i=0;
  int j=0;
  
  double dSum = 0.0;
  
  // convert Y to tabulated count
  std::vector<double> vdX(m0, 0.0);
  for(i=0; i<n0; i++)
  {
    vdX[ivY(i)] += 1.0;
  }
  
  // convert lambda and s to Poisson-Multinomial prob matrix
  NumericMatrix nmP(n0,m0);
  double dLogK = 0.0;
  double dMaxE = 0.0;
  double dSumP = 0.0;
  for(i=0; i<n0; i++)
  { // for numerical stability subtract max value before exp()
    dMaxE = -std::numeric_limits<double>::infinity();
    dSumP = 0.0;
    for(j=0; j<m0; j++)
    {
      nmP(i,j) = nvLambda(i) * nvS(j);
      if(nmP(i,j)>dMaxE) dMaxE = nmP(i,j);
    }
    for(j=0; j<m0; j++)
    {
      nmP(i,j) = exp(nmP(i,j) - dMaxE);
      dSumP += nmP(i,j);
    }
    for(j=0; j<m0; j++)
    {
      nmP(i,j) /= dSumP;
    }    
    // dLogK is for scaling final result
    dLogK += dMaxE + log(dSumP);
  }
  
  // for generating the l vectors
  bool fCarry;
  bool fComplete;
  std::vector<double> vdL(m0-1, 0.0);
  
  // for calculation of q(l) and p(x)
  std::complex<double> cIW(0.0, 2.0*3.14159265358979323846/(n0+1.0));
  std::vector<std::complex<double>> vcQLtemp(m0-1);
  std::complex<double> cTemp(0.0,0.0);
  std::complex<double> cQL(1.0,0.0);
  
  fComplete = false;
  while(!fComplete)
  {
    for(i=0; i<m0-1; i++) 
    {
      vcQLtemp[i] = exp(cIW * vdL[i]);
    }
    
    // iQL = prod(P[,m0] + P[,1:(m0-1)] * iQLtemp1);
    cQL.real(1.0); cQL.imag(0.0); // 1+0i
    for(i=0; i<n0; i++)
    {
      cTemp.real(nmP(i,m0-1)); 
      cTemp.imag(0.0);
      for(j=0; j<m0-1; j++)
      {
        cTemp += nmP(i,j) * vcQLtemp[j];
      }
      cQL *= cTemp;
    }
    
    // dSum += get_real(iQL * exp(-1i*w * (vL' * vX[1:(m0-1)]))); 
    cTemp.real(0.0); cTemp.imag(0.0); // 0+0i
    for(i=0; i<m0-1; i++)
    {
      cTemp += vdL[i] * vdX[i];
    }
    cTemp = cQL * exp(-cIW*cTemp);
    dSum += cTemp.real();
    
    // generate next l: (0,0,0), (0,0,1), ..., (n0,n0,n0)
    j = m0-2;
    fCarry = true;
    while(fCarry)
    {
      vdL[j] = vdL[j] + 1.0;
      if(vdL[j] > n0)
      { // reset digit j to 0, increment digit j-1
        vdL[j] = 0.0;
        j = j-1;
        if(j<0) // signals completion
        {
          fCarry = false;
          fComplete = true;
        }
      } 
      else
      {
        fCarry = false;
      }
    }
  }
  
  if (dSum <= 0.0) {
    Rcout << "⚠️ dSum = " << dSum << " — clipped to small epsilon to avoid log error." << std::endl;
    dSum = 1e-9;
  }
  
  return dLogK - (m0-1.0) * std::log(n0+1) + std::log(dSum);
} // end logDenomDFT()



// log conditional likelihood log L(lambda,s)
double logCL(const std::vector<double>& vdLambda,
             int n0, // number of officers
             const std::vector<double>& vdS,
             const std::vector<int>& viY,
             int iHeapsLimit)
{
  double dReturnVal = 0.0;
  
  if(n0==1)
  {
    dReturnVal = 0.0;  
  }
  else if((n0==2) && (viY[0]==viY[1]))
  {
    dReturnVal = 0.0;
  }
  // Hard code for n0=2, 3x faster than DFT
  else if((n0==2) && (viY[0]!=viY[1])) 
  {
    dReturnVal = 
      vdS[viY[0]]*vdLambda[0] + vdS[viY[1]]*vdLambda[1] -
      logSumExp(vdS[viY[0]]*vdLambda[0] + vdS[viY[1]]*vdLambda[1],
                vdS[viY[1]]*vdLambda[0] + vdS[viY[0]]*vdLambda[1]);
  }
  // recursive no-repeat Heaps for n0s where Heaps faster than DFT
  else if((n0>2) && (n0<=iHeapsLimit))
  {
    dReturnVal = 0.0;
    // compute numerator log sum
    for(int i=0; i<n0; i++)
    {
      dReturnVal += vdS[viY[i]]*vdLambda[i];
    }
    dReturnVal -= std::log(denomHeapsWrapper(viY, vdLambda, vdS));
  }
  // discrete Fourier transform for larger n
  else
  {
    int m0 = vdS.size();      // number of force types
    int i=0;
    int j=0;
    
    double dSum = 0.0;
    double dNumLogSum = 0.0;
    
    // convert Y to tabulated count
    std::vector<double> vdX(m0, 0.0);
    for(i=0; i<n0; i++)
    {
      vdX[viY[i]] += 1.0;
    }
    
    // convert lambda and s to Poisson-Multinomial prob matrix
    std::vector<double> vdP;
    vdP.resize(n0*m0);
    double dLogK = 0.0;
    double dMaxE = 0.0;
    double dSumP = 0.0;
    for(i=0; i<n0; i++)
    { // for numerical stability subtract max value before exp()
      dMaxE = -std::numeric_limits<double>::infinity();
      dSumP = 0.0;
      for(j=0; j<m0; j++)
      {
        vdP[i*m0+j] = vdLambda[i] * vdS[j];
        if(vdP[i*m0+j]>dMaxE) dMaxE = vdP[i*m0+j];
      }
      for(j=0; j<m0; j++)
      {
        vdP[i*m0+j] = exp(vdP[i*m0+j] - dMaxE);
        dSumP += vdP[i*m0+j];
      }
      for(j=0; j<m0; j++)
      {
        vdP[i*m0+j] /= dSumP;
      }    
      // dLogK is for scaling final result
      dLogK += dMaxE + log(dSumP);
    }
    
    // compute numerator log sum
    for(i=0; i<n0; i++)
    {
      dNumLogSum += vdS[viY[i]]*vdLambda[i];
    }
    
    // for generating the l vectors
    bool fCarry;
    bool fComplete;
    std::vector<double> vdL(m0-1, 0.0);
    
    // for calculation of q(l) and p(x)
    std::complex<double> cIW(0.0, 2.0*3.14159265358979323846/(n0+1.0));
    std::vector<std::complex<double>> vcQLtemp(m0-1);
    std::complex<double> cTemp(0.0,0.0);
    std::complex<double> cQL(1.0,0.0);
    
    fComplete = false;
    while(!fComplete)
    {
      for(i=0; i<m0-1; i++) 
      {
        vcQLtemp[i] = exp(cIW * vdL[i]);
      }
      
      // iQL = prod(P[,m0] + P[,1:(m0-1)] * iQLtemp1);
      cQL.real(1.0); cQL.imag(0.0); // 1+0i
      for(i=0; i<n0; i++)
      {
        cTemp.real(vdP[i*m0+(m0-1)]); cTemp.imag(0.0); // P[i,m0] + 0i
        for(j=0; j<m0-1; j++)
        {
          cTemp += vdP[i*m0+j] * vcQLtemp[j];
        }
        cQL *= cTemp;
      }
      
      // dSum += get_real(iQL * exp(-1i*w * (vL' * vX[1:(m0-1)]))); 
      cTemp.real(0.0); cTemp.imag(0.0); // 0+0i
      for(i=0; i<m0-1; i++)
      {
        cTemp += vdL[i] * vdX[i];
      }
      cTemp = cQL * exp(-cIW*cTemp);
      dSum += cTemp.real();
      
      // generate next l: (0,0,0), (0,0,1), ..., (n0,n0,n0)
      j = m0-2;
      fCarry = true;
      while(fCarry)
      {
        vdL[j] = vdL[j] + 1.0;
        if(vdL[j] > n0)
        { // reset digit j to 0, increment digit j-1
          vdL[j] = 0.0;
          j = j-1;
          if(j<0) // signals completion
          {
            fCarry = false;
            fComplete = true;
          }
        } 
        else
        {
          fCarry = false;
        }
      }
    }
    
    if (dSum < 0.0) {
      Rcout << "⚠️ dSum = " << dSum << " — clipped to small epsilon to avoid log error." << std::endl;
      dSum = 1e-9;
    }
    
    // need to include the constant term (m0-1)log(n0+1)
    //    since logCL calcs mix direct, Heaps, and DFT
    dReturnVal = dNumLogSum - dLogK - std::log(dSum) + (m0-1.0) * std::log(n0+1);
  }
  return dReturnVal;
} // end logCL()


// R interface to logCL()
//   useful for testing
// [[Rcpp::export]]
double logCL_R(NumericVector nvLambda,
               int n0, // number of officers
               NumericVector nvS,
               IntegerVector ivY,
               int iHeapsLimit)
{
  std::vector<double> vLambda = Rcpp::as<std::vector<double>>(nvLambda);
  std::vector<double> vS      = Rcpp::as<std::vector<double>>(nvS);
  std::vector<int>    vY      = Rcpp::as<std::vector<int>>(ivY);
  
  return logCL(vLambda, n0, vS, vY, iHeapsLimit);
}


struct LogCLWorker : public Worker {
  
  const RVector<double> nvLambdaAll;
  const RVector<double> nvS;
  const RVector<int> ivIDOff;
  const RVector<int> ivY;
  const RVector<int> ivStart;
  const RVector<int> ivNOff;
  
  std::vector<double> vdStemp;

  double value;
  
  LogCLWorker(const NumericVector nvLambdaAllR, 
              const NumericVector nvSR, 
              const IntegerVector ivIDOffR, 
              const IntegerVector ivYR, 
              const IntegerVector ivStartR, 
              const IntegerVector ivNOffR)
    : nvLambdaAll(nvLambdaAllR),
      nvS(nvSR),
      ivIDOff(ivIDOffR),
      ivY(ivYR),
      ivStart(ivStartR),
      ivNOff(ivNOffR),
      value(0.0) 
  {
    vdStemp.resize(nvS.length());
    std::copy(nvS.begin(), nvS.end(), vdStemp.begin());
  }
  
  LogCLWorker(const LogCLWorker& lclw, Split)
    : nvLambdaAll(lclw.nvLambdaAll), 
      nvS(lclw.nvS), 
      ivIDOff(lclw.ivIDOff), 
      ivY(lclw.ivY), 
      ivStart(lclw.ivStart), 
      ivNOff(lclw.ivNOff), 
      value(0.0) 
  {
    vdStemp.resize(lclw.nvS.length());
    std::copy(lclw.nvS.begin(), lclw.nvS.end(), vdStemp.begin());
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    double dLocalSum = 0.0;
    std::vector<double> vdLambda;
    std::vector<int>    viYTemp;
    
    for (std::size_t i = begin; i < end; ++i) {
      int nOff = ivNOff[i];
      
      vdLambda.resize(nOff);
      viYTemp.resize(nOff);
      
      for (int j = 0; j < nOff; j++) {
        vdLambda[j] = nvLambdaAll[ivIDOff[ivStart[i] + j]];
        viYTemp[j]   = ivY[ivStart[i] + j];
      }
      
      dLocalSum += logCL(vdLambda, nOff, 
                          vdStemp, viYTemp,
                          7); // max/optimal for Heaps
    }
    value += dLocalSum;
  }
  
  // Join the results from different threads
  void join(const LogCLWorker& rhs) {
    value += rhs.value;
  }
};



// [[Rcpp::export]]
// compute logCL for each incident in parallel
double logCLfull(NumericVector nvLambdaAll,
                 NumericVector nvS,
                 IntegerVector ivIDOff,
                 IntegerVector ivY,
                 IntegerVector ivStart,
                 IntegerVector ivNOff)
{
  LogCLWorker worker(nvLambdaAll, nvS, ivIDOff, ivY, ivStart, ivNOff);
  
  parallelReduce(0, ivStart.size(), worker);
  
  return worker.value;
}
