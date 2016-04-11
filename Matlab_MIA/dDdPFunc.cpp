#include <math.h>
#include <matrix.h>
#include <mex.h>

class interp {

private:
  double* d ;
  double* dfdp ;
  double*	dDdP;
  int* idx ;
  int N;
  int Np;
 
public:
  interp(double* vd,double* vdfdp,double* vdDdP,int* vidx, int vN, int vNp) {
    d=vd;
    dfdp=vdfdp;
    dDdP=vdDdP;
    idx=vidx;
    N=vN;
    Np=vNp;  
  }
  
  void operator()(int r1, int r2) const {
    
    double* pdDdp=new double[(int)3*Np];
    
    int s=64;
    for(int i=0;i<3*Np;i++)pdDdp[i]=0;
    int j;

    for(int i=r1;i!=r2;i++) {
      j=i>>6;
		
		//mexprintf("%f",j);
      
      pdDdp[idx[i]]+=dfdp[i]*d[j];
      pdDdp[idx[i]+Np]+=dfdp[i]*d[j+N];
      pdDdp[idx[i]+2*Np]+=dfdp[i]*d[j+2*N];    
    }
    
    for(int j=0;j<3*Np;j++)
      dDdP[j]=dDdP[j]+pdDdp[j];
    
    delete pdDdp;
  }
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  double* d = static_cast<double*>(mxGetData(prhs[0]));
  double* dfdp = static_cast<double*>(mxGetData(prhs[1]));
  int* idx = static_cast<int*>(mxGetData(prhs[2]));
  double* dim_p = static_cast<double*>(mxGetData(prhs[3]));
  
  int* dim_d =(int*)mxGetDimensions(prhs[0]);
  int N=dim_d[0];
  int Np=(int)dim_p[0];

  ////Allocate Tc
  mwSize dims[2];
  dims[0] = (int)dim_p[0]; dims[1] = 3;
  
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  double* dDdP = static_cast<double*>(mxGetData(plhs[0]));
  for(int i=0;i<3*Np;i++)
    dDdP[i]=0;
  
  interp interpO = interp(d,dfdp,dDdP,idx,N,Np);
  interpO(0,N*64);
  
  return; 
};
