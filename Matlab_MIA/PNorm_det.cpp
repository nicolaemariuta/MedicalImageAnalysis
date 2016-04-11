//**************************************************************
//PNorm.cpp is authored by Sune Darkner and is provided as is 
//The file may not be re-distributed without prior permission. 
//The auther cannot in any circumstances be held reliable of 
//damages, financial or other types of loses cased by the use 
//or the use of derivative products from this file. The author 
//cannot be held liable in any way under any circumstance
//if used at all the work must cite:
//Darkner, S & Sporring, J 2011, ' Generalized Partial Volume: 
//An Inferior Density Estimator to Parzen Windows for Normalized 
//Mutual Information ', i G Székely & HK Hahn (red) , 
//Information Processing in Medical Imaging - 22nd 
//International Conference, IPMI 2011: Lecture Notes in Computer 
//Science vol. 6801, Springer Publishing Company, s. 436-447. 
//**************************************************************
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
using namespace std;




void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
const mxArray *prhs[])
{

  if(nrhs==0)
			mexPrintf("PNorm takes 6 arguments,pts,ref_data, Image, range, nobins, offset, scale (optional) ");
		if(nrhs<6)
			mexErrMsgTxt("Number of arguments must be 6");
		//mexPrintf("Spline interpolation takes 4 arguments ");
		double* pts = static_cast<double*>(mxGetData(prhs[0]));
		double* dataR = static_cast<double*>(mxGetData(prhs[1]));
		double* dataI = static_cast<double*>(mxGetData(prhs[2]));
		double* range = static_cast<double*>(mxGetData(prhs[3]));
		double* no_bins = static_cast<double*>(mxGetData(prhs[4]));
		double* offset = static_cast<double*>(mxGetData(prhs[5]));
        double* det = static_cast<double*>(mxGetData(prhs[8]));
		double* scale=new double[2];
		scale[0]=1;
		scale[1]=1;

		scale = static_cast<double*>(mxGetData(prhs[6]));
		double* p = static_cast<double*>(mxGetData(prhs[7]));

		//double* det = static_cast<double*>(mxGetData(prhs[7]));
		
		int ndim =(int)mxGetNumberOfDimensions(prhs[2]);
		int* dim_pts =(int*)mxGetDimensions(prhs[0]);
		int* dim_img=(int*)mxGetDimensions(prhs[2]);
		int* dim_offset=(int*)mxGetDimensions(prhs[5]);
		int* dim_scale=(int*)mxGetDimensions(prhs[6]);
		int N=dim_pts[0];
		////Allocate Tc
		mwSize dims[2];
		mwSize dims2[1];

		dims2[0] = 1; 
		

		
		dims[0] = N; dims[1] = 1;
		
		double* diff1 ;
		double* diff2 ;
		double* diff3 ;
		bool do_deriv=false;
		if(nlhs>3)do_deriv=true;
		plhs[0] = mxCreateNumericArray(1,dims2,mxDOUBLE_CLASS, mxREAL);
		double* val = static_cast<double*>(mxGetData(plhs[0]));
		plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
		diff1 = static_cast<double*>(mxGetData(plhs[1]));
		plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
		diff2 = static_cast<double*>(mxGetData(plhs[2]));
		
		plhs[3] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
		diff3 = static_cast<double*>(mxGetData(plhs[3]));
		
	
	int s=64;
	int idx_idx=0;
	for(int i=0;i<N;i++){
		//val[i]=0;
		if(do_deriv){
		diff1[i]=0;
		diff2[i]=0;
		diff3[i]=0;
		}
		double t=(pts[i]-offset[0])/scale[0]-floor((pts[i]-offset[0])/scale[0]);
		double t2=t*t;
		double t3=t2*t;

		/*double* x={-pow(t,3)+3*pow(t,2)-3*t+1, 3*pow(t,3)-6*pow(t,2)+4, -3*pow(t,3)+3*pow(t,2)+3*t+1,  pow(t,3)};
		double* dx={-pow(t,3)+3*pow(t,2)-3*t+1, 3*pow(t,3)-6*pow(t,2)+4, -3*pow(t,3)+3*pow(t,2)+3*t+1,  pow(t,3)};
		*/
		double x[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
		double dx[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
		t=floor((pts[i]-offset[0])/scale[0])-1;
		int px[4]={(int)min<double>(max<double>(t-1,0),dim_img[0]-1),(int) max<double>(min<double>(t,dim_img[0]-1),0),(int)min<double>(max<double>(t+1,0),dim_img[0]-1),(int) max<double>(min<double>(t+2,dim_img[0]-1),0)};
		
		t=(pts[i+N]-offset[1])/scale[1]-floor((pts[i+N]-offset[1])/scale[1]);
		t2=t*t;
		t3=t2*t;

		double y[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
		double dy[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
		t=floor((pts[i+N]-offset[1])/scale[1])-1;
		int py[4]={(int)min<double>(max<double>(t-1,0),dim_img[1]-1),(int) max<double>(min<double>(t,dim_img[1]-1),0),(int)min<double>(max<double>(t+1,0),dim_img[1]-1),(int) max<double>(min<double>(t+2,dim_img[1]-1),0)};
		
		t=(pts[i+2*N]-offset[2])/scale[2]-floor((pts[i+2*N]-offset[2])/scale[2]);
		t2=t*t;
		t3=t2*t;

		double z[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
		double dz[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
		t=floor((pts[i+2*N]-offset[2])/scale[2])-1;
		int pz[4]={(int)min<double>(max<double>(t-1,0),dim_img[2]-1),(int) max<double>(min<double>(t,dim_img[2]-1),0),(int)min<double>(max<double>(t+1,0),dim_img[2]-1),(int) max<double>(min<double>(t+2,dim_img[2]-1),0)};
		int index=0;
      
        int idxR=(int)floor((dataR[i]-range[2]));	
			t=(dataR[i])-floor(dataR[i]);
			t2=t*t;
			t3=t2*t;

			double tr_val[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
			double lval=0;
		for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					for(int l=0;l<4;l++){

						int idx=px[l]+dim_img[0]*py[k]+dim_img[0]*dim_img[1]*pz[j];
						double dfdp=x[l]*y[k]*z[j]/216;
						lval+=dataI[idx]*dfdp;
						diff1[i]+=dx[l]*y[k]*z[j]/216*dataI[idx];
						diff2[i]+=x[l]*dy[k]*z[j]/216*dataI[idx];
						diff3[i]+=x[l]*y[k]*dz[j]/216*dataI[idx];
						
				
					}
				}
			}
			
			t=lval-floor(lval);
			t2=t*t;
			t3=t2*t;
			int dd=(int)(floor(lval)-range[0]);
			//			
			double t_val[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
			double dt_val[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
			double tmp;
			double dNMIdW=0;			
			for(int m=0;m<4;m++)
						{
								for (int nn=0;nn<4;nn++){
									tmp=pow(abs((idxR+nn-1)-(dd+m-1)),p[0])/(N*36)*tr_val[nn];
									val[0]+=t_val[m]*tmp*det[i];
									dNMIdW+=dt_val[m]*tmp;
								}
						}
			diff1[i]=dNMIdW*diff1[i];
			diff2[i]=dNMIdW*diff2[i];
			diff3[i]=dNMIdW*diff3[i];
            }
};
	
