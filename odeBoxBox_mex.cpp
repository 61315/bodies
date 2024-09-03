// To compile on my machine:
// mex COPTIMFLAGS='-O3 -DNDEBUG' -I/Users/sueda/lib/eigen-3.4.0 odeBoxBox_mex.cpp odeBoxBox.cpp

#include <vector>
#include "mex.h"
#include "odeBoxBox.h"

using namespace Eigen;

#define A(i,j) A[i+j*M]

void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
	// check for proper number of arguments
	if(nrhs != 4) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","4 inputs required.");
	}
	if(nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
	}

	enum
	{
		RHS_E1 = 0,
		RHS_WHD1,
		RHS_E2,
		RHS_WHD2,
		RHS_COUNT
	};

	if(!mxIsDouble(prhs[RHS_E1]) || mxGetNumberOfElements(prhs[RHS_E1]) != 16) {
		mexErrMsgTxt("E1 must be a mat4.");
	}
	if(!mxIsDouble(prhs[RHS_WHD1]) || mxGetNumberOfElements(prhs[RHS_WHD1]) != 3) {
		mexErrMsgTxt("whd1 must be a vec3.");
	}
	if(!mxIsDouble(prhs[RHS_E2]) || mxGetNumberOfElements(prhs[RHS_E2]) != 16) {
		mexErrMsgTxt("E1 must be a mat4.");
	}
	if(!mxIsDouble(prhs[RHS_WHD2]) || mxGetNumberOfElements(prhs[RHS_WHD2]) != 3) {
		mexErrMsgTxt("whd1 must be a vec3.");
	}

	mwSize M;
	double *A;

	// Convert E1
	M = mxGetM(prhs[RHS_E1]); A = mxGetPr(prhs[RHS_E1]);
	Matrix4d E1;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			E1(i,j) = A(i,j);
		}
	}

	// Convert whd1
	M = mxGetM(prhs[RHS_WHD1]); A = mxGetPr(prhs[RHS_WHD1]);
	Vector3d whd1;
	for(int i = 0; i < 3; ++i) {
		whd1(i) = A(i,0);
	}

	// Convert E2
	M = mxGetM(prhs[RHS_E2]); A = mxGetPr(prhs[RHS_E2]);
	Matrix4d E2;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			E2(i,j) = A(i,j);
		}
	}

	// Convert whd2
	M = mxGetM(prhs[RHS_WHD2]); A = mxGetPr(prhs[RHS_WHD2]);
	Vector3d whd2;
	for(int i = 0; i < 3; ++i) {
		whd2(i) = A(i,0);
	}

	// Call collision detector
	Contacts contacts = odeBoxBox(E1, whd1, E2, whd2);

	// Convert collisions
	int nfields = 5;
	const char **fnames;
	fnames = (const char **)mxCalloc(nfields, sizeof(*fnames));
	fnames[0] = "count";
	fnames[1] = "depthMax";
	fnames[2] = "nor";
	fnames[3] = "pos";
	fnames[4] = "depth";
	mxArray *c = mxCreateStructMatrix(1, 1, nfields, fnames);
	mxSetField(c, 0, "count", mxCreateDoubleScalar(contacts.count));
	mxSetField(c, 0, "depthMax", mxCreateDoubleScalar(contacts.depthMax));
	mxArray *mat;
	mat = mxCreateDoubleMatrix(3, 1, mxREAL); M = mxGetM(mat); A = mxGetPr(mat);
	for(int i = 0; i < 3; ++i) {
		A(i,0) = contacts.normal(i);
	}
	mxSetField(c, 0, "nor", mat);
	mat = mxCreateDoubleMatrix(3, contacts.count, mxREAL); M = mxGetM(mat); A = mxGetPr(mat);
	for(int k = 0; k < contacts.count; ++k) {
		for(int i = 0; i < 3; ++i) {
			A(i,k) = contacts.positions[k](i);
		}
	}
	mxSetField(c, 0, "pos", mat);
	mat = mxCreateDoubleMatrix(1, contacts.count, mxREAL); M = mxGetM(mat); A = mxGetPr(mat);
	for(int k = 0; k < contacts.count; ++k) {
		A(0,k) = contacts.depths[k];
	}
	mxSetField(c, 0, "depth", mat);
	mxFree((void *)fnames);
	plhs[0] = c;
}
