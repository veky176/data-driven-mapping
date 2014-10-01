#include "cholesky.h"
#include "utility.h"

#define SAVE_RETMATRIX 0

const double epsilon = 0.00001; /** error tolerance **/

/** check if a matrix is an identity matrix **/
bool isEyeMatrix (double * matrix, int size, double eps)
{
	bool chkret = true;
	for (int row = 0; (row < size) && chkret; row ++)
		for (int col = 0; (col < size) && chkret; col ++)
			if (row == col)
				chkret = (abs (RowMajorElem(matrix, row, col, size) - 1.0) < eps);
			else
				chkret = (abs (RowMajorElem (matrix, row, col, size)) < eps);
	return chkret;
} /** isEyeMatrix **/

/** check if two matrices are equal to each other **/
bool isEqual (double * mat1, double * mat2, int size, double eps)
{
	bool chkret = true;
	int total = size * size;
	for (int elem = 0; (elem < total) && chkret; elem ++) {
			chkret = (abs (mat1[elem] - mat2 [elem]) < eps);
			if (! chkret) {
				cout << "Error: mat1 [" << elem << "] = " << mat1 [elem]
				<< " is not equal to mat2 [" << elem << "] = " << mat2 [elem]
				<< endl; cout.flush ();
			}
	}
	return chkret;
} /** isEqual **/


/** test A = LU where A is a square matrix **/
void test_lu (int n)
{
	SimpleTimer timer (__FUNCTION__);

	int matrix_size = n;
	std::stringstream infilename;
	std::stringstream outfilenameL;
	std::stringstream outfilenameU;

	infilename << "input_rand_" << matrix_size << ".txt";
	outfilenameL << "lu_retL_" << matrix_size << ".txt";
	outfilenameU << "lu_retU_" << matrix_size << ".txt";

	int matrix_elems = matrix_size * matrix_size;

	double * A = Cholesky::scalloc <double, int> (matrix_elems);
	double * L = Cholesky::scalloc <double, int> (matrix_elems);
	double * U = Cholesky::scalloc <double, int> (matrix_elems);

	readMatrixFromFile <double> (A, matrix_size, infilename.str().c_str());
	
#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") input matrix is "
		<< ((Cholesky::_is_symmetric <double, int> (A, matrix_size, true)) ? ("") : ("not "))
		<< "a symmetric matrix" << endl;
	cout.flush ();
#endif

	Cholesky::lu <double, int> (A, L, U, matrix_size);

#if SAVE_RETMATRIX
	SaveMatrix <double, int> (L, matrix_size, outfilenameL.str().c_str ()).save ();
	SaveMatrix <double, int> (U, matrix_size, outfilenameU.str().c_str ()).save ();
#endif

#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") resultant matrix L is "
		<< ((Cholesky::_is_lower_triangular <double, int> (L, matrix_size, true)) ? ("") : ("not "))
		<< "a lower triangular matrix" << endl;
	cout.flush ();

	cout << "(" << ProcInfo::GetThreadId () << ") resultant matrix U is "
		<< ((Cholesky::_is_upper_triangular <double, int> (U, matrix_size, true)) ? ("") : ("not "))
		<< "a upper triangular matrix" << endl;
	cout.flush ();
#endif

	/** check the resultant R **/
	bool chkret = true; /** check result **/
	int elem;
	double * retA = Cholesky::scalloc <double, int> (matrix_elems);

	Cholesky::multi2 <double, int> (L, U, retA, matrix_size, true);

	chkret = isEqual (A, retA, matrix_size, epsilon);

	Cholesky::sfree <double, int> (L, matrix_elems);
	Cholesky::sfree <double, int> (U, matrix_elems);
	Cholesky::sfree <double, int> (A, matrix_elems);
	Cholesky::sfree <double, int> (retA, matrix_elems);

	cout << "(" << ProcInfo::GetThreadId () << ")  ===> " << __FUNCTION__
		<< " <=== is " << ((chkret) ? ("passed") : ("failed")) << endl;
	cout.flush ();
} /** test_lu **/


/** test R = A / trans (B) where B is symmetric **/
void test_trsm (int n)
{
	SimpleTimer timer (__FUNCTION__);

	int matrix_size = n;
	std::stringstream infilename1;
	std::stringstream infilename2;
	std::stringstream outfilename;

	infilename1 << "input_" << matrix_size << ".txt";
	infilename2 << "input2_" << matrix_size << ".txt";
	outfilename << "trsm_" << matrix_size << "_ret.txt";

	int matrix_elems = matrix_size * matrix_size;

	double * A = Cholesky::scalloc <double, int> (matrix_elems);
	double * B = Cholesky::scalloc <double, int> (matrix_elems);
	double * R = Cholesky::scalloc <double, int> (matrix_elems);

	readMatrixFromFile <double> (A, matrix_size, infilename1.str().c_str());
	readMatrixFromFile <double> (B, matrix_size, infilename2.str().c_str());

	Cholesky::trsm <double, int> (A, B, R, matrix_size, true);

#if SAVE_RETMATRIX
	SaveMatrix <double, int> (R, matrix_size, outfilename.str().c_str ()).save ();
#endif

	double * ret = Cholesky::scalloc <double, int> (matrix_elems);

	Cholesky :: trmm4 <double, int> (R, B, ret, matrix_size, true);

	cout << "(" << ProcInfo::GetThreadId () << ")  ===> " << __FUNCTION__
		<< " <=== is " << ((isEqual (A, ret, matrix_size, epsilon)) ? ("passed") : ("failed"))
		<< endl;
	cout.flush ();

	Cholesky::sfree <double, int> (ret, matrix_elems);
	Cholesky::sfree <double, int> (A, matrix_elems);
	Cholesky::sfree <double, int> (B, matrix_elems);
	Cholesky::sfree <double, int> (R, matrix_elems);
} /** test_trsm **/


/** test Cholesky::inverse_lo_triangular : R = inv (L) **/
void test_inverse_lo_triangular ()
{
	SimpleTimer timer (__FUNCTION__);

	int matrix_size = 500;
	std::string infilename ("matrix_LoTri_500.txt");
	std::string outfilename ("matrix_LoTri_ret.txt");

	int matrix_elems = matrix_size * matrix_size;

	double * L = Cholesky::scalloc <double, int> (matrix_elems);
	double * R = Cholesky::scalloc <double, int> (matrix_elems);

	readMatrixFromFile <double> (L, matrix_size, infilename.c_str());
	
#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") input matrix is "
		<< ((Cholesky::_is_lower_triangular <double, int> (L, matrix_size, true)) ? ("") : ("not "))
		<< "a lower triangular matrix (" << matrix_size << " x " << matrix_size << ")" << endl;
	cout.flush ();
#endif

	Cholesky::inverse_lo_triangular <double, int> (L, R, matrix_size, true);

#if SAVE_RETMATRIX
	SaveMatrix <double, int> (R, matrix_size, outfilename.c_str ()).save ();
#endif

#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") resultant matrix is "
		<< ((Cholesky::_is_lower_triangular <double, int> (R, matrix_size, true)) ? ("") : ("not "))
		<< "a lower triangular matrix" << endl;
	cout.flush ();
#endif

	/** check the inverse matrix **/
	double * ret = Cholesky::scalloc <double, int> (matrix_elems);
	int row, col;
	bool chkret = true;

	Cholesky::multi2 <double, int> (L, R, ret, matrix_size);

	cout << "(" << ProcInfo::GetThreadId () << ")  ===> " << __FUNCTION__
		<< " <=== is " << ((isEyeMatrix (ret, matrix_size, epsilon)) ? ("passed") : ("failed"))
		<< endl;
	cout.flush ();
	
	Cholesky::sfree <double, int> (ret, matrix_elems);
	Cholesky::sfree <double, int> (R, matrix_elems);
	Cholesky::sfree <double, int> (L, matrix_elems);
} /** test_inverse_lo_triangular **/


/** test Cholesky::inverse_up_triangular **/
void test_inverse_up_triangular ()
{
	SimpleTimer timer (__FUNCTION__);

	int matrix_size = 500;
	std::string infilename ("matrix_UpTri_500.txt");
	std::string outfilename ("matrix_UpTri_ret.txt");

	int matrix_elems = matrix_size * matrix_size;

	double * U = Cholesky::scalloc <double, int> (matrix_elems);
	double * R = Cholesky::scalloc <double, int> (matrix_elems);

	readMatrixFromFile <double> (U, matrix_size, infilename.c_str());

#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") input matrix is "
		<< ((Cholesky::_is_upper_triangular <double, int> (U, matrix_size, true)) ? ("") : ("not "))
		<< "a upper triangular matrix (" << matrix_size << " x " << matrix_size << ")" << endl;
	cout.flush ();
#endif

	Cholesky::inverse_up_triangular <double, int> (U, R, matrix_size, true);

#if SAVE_RETMATRIX
	SaveMatrix <double, int> (R, matrix_size, outfilename.c_str ()).save ();
#endif

#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") resultant matrix is "
		<< ((Cholesky::_is_upper_triangular <double, int> (R, matrix_size, true)) ? ("") : ("not "))
		<< "a upper triangular matrix" << endl;
	cout.flush ();
#endif

	/** check the inverse matrix **/
	double * ret = Cholesky::scalloc <double, int> (matrix_elems);
	int row, col;
	bool chkret = true;

	Cholesky::multi2 <double, int> (U, R, ret, matrix_size);

	cout << "(" << ProcInfo::GetThreadId () << ")  ===> " << __FUNCTION__
		<< " <=== is " << ((isEyeMatrix (ret, matrix_size, epsilon)) ? ("passed") : ("failed"))
		<< endl;
	cout.flush ();

	Cholesky::sfree <double, int> (ret, matrix_elems);
	Cholesky::sfree <double, int> (R, matrix_elems);
	Cholesky::sfree <double, int> (U, matrix_elems);
} /** test_inverse_up_triangular **/


/** test the sequential Cholesky factorization algorithm : A = L * trans (L) **/
void test_potrf (int n)
{
	SimpleTimer timer (__FUNCTION__);

	int matrix_size = n;
	std::stringstream infilename;
	std::stringstream outfilename;

	infilename << "input_" << matrix_size << ".txt";
	outfilename << "chol_" << matrix_size << ".txt";

	int matrix_elems = matrix_size * matrix_size;

	double * A = Cholesky::scalloc <double, int> (matrix_elems);
	double * L = Cholesky::scalloc <double, int> (matrix_elems);

	readMatrixFromFile <double> (A, matrix_size, infilename.str().c_str());

#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") input matrix is "
		<< ((Cholesky::_is_symmetric <double, int> (A, matrix_size, true)) ? ("") : ("not "))
		<< "a symmetric matrix" << endl;
	cout.flush ();
#endif

	Cholesky::potrf <double, int> (A, L, matrix_size, true);

#if SAVE_RETMATRIX
	SaveMatrix <double, int> (L, matrix_size, outfilename.str().c_str ()).save ();
#endif

#ifdef CHECK_ASSUMPTIONS
	cout << "(" << ProcInfo::GetThreadId () << ") resultant matrix is "
		<< ((Cholesky::_is_lower_triangular <double, int> (L, matrix_size, true)) ? ("") : ("not "))
		<< "a lower triangular matrix" << endl;
	cout.flush ();
#endif

	/** check the resultant R **/
	bool chkret = true; /** check result **/
	int elem;
	double * trans_L = Cholesky::scalloc <double, int> (matrix_elems);
	double * retA = Cholesky::scalloc <double, int> (matrix_elems);
	Cholesky::transpose <double, int> (L, trans_L, matrix_size, false);
	Cholesky::multi2 <double, int> (L, trans_L, retA, matrix_size, true);
	for (elem = 0; elem < matrix_elems && chkret; elem ++)  
		chkret = (abs (A [elem] - retA [elem]) < epsilon);


	Cholesky::sfree <double, int> (trans_L, matrix_elems);
	Cholesky::sfree <double, int> (retA, matrix_elems);
	Cholesky::sfree <double, int> (L, matrix_elems);
	Cholesky::sfree <double, int> (A, matrix_elems);

	cout << "(" << ProcInfo::GetThreadId () << ")  ===> " << __FUNCTION__
		<< " <=== is " << ((chkret) ? ("passed") : ("failed")) << endl;
	cout.flush ();
} /** test_potrf **/


/** test R = inv (A) where A is a triangular matrix **/
void test_trtri (int lo_tri, int n)
{
	SimpleTimer timer (__FUNCTION__);

	int matrix_size = n;
	std::stringstream infilename;
	std::stringstream outfilename;

	if (lo_tri)  {
		infilename << "matrix_LoTri_" << matrix_size << ".txt";
		outfilename << "trtri_LoTri_" << matrix_size << "_ret.txt";
	} else {
		infilename << "matrix_UpTri_" << matrix_size << ".txt";
		outfilename << "trtri_UpTri_" << matrix_size << "_ret.txt";
	}
	
	int matrix_elems = matrix_size * matrix_size;

	double * A = Cholesky::scalloc <double, int> (matrix_elems);
	double * R = Cholesky::scalloc <double, int> (matrix_elems);

	readMatrixFromFile <double> (A, matrix_size, infilename.str().c_str());

	Cholesky::trtri <double, int> (A, R, matrix_size, lo_tri, true);

#if SAVE_RETMATRIX
	SaveMatrix <double, int> (R, matrix_size, outfilename.str().c_str ()).save ();
#endif

	double * ret = Cholesky::scalloc <double, int> (matrix_elems);

	Cholesky::multi2 <double, int> (A, R, ret, matrix_size, true);
	
	cout << "(" << ProcInfo::GetThreadId () << ")  ===> " << __FUNCTION__
		<< " <=== is " << ((isEyeMatrix (ret, matrix_size, epsilon)) ? ("passed") : ("failed"))
		<< endl;
	cout.flush ();

	Cholesky::sfree <double, int> (ret, matrix_elems);
	Cholesky::sfree <double, int> (A, matrix_elems);
	Cholesky::sfree <double, int> (R, matrix_elems);
} /** test_trtri **/