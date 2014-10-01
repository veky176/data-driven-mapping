#ifndef _CHOLESKY_H
#define _CHOLESKY_H

#include "../strassen/strassen.h"

// #define CHECK_ASSUMPTIONS

/** test functions **/
void test_lu (int n = 1000); /** test A = LU where A is a square matrix **/
void test_trsm (int n = 1000);  /** test R = A / trans (B) where B is symmetric **/
void test_inverse_lo_triangular (); /** test R = inv (L) where L is a lower triangular matrix **/
void test_inverse_up_triangular (); /** test R = inv (U) where U is a upper triangular matrix **/
void test_potrf (int n = 1000);  /** test R = chol (A) where is A is SPD **/
void test_trtri (int lo_tri = true, int n = 500); /** test R = inv (A) where A is a triangular matrix **/

class Cholesky : public Strassen
{
public:

	Cholesky () {}
	virtual ~Cholesky () {}

	/********************/

	template <typename DataType, typename SizeType>
	static bool _is_symmetric (DataType * A, SizeType size, bool row_major = true);

	template <typename DataType, typename SizeType>
	static bool _is_lower_triangular (DataType * A, SizeType size, bool row_major = true);

	template <typename DataType, typename SizeType>
	static bool _is_upper_triangular (DataType * A, SizeType size, bool row_major = true);

#ifdef CHECK_ASSUMPTIONS
	template <typename DataType, typename SizeType>
	static bool _is_spd_matrix (DataType * A, SizeType size, bool row_major = true);

	template <typename DataType, typename SizeType>
	static bool _is_invertible (DataType * A, SizeType size);

	/** return the determinate of a square matrix (size-by-size) (todo-list) **/
	template <typename DataType, typename SizeType>
	static DataType _det_matrix (DataType * A, SizeType size);
#endif
	/********************/

	/** Here, we always assume that A is a row-major and symetric matrix **/

	/** partition a matrix into a set of submatrices **/
	template <typename DataType, typename SizeType>
	static void partition (DataType * A, SizeType size, int partition_size, int start_id,
		DataModuleContainer <DataType> & buffer);

	/** copy back submatrices to a matrix (symmetric) **/
	template <typename DataType, typename SizeType>
	static void copyback (DataType * A, SizeType size, int partition_size, int start_id,
		DataModuleContainer <DataType> & buffer);

	/** copy back submatrices to a matrix (lower triangular) **/
	template <typename DataType, typename SizeType>
	static void copyback2 (DataType * A, SizeType size, int partition_size, int start_id,
		DataModuleContainer <DataType> & buffer);

	/** transpose a matrix : B = transpose of A **/
	template <typename DataType, typename SizeType>
	static void transpose (DataType * A, DataType * B, SizeType size, bool symmetry = true);

	/** copy a matrix : A to B **/
	template <typename DataType, typename SizeType>
	static void copy (DataType * A, DataType * B, SizeType size);

	/** assume: A/B/C are both raw-major matrices, C = A * B (refresh, default) or C+ = A * B (not refresh) **/
	template <typename DataType, typename SizeType>
	static void multi2 (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh = true);

	/** R = inv (L) and L is a lower triangular matrix, R should be also a lower triangular matrix **/
	template <typename DataType, typename SizeType>
	static void inverse_lo_triangular (DataType * L, DataType * R, SizeType size, bool refresh = true);

	/** R = inv (U) and U is a upper triangular matrix, R should be also a upper triangular matrix **/
	template <typename DataType, typename SizeType>
	static void inverse_up_triangular (DataType * U, DataType * R, SizeType size, bool refresh = true);

	/** LU decomposition : A = LU **/
	template <typename DataType, typename SizeType>
	static void lu (DataType * A, DataType * L, DataType * U, SizeType size);

	/**
	Essential algorithms for Cholesky are as follows.
	In the following algorithms, we always assume A/B/C/R are in row-major (no exception).
	**/

	/** POTRF : compute the Cholesky factorization (lower triangular)
	i.e., it can compute lower triangular A such that A = R * trans (R) **/
	template <typename DataType, typename SizeType>
	static void potrf (DataType * A, DataType * R, SizeType size, bool refresh = true);

	/** TRSM : solve the matrix equation : R * trans (B) = A or R = A / trans(B) **/
	template <typename DataType, typename SizeType>
	static void trsm (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** TRSM2 : solve the matrix equation : R * B = A or R = A / B **/
	template <typename DataType, typename SizeType>
	static void trsm2 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** TRSM3 : solve the matrix equation : R = - A / B **/
	template <typename DataType, typename SizeType>
	static void trsm3 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** SYRK : R = A - B * trans (B) **/
	template <typename DataType, typename SizeType>
	static void syrk (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** SYRK2 : R = A + trans (B) * B **/
	template <typename DataType, typename SizeType>
	static void syrk2 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** GEMM : R = A - B * trans (C) **/
	template <typename DataType, typename SizeType>
	static void gemm (DataType * A, DataType * B, DataType * C, DataType * R, SizeType size, bool refresh = true);

	/** GEMM2 : R = A - B * C **/
	template <typename DataType, typename SizeType>
	static void gemm2 (DataType * A, DataType * B, DataType * C, DataType * R, SizeType size, bool refresh = true);

	/** GEMM3 : R = A + trans (B) * C **/
	template <typename DataType, typename SizeType>
	static void gemm3 (DataType * A, DataType * B, DataType * C, DataType * R, SizeType size, bool refresh = true);

	/** TRTRI : compute the inverse of a upper or lower triangular matrix A
	R = inv (A) where A has to be a upper or lower triangular matrix,
	default : A is a lower triangular matrix (lo_triangular = true) **/
	template <typename DataType, typename SizeType>
	static void trtri (DataType * A, DataType * R, SizeType size, bool lo_triangular = true, bool refresh = true);

	/** LAUUM : R = trans (A) * A **/
	template <typename DataType, typename SizeType>
	static void lauum (DataType * A, DataType * R, SizeType size, bool refresh = true);

	/** TRMM : R = A * B **/
	template <typename DataType, typename SizeType>
	static void trmm (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** TRMM2 : R = - A * B **/
	template <typename DataType, typename SizeType>
	static void trmm2 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** TRMM3 : R = trans (A) * B **/
	template <typename DataType, typename SizeType>
	static void trmm3 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

	/** TRMM4 : R = A * trans (B) for Cholesky :: trsm (only) **/
	template <typename DataType, typename SizeType>
	static void trmm4 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh = true);

}; /** class Cholesky **/


/** implementation **/
/** partition a matrix into a set of submatrices : A is assumed to be a row-major **/
template <typename DataType, typename SizeType>
static void Cholesky :: partition (DataType * A, SizeType size,
	int partition_size, int start_id, DataModuleContainer <DataType> & buffer) 
{
	SizeType block_size = static_cast <SizeType> (size / partition_size);
	bool row_major = true;

	if ( (size % partition_size ) != 0) {
		cout << "Error : " << size << " must be a multiple of " << partition_size << endl;
		cout.flush (); 
		press_key ();
		return;
	} else {
		cout << "(" << ProcInfo::GetThreadId () << ") data partition "
			<< size << " x " << size << " matrix (ptr = " << A << ") "
			<< ((row_major) ? "row_major" : "column_major")
			<< " to " << partition_size << " x " << partition_size << " submatrices (each = "
			<< block_size << " x " << block_size << ") start_id = " << start_id
			<< endl; cout.flush ();
	}

#ifdef CHECK_ASSUMPTIONS
	if (! Cholesky::_is_symmetric <DataType, SizeType> (A, size, true)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Input Error: (ptr = " << A << ") should be a symmetric matrix so "
			<< __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str (""); 
		press_key ();
		return;
	}
#endif

	DataModule <DataType> * block = NULL;
	int id = start_id;

	for (int i = 0; i < partition_size; i++) { /** col **/
		for (int j = i; j < partition_size; j++) { /** row **/
			block = new DataModule <DataType> (id, block_size, -1);
			if (block == NULL) {
				cout << "Error: memory not enough" << endl; cout.flush ();
				press_key ();
				return;
			}
			DataType * ptr = block->data ();
			/** data copy **/
			for (int row = 0; row < block_size; row ++) {
				for (int col = 0; col < block_size; col ++) {
					RowMajorElem (ptr, row, col, block_size) =
					RowMajorElem (A, row + j * block_size, col + i * block_size, size);
				}
			}
			buffer.add (block);
			id += 1;
		} /** end of j-loop **/
	} /** end of i-loop **/

} /** Cholesky :: partition **/


/** copy back submatrices to a matrix **/
template <typename DataType, typename SizeType>
static void Cholesky :: copyback (DataType * A, SizeType size, int partition_size, int start_id,
	DataModuleContainer <DataType> & buffer)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType block_size = static_cast <SizeType> (size / partition_size);
	bool row_major = true;

	if ( (size % partition_size ) != 0) {
		cout << "Error : " << size << " must be a multiple of " << partition_size << endl;
		cout.flush ();
		press_key ();
		return;
	} else {
		cout << "(" << ProcInfo::GetThreadId () << ") data copyback "
			<< partition_size << " x " << partition_size 
			<< " submatrices (each = " << block_size << " x " << block_size << ") to "
			<< size << " x " << size << " matrix (ptr = " << A << ") "
			<< ((row_major) ? "row_major" : "column_major")
			<< " start_id = " << start_id << endl;
		cout.flush ();
	}

	DataModule <DataType> * block = NULL;
	int id = start_id, row, col;

	for (int i = 0; i < partition_size; i++) { /** col **/
		for (int j = 0; j < partition_size; j++) { /** row **/
			id = start_id + max (i, j) - min (i, j);
			for (int k = 0; k < min (i, j); k++) id += (partition_size  - k);
			block = buffer [id];
			if (block == NULL) {
				cout << "Error: cannot find the data module, id = " << id << endl;
				cout.flush ();
				press_key ();
				return;
			}
			DataType * ptr = block->data ();
			bool trans = (j < i);
			/** data copy **/
			cout << "(" << ProcInfo::GetThreadId () << ") copy " << ((trans) ? "transpose " : "")
				<< "(data " << block->id () << ") to submatrix (" << j << ", " << i << ")" << endl;
			cout.flush ();

			if (trans) {
				for (int row = 0; row < block_size; row ++) {
					for (int col = 0; col < block_size; col ++) {
						RowMajorElem (A, row + j * block_size, col + i * block_size, size) = 
							RowMajorElem (ptr, col, row, block_size);
					}
				}
			} else {
				for (int row = 0; row < block_size; row ++) {
					for (int col = 0; col < block_size; col ++) {
						RowMajorElem (A, row + j * block_size, col + i * block_size, size) = 
							RowMajorElem (ptr, row, col, block_size);
					}
				}
			}
			/** data copy **/

		} /** end of j-loop **/
	} /** end of i-loop **/

	cout << "(" << ProcInfo::GetThreadId() << ") output matrix is "
		<< (Cholesky::_is_symmetric <DataType, SizeType> (A, size, true) ? "" : "not ")
		<< "a symmetric matrix" << endl;
	cout.flush ();

} /** Cholesky :: copyback **/


/** copy back submatrices to a matrix (lower triangular) **/
template <typename DataType, typename SizeType>
static void Cholesky :: copyback2 (DataType * A, SizeType size, int partition_size, int start_id,
		DataModuleContainer <DataType> & buffer)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType block_size = static_cast <SizeType> (size / partition_size);
	bool row_major = true;

	if ( (size % partition_size ) != 0) {
		cout << "Error : " << size << " must be a multiple of " << partition_size << endl;
		cout.flush ();
		press_key ();
		return;
	} else {
		cout << "(" << ProcInfo::GetThreadId () << ") data copyback "
			<< partition_size << " x " << partition_size 
			<< " submatrices (each = " << block_size << " x " << block_size << ") to "
			<< size << " x " << size << " lower-triangular matrix (ptr = " << A << ") "
			<< ((row_major) ? "row_major" : "column_major")
			<< " start_id = " << start_id << endl;
		cout.flush ();
	}

	DataModule <DataType> * block = NULL;
	int id = start_id, row, col;

	for (int i = 0; i < partition_size; i++) { /** col **/
		for (int j = 0; j < partition_size; j++) { /** row **/
			id = start_id + max (i, j) - min (i, j);
			for (int k = 0; k < min (i, j); k++) id += (partition_size  - k);
			block = buffer [id];
			if (block == NULL) {
				cout << "Error: cannot find the data module, id = " << id << endl;
				cout.flush ();
				press_key ();
				return;
			}

			if (j < i)  continue; /** only copy the lower triangular part **/

			DataType * ptr = block->data ();

			/** data copy **/
			cout << "(" << ProcInfo::GetThreadId () << ") copy (data "
				<< block->id () << ") to submatrix (" << j << ", " << i << ")" << endl;
			cout.flush ();

			for (int row = 0; row < block_size; row ++) {
				for (int col = 0; col < block_size; col ++) {
					RowMajorElem (A, row + j * block_size, col + i * block_size, size) =
						RowMajorElem (ptr, row, col, block_size);
				}
			} /** data copy **/

		} /** end of j-loop **/
	} /** end of i-loop **/

	cout << "(" << ProcInfo::GetThreadId() << ") output matrix is "
		<< (Cholesky::_is_lower_triangular <DataType, SizeType> (A, size, true) ? "" : "not ")
		<< "a lower triangular matrix" << endl;
	cout.flush ();

} /** Cholesky :: copyback2 **/


/** transpose a matrix : B = transpose of A **/
template <typename DataType, typename SizeType>
static void Cholesky :: transpose (DataType * A, DataType * B, SizeType size, bool symmetry)
{
	SizeType row, col;

#ifdef CHECK_ASSUMPTIONS
	symmetry = Cholesky::_is_symmetric <DataType, SizeType> (A, size, true);
#endif

	if (symmetry) { /** for a symmetric matrix **/
		Cholesky::copy <DataType, SizeType> (A, B, size);
	} 
	else { /** for a general dense matrix **/
		if (A == B) {
			for (row = 0; row < size; row ++) {
				for (col = 0; col < size; col ++) {
					DataType t = RowMajorElem (B, col, row, size);
					RowMajorElem (B, col, row, size) = RowMajorElem (A, row, col, size);
					RowMajorElem (A, row, col, size) = t;
				}
			}
		} /** if (A == B) **/
		else {
			for (row = 0; row < size; row ++) {
				for (col = 0; col < size; col++)
					RowMajorElem (B, col, row, size) = RowMajorElem (A, row, col, size);
			}
		} /** else for if (A == B) **/
	} /** for a general matrix **/

} /** Cholesky :: transpose **/


/** copy a matrix : A to B **/
template <typename DataType, typename SizeType>
static void Cholesky :: copy (DataType * A, DataType * B, SizeType size)
{
	if (A == B)  return;

	for (SizeType elem = 0; elem < size * size; elem ++)  B [elem] = A [elem];
} /** Cholesky :: copy **/


/** assume: A/B/C are both raw-major matrices, 
	C = A * B (refresh, default) or C+ = A * B (not refresh) **/
template <typename DataType, typename SizeType>
static void Cholesky :: multi2 (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType r, c, k, c1, elem;
	DataType * ptr1, * ptr2, ret;

	for (r = 0; r < size; r++) {
		for (c = 0; c < size; c++) 
		{
			for (ret = 0, k = 0; k < size; k ++)
				ret += RowMajorElem (A, r, k, size) * RowMajorElem (B, k, c, size);
			if (refresh)
				RowMajorElem (C, r, c, size) = ret;
			else
				RowMajorElem (C, r, c, size) += ret;
		}
	}

#if 0
	if (refresh) {
		for (r = 0; r < size; r++) {
			ptr1 = A + r * size;
			for (c = 0; c < size; c++) {
				ptr2 = B + c;
				for (ret = 0, elem = 0, c1 = 0; elem < size; elem++, c1 += size)
					ret += ptr1 [elem] * ptr2 [c1];
				C [r * size + c] = ret; /** refresh **/
			}
		}
	}
	else {
		for (r = 0; r < size; r++) {
			ptr1 = A + r * size;
			for (c = 0; c < size; c++) {
				ptr2 = B + c;
				for (ret = 0, elem = 0, c1 = 0; elem < size; elem++, c1 += size)
					ret += ptr1 [elem] * ptr2 [c1];
				C [r * size + c] += ret; /** not refresh **/
			}
		}
	}
#endif

} /** Cholesky :: multi2 **/


/** check if a matrix is symmetric **/
template <typename DataType, typename SizeType>
static bool Cholesky :: _is_symmetric (DataType * A, SizeType size, bool row_major)
{
	DataType val1, val2;
	if (row_major) {
		for (SizeType row = 0; row < size; row ++) {
			for (SizeType col = row + 1; col < size; col ++) {
				val1 = RowMajorElem (A, row, col, size);
				val2 = RowMajorElem (A, col, row, size);
				if (val1 != val2)  return false;
			}
		}
	} else {
		for (SizeType row = 0; row < size; row ++) {
			for (SizeType col = row + 1; col < size; col ++) {
				val1 = ColMajorElem (A, row, col, size);
				val2 = ColMajorElem (A, col, row, size);
				if (val1 != val2)  return false;
			}
		}
	}
	return true;
} /** Cholesky :: _is_symmetric **/


/** check if a matrix is a lower triangular matrix **/
template <typename DataType, typename SizeType>
static bool Cholesky :: _is_lower_triangular (DataType * A, SizeType size, bool row_major)
{
	DataType val;

	if (row_major) {
		for (SizeType row = 0; row < size; row ++) {
			for (SizeType col = row + 1; col < size; col ++) {
				val = RowMajorElem (A, row, col, size);
				if (abs (val) > epsilon)  return false;
			}
		}
	} else {
		for (SizeType row = 0; row < size; row ++) {
			for (SizeType col = row + 1; col < size; col ++) {
				val = ColMajorElem (A, row, col, size);
				if (abs (val) > epsilon)  return false;
			}
		}
	}

	return true;
} /** Cholesky :: _is_lower_triangular **/


/** check if a matrix is a upper triangular matrix **/
template <typename DataType, typename SizeType>
static bool Cholesky :: _is_upper_triangular (DataType * A, SizeType size, bool row_major)
{
	DataType val;

	if (row_major) {
		for (SizeType row = 0; row < size; row ++) {
			for (SizeType col = 0; col < row; col ++) {
				val = RowMajorElem (A, row, col, size);
				if (abs (val) > epsilon)  return false;
			}
		}
	} else {
		for (SizeType row = 0; row < size; row ++) {
			for (SizeType col = 0; col < row; col ++) {
				val = ColMajorElem (A, row, col, size);
				if (abs (val) > epsilon)  return false;
			}
		}
	}

	return true;
} /** Cholesky :: _is_upper_triangular **/


#ifdef CHECK_ASSUMPTIONS

/** check if a matrix is a SPD matrix (note: this is a sufficient condition but not necessary.) **/
template <typename DataType, typename SizeType>
static bool Cholesky :: _is_spd_matrix (DataType * A, SizeType size, bool row_major)
{
	bool chkret = false;
	DataType row_sum, col_sum, value;

	if (! Cholesky::_is_symmetric <DataType, SizeType> (A, size, row_major)) {
#if 1
		cout << "Error: it is not a symmetric matrix" << endl;  cout.flush ();
#endif
		return false;
	}

	for (SizeType elem = 0; elem < size; elem ++) {
		row_sum = 0;
		col_sum = 0;
		value = RowMajorElem (A, elem, elem, size);
		if (value < 0) {
#if 1
			cout << "Error: matrix (" << elem << ", " << elem << ") = " << value << endl; cout.flush ();
#endif
			return false;
		}
		if (row_major) {
			for (SizeType k = 0; k < size; k ++)  {
				if (elem == k) continue;
				row_sum += static_cast <DataType> (abs (RowMajorElem (A, elem, k, size)));
			}
			for (SizeType k = 0; k < size; k ++)  {
				if (elem == k) continue;
				col_sum += static_cast <DataType> (abs (RowMajorElem (A, k, elem, size)));
			}
		} else {
			for (SizeType k = 0; k < size; k ++)  {
				if (elem == k) continue;
				row_sum += static_cast <DataType> (abs (ColMajorElem (A, elem, k, size)));
			}
			for (SizeType k = 0; k < size; k ++)  {
				if (elem == k) continue;
				col_sum += static_cast <DataType> (abs (ColMajorElem (A, k, elem, size)));
			}
		}
		if (value < row_sum ||  value < col_sum)  {
#if 0
			if (value < row_sum)
				cout << "Error: matrix (" << elem << ", " << elem << ") = " << value << ", row_sum = " << row_sum << endl;
			if (value < col_sum)
				cout << "Error: matrix (" << elem << ", " << elem << ") = " << value << ", col_sum = " << col_sum << endl;
			cout.flush ();
#endif
			return false;
		}
		if (value > row_sum && value > col_sum)  chkret = true;
	}

	return chkret;
} /** Cholesky :: _is_spd_matrix **/


/** check if a matrix is inversible **/
template <typename DataType, typename SizeType>
static bool Cholesky :: _is_invertible (DataType * A, SizeType size)
{
	bool chkret = true;

	/** for a triangular matrix **/
	if (Cholesky::_is_lower_triangular <DataType, SizeType> (A, size) || 
		Cholesky::_is_upper_triangular <DataType, SizeType> (A, size))
	{
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") checking invertibility of a "
			<< ((Cholesky::_is_lower_triangular <DataType, SizeType> (A, size)) ? "lower" : "upper")
			<< " matrix" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");

		for (SizeType elem = 0; elem < size && chkret; elem ++) {
			DataType value = RowMajorElem (A, elem, elem, size);
			chkret = (abs (value) > epsilon); /** non-zero diagonal elements **/
			if (! chkret) {
				ss << "(" << ProcInfo::GetThreadId() << ") matrix (ptr = " << A
					<< ") has a too small diagonal element (" << elem 	<< ", " 
					<< elem << ") = [" << value << "]" << endl;
				cout << ss.str (); cout.flush (); ss.str ("");
			}
		}
		return chkret;
	}

	/** for a general matrix, LU decomposition **/
	SizeType matrix_elems = size * size;
	DataType * L = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * U = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

	Cholesky::lu <DataType, SizeType> (A, L, U, size); /** A = L * U **/

	for (SizeType elem = 0; elem < size && chkret; elem ++) 
	{
		DataType value1 = RowMajorElem (L, elem, elem, size);
		DataType value2 = RowMajorElem (U, elem, elem, size);
		chkret = (abs (value1) > epsilon); /** non-zero diagonal elements **/
		chkret = (abs (value2) > epsilon); /** non-zero diagonal elements **/
	}

	Cholesky::sfree <DataType, SizeType> (L, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (U, matrix_elems);

	return chkret;

#if 0 /** risky : overflown **/
	DataType det = Cholesky::_det_matrix <DataType, SizeType> (A, size);

	std::stringstream ss;
	ss << "(" << ProcInfo::GetThreadId() << ") determint of the matrix is " << det << endl;
	std::cout << ss.str (); cout.flush (); ss.str ("");

	return (abs (det) > epsilon);
#endif
} /** Cholesky :: _is_invertible **/


/** return the determinate of a square matrix (size-by-size) **/
template <typename DataType, typename SizeType>
static DataType Cholesky :: _det_matrix (DataType * A, SizeType size)
{
	DataType det = 1.0;
	SizeType matrix_elems= size * size;

	if (Cholesky::_is_lower_triangular <DataType, SizeType> (A, size) || 
		Cholesky::_is_upper_triangular <DataType, SizeType> (A, size))
	{
		for (SizeType elem = 0; elem < size; elem ++) {
			DataType value = RowMajorElem (A, elem, elem, size);
			det = det * value;
			if (abs (value) < epsilon) {
				std::stringstream ss;
				ss << "(" << ProcInfo::GetThreadId() << ") matrix (ptr = " << A
					<< ") has a too small diagonal element (" << elem
					<< ", " << elem << ") = [" << value << "]"
					<< endl;
				cout << ss.str (); cout.flush (); ss.str ("");
			}
		}
		return det;
	}
	
	DataType * L = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * U = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

	Cholesky::lu <DataType, SizeType> (A, L, U, size); /** A = L * U **/

	for (SizeType elem = 0; elem < size; elem ++) 
	{
		det = det * RowMajorElem (L, elem, elem, size) * RowMajorElem (U, elem, elem, size);
	}

	Cholesky::sfree <DataType, SizeType> (L, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (U, matrix_elems);

	return det;
}

#endif // #ifdef CHECK_ASSUMPTIONS

/** --------------------------------------------- **/

/** R = inv (L) and L is a lower triangular matrix, R should be also a lower triangular matrix **/
template <typename DataType, typename SizeType>
static void Cholesky :: inverse_lo_triangular (DataType * L, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	
#ifdef CHECK_ASSUMPTIONS
	/** L has to be a lower triangular matrix **/
	if (! Cholesky::_is_lower_triangular <DataType, SizeType> (L, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Input Error: input matrix (ptr = " << L 
			<< ") should be a lower triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

	/** forward substitution **/
	for (SizeType col = 0; col < size; col ++) {
		RowMajorElem (R, col, col, size) = 1.0 / (RowMajorElem (L, col, col, size));
		for (SizeType row = col + 1; row < size; row ++) {
			DataType value = 0;
			for (SizeType k = col; k < row; k++) 
				value += RowMajorElem (L, row, k, size) * RowMajorElem (R, k, col, size);
			value = (-1) * value / (RowMajorElem (L, row, row, size));
			RowMajorElem (R, row, col, size) = value;
		} /** row-loop **/
	} /** col-loop **/
	
#ifdef CHECK_ASSUMPTIONS
	/** R has to be a lower triangular matrix **/
	if (! Cholesky::_is_lower_triangular <DataType, SizeType> (R, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Output Error: output matrix (ptr = " << R
			<< ") should be a lower triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

} /** Cholesky :: inverse_lo_triangular **/

/** --------------------------------------------- **/

/** R = inv (U) and U is a upper triangular matrix, R should be also a upper triangular matrix **/
template <typename DataType, typename SizeType>
static void Cholesky :: inverse_up_triangular (DataType * U, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	
#ifdef CHECK_ASSUMPTIONS
	/** U has to be a upper triangular matrix **/
	if (! Cholesky::_is_upper_triangular <DataType, SizeType> (U, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Input Error: input matrix (ptr = " << U
			<< ") should be a upper triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

	/** backward substitution **/
	for (SizeType col = 0; col < size; col ++) {
		RowMajorElem (R, col, col, size) = 1.0 / (RowMajorElem (U, col, col, size));
		for (SizeType row = col + 1; row < size; row ++) {
			DataType value = 0;
			for (SizeType k = col; k < row; k++) 
				value += RowMajorElem (U, k, row, size) * RowMajorElem (R, col, k, size);
			value = (-1) * value / (RowMajorElem (U, row, row, size));
			RowMajorElem (R, col, row, size) = value;
		} /** row-loop **/
	} /** col-loop **/

#ifdef CHECK_ASSUMPTIONS
	/** R has to be a upper triangular matrix **/
	if (! Cholesky::_is_upper_triangular <DataType, SizeType> (R, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Output Error: output matrix (ptr = " << R
			<< ") should be a upper triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

	return;
} /** Cholesky :: inverse_up_triangular **/

/** --------------------------------------------- **/

/** LU decomposition : A = LU **/
template <typename DataType, typename SizeType>
static void Cholesky :: lu (DataType * A, DataType * L, DataType * U, SizeType size)
{
	SimpleTimer timer (__FUNCTION__);

	SizeType k, i, j, p;
	DataType sum = 0;

	/** LU decomposition :
	Assume that L and U are now two zero matrices.
	Diagonal elements of U are all one's. **/
	for (k = 0; k < size; k++) RowMajorElem (U, k, k, size) = 1;
	L [0] = A [0];
	if (abs (L[0]) > epsilon) {
		for (k = 1; k < size; k++)
			RowMajorElem (U, 0, k, size) = RowMajorElem (A, 0, k, size) / L [0];
	}

	for (k = 1; k < size; k ++) {

		for (i = 0; i <= k; i++) {
			sum = 0;
			for (p = 0; p < i; p++) 
				sum += RowMajorElem (L, k, p, size) * RowMajorElem (U, p, i, size);
			RowMajorElem (L, k, i, size) = RowMajorElem (A, k, i, size) - sum;
		}

		for (j = k + 1; j < size; j++) {
			sum = 0;
			for (p = 0; p < k; p++)
				sum += RowMajorElem (L, k, p, size) * RowMajorElem (U, p, j, size);
			sum = RowMajorElem (A, k, j, size) - sum;
			RowMajorElem (U, k, j, size) = sum / (RowMajorElem (L, k, k, size));
		}
	}

#ifdef CHECK_ASSUMPTIONS
	/** L has to be a upper triangular matrix **/
	if (! Cholesky::_is_lower_triangular <DataType, SizeType> (L, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Output Error: output matrix (ptr = " << L
			<< ") should be a lower triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
	/** U has to be a upper triangular matrix **/
	if (! Cholesky::_is_upper_triangular <DataType, SizeType> (U, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Output Error: output matrix (ptr = " << U
			<< ") should be a upper triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif
} /** Cholesky :: lu **/


/**************************************************/

/** essential algorithms for Cholesky (as below) **/

/** POTRF : compute the Cholesky factorization (lower triangular: R)
	i.e., it can compute lower triangular A such that A = R * trans (R) **/
template <typename DataType, typename SizeType>
static void Cholesky :: potrf (DataType * A, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	DataType ret, temp;
	SizeType row, col, k;

#ifdef CHECK_ASSUMPTIONS
	/** A has to be a SPD matrix **/
	if (! Cholesky::_is_symmetric <DataType, SizeType> (A, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Input Error: input matrix (ptr = " << A
			<< ") should be a symmetric matrix" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

	/** set zeros to R **/
	Cholesky::zeros <DataType, SizeType> (R, size);

	/** TODO :: chol (A) = R * trans (R) **/
	for (row = 0; row < size; row ++) {
		ret = 0.0;
		/** compute R (row, col) for row > col **/
		for (col = 0; col < row; col ++) 
		{
			temp = 1.0 / RowMajorElem (R, col, col, size);
			ret = 0.0;
			for (k = 0; k < col; k++) {
				DataType d1 = RowMajorElem (R, row, k, size);
				DataType d2 = RowMajorElem (R, col, k, size);
				ret = ret + d1 * d2;
			}
			RowMajorElem (R, row, col, size) = temp * ((RowMajorElem (A, row, col, size)) - ret);
		}
		/** compute R (row, row) **/
		ret = 0.0;
		for (k = 0; k < row; k++) {
			DataType d1 = RowMajorElem (R, row, k, size);
			ret += d1 * d1;
		}
		ret = RowMajorElem (A, row, row, size) - ret;
		ret = static_cast <DataType> (sqrt (ret));
		RowMajorElem (R, row, row, size) = ret;
	} /** for-row **/


#ifdef CHECK_ASSUMPTIONS
	/** R has to be a lower triangular matrix **/
	if (! Cholesky::_is_lower_triangular <DataType, SizeType> (R, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Output Error: output matrix (ptr = " << R
			<< ") should be a lower triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif
	
	return;

} /** Cholesky :: potrf **/

/** --------------------------------------------- **/

/** TRSM : R = A / trans(B) **/
template <typename DataType, typename SizeType>
static void Cholesky :: trsm (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	SizeType matrix_elems = size * size;

	/** B is a lower triangular matrix **/
	if (Cholesky::_is_lower_triangular <DataType, SizeType> (B, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") input matrix B is a lower triangular matrix"<< endl;
		cout << ss.str (); cout.flush (); ss.str ("");

		DataType * transB = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
		DataType * invB = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

		Cholesky :: transpose <DataType, SizeType> (B, transB, size, false); /** transB = trans (B) **/
		Cholesky :: trtri <DataType, SizeType> (transB, invB, size, false, true); /** invB = inv (transB) **/

		Cholesky :: multi2 <DataType, SizeType> (A, invB, R, size, true); /** R = A * invB **/

		Cholesky::sfree <DataType, SizeType> (transB, matrix_elems);
		Cholesky::sfree <DataType, SizeType> (invB, matrix_elems);

		return;
	}
	else {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") input matrix B is NOT a lower triangular matrix"<< endl;
		cout << ss.str (); cout.flush (); ss.str ("");
	}
		
	/** B is a upper triangular matrix **/
	if (Cholesky::_is_upper_triangular <DataType, SizeType> (B, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") input matrix B is an upper triangular matrix"<< endl;
		cout << ss.str (); cout.flush (); ss.str ("");

		DataType * transB = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
		DataType * invB = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

		Cholesky :: transpose <DataType, SizeType> (B, transB, size, false); /** transB = trans (B) **/
		Cholesky :: trtri <DataType, SizeType> (transB, invB, size, true, true); /** invB = inv (transB) **/

		Cholesky :: multi2 <DataType, SizeType> (A, invB, R, size, true); /** R = A * invB **/

		Cholesky::sfree <DataType, SizeType> (transB, matrix_elems);
		Cholesky::sfree <DataType, SizeType> (invB, matrix_elems);

		return;
	}
	else {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") input matrix B is NOT an upper triangular matrix"<< endl;
		cout << ss.str (); cout.flush (); ss.str ("");
	}

	/** B is a symmetric matrix (Cholesky method) **/
	if (Cholesky::_is_symmetric <DataType, SizeType> (B, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") input matrix B is a symmetric matrix"<< endl;
		cout << ss.str (); cout.flush (); ss.str ("");

		DataType * L = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
		DataType * invL = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
		
		Cholesky :: potrf <DataType, SizeType> (B, L, size, true); /** B = L * trans (L) **/
		Cholesky :: trtri <DataType, SizeType> (L, invL, size, true, true); /** invL = inv(L) **/
		Cholesky :: trmm4 <DataType, SizeType> (A, invL, L, size, true); /** L = A * trans (invL) **/
		Cholesky :: multi2 <DataType, SizeType> (L, invL, R, size, true); /** R = L * invL **/
		
		Cholesky::sfree <DataType, SizeType> (invL, matrix_elems);	
		Cholesky::sfree <DataType, SizeType> (L, matrix_elems);

		return;
	}
	else {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") input matrix B is NOT a symmetric matrix"<< endl;
		cout << ss.str (); cout.flush (); ss.str ("");
	}

	/** general dense matrix **/

	DataType * transB = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * L = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * U = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * invL = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * invU = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * temp = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

	Cholesky :: transpose <DataType, SizeType> (B, transB, size, false); /** transB = trans (B) **/
	Cholesky :: lu <DataType, SizeType> (transB, L, U, size);  /** transB = L * U **/

	Cholesky :: trtri <DataType, SizeType> (L, invL, size, true, true); /** invL = inv (L) **/
	Cholesky :: trtri <DataType, SizeType> (U, invU, size, false, true); /** invU = inv (U) **/

	Cholesky :: multi2 <DataType, SizeType> (A, invU, temp, size, true); /** temp = A * inv(U) **/
	Cholesky :: multi2 <DataType, SizeType> (temp, invL, R, size, true); /** R = A * inv(U) * inv(L) **/

	Cholesky::sfree <DataType, SizeType> (transB, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (L, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (U, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (invL, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (invU, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (temp, matrix_elems);

	return;
}	/** Cholesky :: trsm **/

/** --------------------------------------------- **/

/** TRSM2 : solve the matrix equation : R * B = A or R = A / B **/
template <typename DataType, typename SizeType>
static void Cholesky :: trsm2 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

#ifdef CHECK_ASSUMPTIONS
	/** B has to be a symmetric matrix **/
	if (! Cholesky::_is_symmetric <DataType, SizeType> (B, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Input Error: input matrix (ptr = " << B
			<< ") should be a symmetric matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
	if (! Cholesky::_is_invertible <DataType, SizeType> (B, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId() << ") Input Error: input matrix (ptr = " << B
			<< ") should be invertible so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

	SizeType matrix_elems = size * size;

#if 1  /** B is a SPD matrix (Cholesky method) **/

	DataType * L = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * invL = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
		
	Cholesky :: potrf <DataType, SizeType> (B, L, size, true); /** B = L * trans (L) **/
	Cholesky :: trtri <DataType, SizeType> (L, invL, size, true, true); /** invL = inv(L) **/
	Cholesky :: trmm4 <DataType, SizeType> (A, invL, L, size, true); /** L = A * trans (invL) **/
	Cholesky :: multi2 <DataType, SizeType> (L, invL, R, size, true); /** R = L * invL **/
		
	Cholesky::sfree <DataType, SizeType> (invL, matrix_elems);	
	Cholesky::sfree <DataType, SizeType> (L, matrix_elems);

#else /** B is a non-SPD matrix (LU method) **/

	DataType * L = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * U = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * invL = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * invU = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * temp = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

	Cholesky :: lu <DataType, SizeType> (B, L, U, size);  /** transB = L * U **/

	Cholesky :: trtri <DataType, SizeType> (L, invL, size, true, true); /** invL = inv (L) **/
	Cholesky :: trtri <DataType, SizeType> (U, invU, size, false, true); /** invU = inv (U) **/

	Cholesky :: multi2 <DataType, SizeType> (A, invU, temp, size, true); /** temp = A * inv(U) **/
	Cholesky :: multi2 <DataType, SizeType> (temp, invL, R, size, true); /** R = A * inv(U) * inv(L) **/

	Cholesky::sfree <DataType, SizeType> (L, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (U, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (invL, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (invU, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (temp, matrix_elems);
#endif

	return;
} /** Cholesky :: trsm2 **/

/** --------------------------------------------- **/

/** TRSM3 : solve the matrix equation : R = - A / B **/
template <typename DataType, typename SizeType>
static void Cholesky :: trsm3 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	Cholesky::trsm2 <DataType, SizeType> (A, B, R, size, refresh);
	
	SizeType matrix_elem = size * size;

	for (SizeType elem = 0; elem < matrix_elem; elem ++)  R [elem] = - R[elem];

} /** Cholesky :: trsm3 **/

/** --------------------------------------------- **/

/** SYRK : R = A - B * trans (B) **/
template <typename DataType, typename SizeType>
static void Cholesky :: syrk (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0;
			for (SizeType k = 0; k < size; k ++)
				ret += RowMajorElem (B, row, k, size) * RowMajorElem (B, col, k, size);
			DataType da = RowMajorElem (A, row, col, size) - ret;
			if (refresh)  	RowMajorElem (R, row, col, size) = da;
			else  RowMajorElem (R, row, col, size) += da;
		}
	}

	return;
} /** Cholesky :: syrk **/

/** --------------------------------------------- **/

/** SYRK2 : R = A + trans (B) * B **/
template <typename DataType, typename SizeType>
static void Cholesky :: syrk2 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0;
			for (SizeType k = 0; k < size; k ++) /** trans (B) * B **/
				ret += RowMajorElem (B, k, row, size) * RowMajorElem (B, k, col, size);
			DataType da = RowMajorElem (A, row, col, size) + ret; /** A + trans (B) * B **/
			if (refresh)  	RowMajorElem (R, row, col, size) = da;
			else  RowMajorElem (R, row, col, size) += da;
		}
	}

	return;
} /** Cholesky :: syrk2 **/

/** --------------------------------------------- **/

/** GEMM : R = A - B * trans (C) **/
template <typename DataType, typename SizeType>
static void Cholesky :: gemm (DataType * A, DataType * B, DataType * C, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = B * trans (C) **/
				ret += RowMajorElem (B, row, k, size) * RowMajorElem (C, col, k, size);
			DataType da = RowMajorElem (A, row, col, size) - ret; /** A - B * trans (C) **/
			if (refresh) RowMajorElem (R, row, col, size) = da;
			else RowMajorElem (R, row, col, size) += da;
		}
	}

} /** Cholesky :: gemm **/

/** --------------------------------------------- **/

/** GEMM2 : R = A - B * C **/
template <typename DataType, typename SizeType>
static void Cholesky :: gemm2 (DataType * A, DataType * B, DataType * C, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = B * C **/
				ret += RowMajorElem (B, row, k, size) * RowMajorElem (C, k, col, size);
			DataType da = RowMajorElem (A, row, col, size) - ret; /** A - B * C **/
			if (refresh)  	RowMajorElem (R, row, col, size) = da;
			else RowMajorElem (R, row, col, size) += da;
		}
	}

	return;
} /** Cholesky :: gemm2 **/

/** --------------------------------------------- **/

/** GEMM3 : R = A + trans (B) * C **/
template <typename DataType, typename SizeType>
static void Cholesky :: gemm3 (DataType * A, DataType * B, DataType * C, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = trans (B) * C **/
				ret += RowMajorElem (B, k, row, size) * RowMajorElem (C, k, col, size);
			DataType da = RowMajorElem (A, row, col, size) + ret; /** A + trans (B) * C **/
			if (refresh)  	RowMajorElem (R, row, col, size) = da;
			else RowMajorElem (R, row, col, size) += da;
		}
	}

	return;
} /** Cholesky :: gemm3 **/

/** --------------------------------------------- **/

/** TRTRI : compute the inverse of a upper or lower triangular matrix A
	R = inv (A) where A has to be a upper or lower triangular matrix,
	default : A is a lower triangular matrix (lo_triangular = true) **/
template <typename DataType, typename SizeType>
static void Cholesky :: trtri (DataType * A, DataType * R, SizeType size, bool lo_triangular, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

#ifdef CHECK_ASSUMPTIONS
	bool is_lo_tri = Cholesky::_is_lower_triangular <DataType, SizeType> (A, size, true);
	bool is_up_tri = Cholesky::_is_upper_triangular <DataType, SizeType> (A, size, true);
	if (lo_triangular == true) {
		if (is_lo_tri == false) {
			if (is_up_tri == true) {
				cout << "Input Warn: (ptr = " << A << ") should be a lower triangular matrix but it is a upper triangular matrix." << endl;
				cout.flush ();
				lo_triangular = false;
			} else {
				cout << "Input Error: (ptr = " << A << ") should be a lower triangular matrix but it is a dense matrix so "
					<< __FUNCTION__ << " aborted" << endl;  cout.flush ();
				press_key ();
				return;
			}
		}
	} else {
		if (is_up_tri == false) {
			if (is_lo_tri == true) {
				cout << "Input Warn: (ptr = " << A << ") should be a upper triangular matrix but it is a lower triangular matrix." << endl;
				cout.flush ();
				lo_triangular = true;
			} else {
				cout << "Input Error: (ptr = " << A << ") should be a upper triangular matrix but it is a dense matrix so "
					<< __FUNCTION__ << " aborted" << endl;  cout.flush ();
				press_key ();
				return;
			}
		}
	}
	if (is_lo_tri == false && is_up_tri == false) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId () << ") Input Error: input matrix (ptr =" << A
			<< ") should be a " << ((lo_triangular) ? "lower" : "upper")
			<< " triangular matrix so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
	if (! Cholesky::_is_invertible <DataType, SizeType> (A, size)) {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId () << ") Input Error: input matrix (ptr =" << A
			<< ") should be invertible so " << __FUNCTION__ << " aborted" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
		press_key ();
		return;
	}
#endif

	if (lo_triangular)
		Cholesky::inverse_lo_triangular <DataType, SizeType> (A, R, size, refresh);
	else
		Cholesky::inverse_up_triangular <DataType, SizeType> (A, R, size, refresh);

	return;
} /** Cholesky :: trtri **/

/** --------------------------------------------- **/

/** LAUUM : R = trans (A) * A **/
template <typename DataType, typename SizeType>
static void Cholesky :: lauum (DataType * A, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = trans (A) * A **/
				ret += RowMajorElem (A, k, row, size) * RowMajorElem (A, k, col, size);
			if (refresh)  	RowMajorElem (R, row, col, size) = ret;
			else RowMajorElem (R, row, col, size) += ret;
		}
	}

	return;
} /** Cholesky :: lauum **/

/** --------------------------------------------- **/

/** TRMM : R = A * B **/
template <typename DataType, typename SizeType>
static void Cholesky :: trmm (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
#if 0
	Cholesky::multi2 <DataType, SizeType> (A, B, R, size, refresh); /** R = A * B **/
#endif

	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = A * B **/
				ret += RowMajorElem (A, row, k, size) * RowMajorElem (B, k, col, size);
			if (refresh)  	RowMajorElem (R, row, col, size) = ret;
			else RowMajorElem (R, row, col, size) += ret;
		}
	}

	return;
} /** Cholesky :: trmm **/

/** --------------------------------------------- **/

/** TRMM2 : R = - A * B **/
template <typename DataType, typename SizeType>
static void Cholesky :: trmm2 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = A * B **/
				ret += RowMajorElem (A, row, k, size) * RowMajorElem (B, k, col, size);
			ret = (-1) * ret; /** minus sign **/
			if (refresh)  	RowMajorElem (R, row, col, size) = ret;
			else RowMajorElem (R, row, col, size) += ret;
		}
	}

	return;
} /** Cholesky :: trmm2 **/

/** --------------------------------------------- **/

/** TRMM3 : R = trans (A) * B **/
template <typename DataType, typename SizeType>
static void Cholesky :: trmm3 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = trans (A) * B **/
				ret += RowMajorElem (A, k, row, size) * RowMajorElem (B, k, col, size);
			if (refresh)  	RowMajorElem (R, row, col, size) = ret;
			else RowMajorElem (R, row, col, size) += ret;
		}
	}

	return;
} /** Cholesky :: trmm3 **/

/** --------------------------------------------- **/

/** TRMM4 : R = A * trans (B) for Cholesky :: trsm (only) **/
template <typename DataType, typename SizeType>
static void Cholesky :: trmm4 (DataType * A, DataType * B, DataType * R, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);

	for (SizeType row = 0; row < size; row ++) {
		for (SizeType col = 0; col < size; col ++) {
			DataType ret = 0.0;
			for (SizeType k = 0; k < size; k ++) /** compute ret = trans (A) * B **/
				ret += RowMajorElem (A, row, k, size) * RowMajorElem (B, col, k, size);
			if (refresh)  	RowMajorElem (R, row, col, size) = ret;
			else RowMajorElem (R, row, col, size) += ret;
		}
	}

	return;
} /** Cholesky :: trmm4 **/


#endif /** _CHOLESKY_H **/