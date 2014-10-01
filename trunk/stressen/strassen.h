#ifndef _STRASSEN_H
#define _STRASSEN_H
/**
Data Mapping Algrotihm based Strassen Algorithm
**/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "boost\atomic.hpp"
#include "procinfo.h"
#include "timer.h"
#include "print_bytes.h"
#include "data_module.h"

using namespace std;

#define RowMajorElem(Matrix,Row,Col,Size) ((Matrix)[(Row)*(Size)+(Col)])
#define ColMajorElem(Matrix,Row,Col,Size) ((Matrix)[(Row)+(Col)*(Size)])

extern boost::mutex counter_mutex;

class Strassen
{
public:
	static size_t _maxMemoryUse;
	static size_t _nowMemoryUse;

public:
	Strassen () {}
	virtual ~Strassen () {}

	/** Here, we always assume that A is a row-major matrix and B is a column-major matrix,
	the resulting C matrix is a row-major matrix. **/

	/** partition a matrix into a set of submatrices **/
	template <typename DataType, typename SizeType>
	static void partition (DataType * A, bool row_major, SizeType size, int partition_size, int start_id,
		DataModuleContainer <DataType> & buffer);

	/** copy back submatrices to a matrix **/
	template <typename DataType, typename SizeType>
	static void copyback (DataType * A, bool row_major, SizeType size, int partition_size, int start_id,
		DataModuleContainer <DataType> & buffer);

	/** C = A * B (refresh, default) or C+ = A * B (not refresh) **/
	template <typename DataType, typename SizeType>
	static void multi (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh = true);

	/** C = A + B (refresh, default) or C+ = A + B (not refresh)
	A/B/C are the row-major (or column-major) matrices.
	**/
	template <typename DataType, typename SizeType>
	static void add (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh = true);

	/** C = A + B (refresh, default) or C+ = A + B (not refresh)
	A/C are row-major matrices and B is a column-major matrix.
	**/
	template <typename DataType, typename SizeType>
	static void add2 (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh = true);

	/** C = A - B (refresh, default) or C+ = A - B (not refresh)
	A/B/C are the row-major (or column-major) matrices.
	**/
	template <typename DataType, typename SizeType>
	static void minus (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh = true);

	/** C = A - B (refresh, default) or C+ = A - B (not refresh)
	A/C are row-major matrices and B is a column-major matrix.
	**/
	template <typename DataType, typename SizeType>
	static void minus2 (DataType *A, DataType *B, DataType *C, SizeType size, bool refresh = true);

	/** P = (A1 + A2) * (B1 + B2) if plusA = true and plusB = true (refresh, default)
	P += (A1 + A2) * (B1 + B2) if plusA = true and plusB = true (not refresh)  **/
	template <typename DataType, typename SizeType>
	static void pmulti (DataType *A1, DataType *A2, DataType *B1, DataType *B2, DataType *P,
		SizeType size, bool plusA, bool plusB, bool refresh = true);

	/** P = A1 + A2 - A3 + A4 (refresh, default) or  P += A1 + A2 - A3 + A4 (not refresh) **/
	template <typename DataType, typename SizeType>
	static void addMinusAdd (DataType *A1, DataType *A2, DataType *A3, DataType *A4, DataType *P,
		SizeType size, bool refresh = true);

	template <typename DataType, typename SizeType>
	static DataType * scalloc (SizeType size) {
		DataType * Ptr = (DataType *) calloc (size, sizeof (DataType));
		if (Ptr == NULL) { 
			cout << "(" << ProcInfo::GetThreadId () << ") cannot allocate memory for "
				<< (size * sizeof (DataType)) << " bytes)" << endl;
			cout.flush ();
			return NULL;
		} else {
#if defined DEBUG && DEBUG > 2
			cout << "(" << ProcInfo::GetThreadId () << ") allocate memory (ptr = "
				<< Ptr << ", size = " << PrintBytes (size * sizeof (DataType)) << ")" << endl;
			cout.flush ();
#endif
		}
		_nowMemoryUse = _nowMemoryUse + size * sizeof (DataType);
		_maxMemoryUse = max (_maxMemoryUse, _nowMemoryUse);
		writeLog (Strassen::_nowMemoryUse, memlog);
		return Ptr;
	} /** scalloc **/

	template <typename DataType, typename SizeType>
	static void sfree (DataType *&Ptr, SizeType total_elems) {
#if defined DEBUG && DEBUG > 2
		cout << "(" << ProcInfo::GetThreadId () << ") deallocate memory (ptr = "
			<< Ptr << ", size = " << PrintBytes (total_elems * sizeof (DataType)) << ")" << endl;
		cout.flush ();
#endif
		if (Ptr != NULL) free (Ptr); 
		Ptr = NULL;
		_nowMemoryUse = max (_nowMemoryUse - total_elems * sizeof (DataType), 0);
		writeLog (Strassen::_nowMemoryUse, memlog);
	} /** sfree **/

	template <typename DataType, typename SizeType>
	static void zeros (DataType * Ptr, SizeType size) {
		if (Ptr == NULL) return;
		memset (Ptr, 0, sizeof (DataType) * size * size);
	} /** zeros **/

	template <typename SizeType>
	static void chgMemoryUse (SizeType size, bool addFlag) {
		boost::mutex::scoped_lock lock (counter_mutex);
		if (addFlag) {
			Strassen::_nowMemoryUse += size;
			Strassen::_maxMemoryUse = max (Strassen::_maxMemoryUse, Strassen::_nowMemoryUse);
		} else {
			if (Strassen::_nowMemoryUse <= size) Strassen::_nowMemoryUse = 0;
			else Strassen::_nowMemoryUse -= size;
		}
	} /** chgMemoryUse **/

	static void resetMemoryCounter () {
		cout << "Reset Strassen::_nowMemoryUse "
			<< PrintBytes (Strassen::_nowMemoryUse) << " to 0" << endl;
		cout << "Reset Strassen::_maxMemoryUse "
			<< PrintBytes (Strassen::_maxMemoryUse) << " to 0" << endl;
		Strassen::_nowMemoryUse = 0;
		Strassen::_maxMemoryUse = 0;
	}

}; /** class Strassen **/


/** implementation **/
template <typename DataType, typename SizeType>
static void Strassen :: partition (DataType * A, bool row_major,
	SizeType size, int partition_size, int start_id,
	DataModuleContainer <DataType> & buffer)
{
	SizeType block_size = static_cast <SizeType> (size / partition_size);

	if ( (size % partition_size ) != 0) {
		cout << "Error : " << size << " must be a multiple of " << partition_size << endl;
		cout.flush ();
		return;
	} else {
		cout << "(" << ProcInfo::GetThreadId () << ") data partition "
			<< size << " x " << size << " matrix (ptr = " << A << ") "
			<< ((row_major) ? "row_major" : "column_major")
			<< " to " << partition_size << " x " << partition_size << " submatrices (each = "
			<< block_size << " x " << block_size << ") start_id = " << start_id
			<< endl; cout.flush ();
	}

	DataModule <DataType> * block = NULL;
	int id = start_id;

	for (int i = 0; i < partition_size; i++) {
		for (int j = 0; j < partition_size; j++) {
			block = new DataModule <DataType> (id, block_size, -1);
			if (block == NULL) {
				cout << "Error: memory not enough" << endl; cout.flush ();
				return;
			}
			DataType * ptr = block->data ();
			/** data copy **/
			for (int row = 0; row < block_size; row ++) {
				for (int col = 0; col < block_size ; col ++) {
					if (row_major) { /** row major matrix **/
						RowMajorElem (ptr, row, col, block_size) = RowMajorElem (A, row + i * block_size, col + j * block_size, size);
					} else {
						ColMajorElem (ptr, row, col, block_size) = ColMajorElem (A, row + i * block_size, col + j * block_size, size);
					}
				} /** end of col-loop **/
			} /** end of row-loop **/
			buffer.add (block); /** push a new data block to the buffer **/
			id += 1; /** next block id **/
		} /** end of j-loop **/
	} /** end of i-loop **/

} /** Strassen :: partition **/


template <typename DataType, typename SizeType>
static void Strassen :: copyback (DataType * A, bool row_major,
	SizeType size, int partition_size, int start_id,
	DataModuleContainer <DataType> & buffer)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType block_size = static_cast <SizeType> (size / partition_size);

	if ( (size % partition_size ) != 0) {
		cout << "Error : " << size << " must be a multiple of " << partition_size << endl;
		cout.flush ();
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
	int id = start_id;

	for (int i = 0; i < partition_size; i++) {
		for (int j = 0; j < partition_size; j++) {
			id = start_id + i * partition_size + j; /** id of the next data module **/
			block = buffer [id];
			if (block == NULL) {
				cout << "Error: cannot find the data module, id = " << id << endl;
				cout.flush ();
				return;
			}
			DataType * ptr = block->data ();
			/** data copy **/
			for (int row = 0; row < block_size; row ++) {
				for (int col = 0; col < block_size ; col ++) {
					if (row_major) { /** row major matrix **/
						RowMajorElem (A, row + i * block_size, col + j * block_size, size) = RowMajorElem (ptr, row, col, block_size);
					} else {
						ColMajorElem (A, row + i * block_size, col + j * block_size, size) = ColMajorElem (ptr, row, col, block_size);
					}
				} /** end of col-loop **/
			} /** end of row-loop **/
		} /** end of j-loop **/
	} /** end of i-loop **/

} /** Strassen :: copyback **/


/** A/C is a row-major matrix and B is a column-major matrix.
C = A * B. A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen::multi (DataType * A, DataType * B, DataType * C, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType r, c, elem;
	DataType * ptr1, * ptr2, ret;

	if (refresh) {
		for (r = 0; r < size; r++) {
			ptr1 = A + r * size;
			for (c = 0; c < size; c++) {
				ptr2 = B + c * size;
				for (ret = 0, elem = 0; elem < size; elem++)
					ret += ptr1 [elem] * ptr2 [elem];
				C [r * size + c] = ret; /** refresh **/
			}
		}
	}
	else {
		for (r = 0; r < size; r++) {
			ptr1 = A + r * size;
			for (c = 0; c < size; c++) {
				ptr2 = B + c * size;
				for (ret = 0, elem = 0; elem < size; elem++)
					ret += ptr1 [elem] * ptr2 [elem];
				C [r * size + c] += ret; /** not refresh **/
			}
		}
	}
} /** Strassen::multi **/

/** A/C is a row-major matrix and B is a row-major matrix.
C = A + B. (refresh = true) or C += A + B. (refresh = false)
A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen::add (DataType * A, DataType * B, DataType * C, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType i, r, c;
	DataType * ptr1 = C, * ptr2 = A, * ptr3 = B;

	if (refresh) {
		for (i = 0; i < size * size; i++)
			ptr1 [i] = ptr2 [i] + ptr3 [i];  /** refresh **/
	}
	else {
		for (i = 0; i < size * size; i++)
			ptr1 [i] += ptr2 [i] + ptr3 [i];  /** not refresh **/
	}
} /** Strassen::add **/


/** A/C is a row-major matrix and B is a column-major matrix.
C = A + B. (refresh = true) or C += A + B. (refresh = false)
A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen::add2 (DataType * A, DataType * B, DataType * C, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType r, c, c1;
	DataType * ptr1, * ptr2, * ptr3;

	if (refresh) {
		for (r = 0; r < size; r++) {
			ptr1 = C + r * size;
			ptr2 = A + r * size;
			ptr3 = B + r;
			for (c = 0, c1 = 0; c < size; c++) {
				ptr1 [c] = ptr2 [c] + ptr3 [c1]; /** refresh **/
				c1+= size;
			}
		}
	}
	else {
		for (r = 0; r < size; r++) {
			ptr1 = C + r * size;
			ptr2 = A + r * size;
			ptr3 = B + r;
			for (c = 0, c1 = 0; c < size; c++) {
				ptr1 [c] += ptr2 [c] + ptr3 [c1]; /** not refresh **/
				c1+= size;
			}
		}
	}
} /** Strassen::add2 **/


/** A/C is a row-major matrix and B is a row-major matrix.
C = A - B. (refresh = true) or C += (A - B). (refresh = false)
A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen::minus (DataType * A, DataType * B, DataType * C, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType i = 0, total_size = size * size;

	if (refresh) {
		for (i = 0; i < total_size; i ++)
			C [i] = A [i] - B [i];  /** refresh **/
	}
	else {
		for (i = 0; i < total_size; i ++)
			C [i] += (A [i] - B [i]);  /** not refresh **/
	}
} /** Strassen::minus **/


/** A/C is a row-major matrix and B is a column-major matrix.
C = A - B. (refresh = true) or C += A - B. (refresh = false)
A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen::minus2 (DataType * A, DataType * B, DataType * C, SizeType size, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType r, c, c1;
	DataType * ptr1, * ptr2, * ptr3;

	if (refresh) {
		for (r = 0; r < size; r++) {
			ptr1 = C + r * size;
			ptr2 = A + r * size;
			ptr3 = B + r;
			for (c = 0, c1 = 0; c < size; c++) {
				ptr1 [c] = ptr2 [c] - ptr3 [c1]; /** refresh **/
				c1+= size;
			}
		}
	}
	else {
		for (r = 0; r < size; r++) {
			ptr1 = C + r * size;
			ptr2 = A + r * size;
			ptr3 = B + r;
			for (c = 0, c1 = 0; c < size; c++) {
				ptr1 [c] += (ptr2 [c] - ptr3 [c1]); /** not refresh **/
				c1+= size;
			}
		}
	}
} /** Strassen::minus2 **/

/** A/C is a row-major matrix and B is a column-major matrix.
C = (A1 + A2) * (B1 + B2). A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen :: pmulti (DataType *A1, DataType *A2, DataType *B1, DataType *B2, DataType *P,
		SizeType size, bool plusA, bool plusB, bool refresh)
{
	SimpleTimer timer (__FUNCTION__);
	SizeType matrix_elems = size * size;

	bool newP1 = false, newP2 = false;
	DataType * P1 = NULL, * P2 = NULL;

	if (A1 == NULL && A2 != NULL) P1 = A2;
	if (A1 != NULL && A2 == NULL) P1 = A1;
	if (A1 != NULL && A2 != NULL) {
		P1 = Strassen::scalloc <DataType, SizeType> (size * size);
		newP1 = true;
		if (plusA)  Strassen::add (A1, A2, P1, size);
		else Strassen::minus (A1, A2, P1, size);
	}

	if (B1 == NULL && B2 != NULL) P2 = B2;
	if (B1 != NULL && B2 == NULL) P2 = B1;
	if (B1 != NULL && B2 != NULL) {
		P2 = Strassen::scalloc <DataType, SizeType> (size * size);
		newP2 = true;
		if (plusB) Strassen::add (B1, B2, P2, size);
		else Strassen::minus (B1, B2, P2, size);
	}

	Strassen::multi (P1, P2, P, size, refresh);

	if (newP1) Strassen::sfree <DataType, SizeType> (P1, size * size);
	if (newP2) Strassen::sfree <DataType, SizeType> (P2, size * size);
} /** Strassen :: pmulti **/

/** A/C is a row-major matrix and B is a column-major matrix.
C = A1 + A2 - A3 + A4. A/B/C is a square matrix of size x size **/
template <typename DataType, typename SizeType>
void Strassen :: addMinusAdd (DataType *A1, DataType *A2, DataType *A3, DataType *A4, DataType *P,
		SizeType size, bool refresh)
{
	if (A1 == NULL || A2 == NULL) {
		cout << "Runtime error: A1 and A2 must be NOT a NULL in " 
			<< __LINE__ << " of " << __FUNCTION__ << endl;
		cout.flush ();
		return;
	}
	Strassen::add <DataType, SizeType> (A1, A2, P, size, refresh);
	if (A3 != NULL)  Strassen::minus <DataType, SizeType> (P, A3, P, size, true);
	if (A4 != NULL)  Strassen::add <DataType, SizeType> (P, A4, P, size, true);
} /** Strassen :: addMinusAdd **/



#endif /** #ifndef _STRASSEN_H **/