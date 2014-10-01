/**
Data Mapping Algrotihm based Strassen Algorithm
**/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <time.h>

#include "boost\lockfree\queue.hpp"
#include "boost\lockfree\stack.hpp"

#include "utility.h"
#include "procinfo.h"
#include "strassen.h"
#include "data_module.h"
#include "print_bytes.h"
#include "thread_pool.h"
#include "tasks.h"
#include "data_mapping.h"

using namespace std;

typedef long double ElemType;
typedef int SizeType;

const double epsilon = 0.00001; /** error tolerance **/

void usage (const char *);

template <typename DataType, typename SizeType>
long double arrayRMSD (DataType *, DataType *, SizeType);

template <typename DataType, typename SizeType>
void tradMatrixMatrixMulti (DataType * A, DataType * B, DataType * C, SizeType matrix_size);

template <typename DataType, typename SizeType>
void tradMatrixMatrixAdd (DataType * A, DataType * B, DataType * C, SizeType matrix_size);

template <typename DataType, typename SizeType>
bool isEqual (DataType * A, DataType * B, SizeType total_size, DataType tolerance = epsilon);

boost::lockfree::queue <int> newAvailModules (100);

LogFile logger ("runtime.log", 0);
LogFile errlog ("runtime.err", 0);
LogFile memlog ("runtime.mem", 0);

#define MAT_MULTI
// #define DO_VERIFY

/**********************/
int main (int argc, char * argv [])
{
	if (argc < 3) {
		usage (argv [0]);
		return 1;
	}

	bool keep_running = true;
	string func_name (argv [0]);
	SizeType matrix_size = (SizeType) atoi ((argv [1])); /** matrix size **/
	string dep_filename (argv [2]); /** filename **/
	SizeType matrix_elems = matrix_size * matrix_size; /** square matrix **/
	string data_type (typeid(ElemType).name ());

	/** number of threads **/
	int threads = (argc > 3) ? (max(1, atoi(argv[3]))) : (1);

	int partition_size = 2; /** partition method (config) **/
	int new_id;

	ThreadPool <ElemType> scheduler (threads); /** scheduler pool **/
	DataModuleContainer <ElemType> buffer; /** memory buffers **/
	DataMapping depd_matrix (dep_filename.c_str (), 2, 100); /** data depd matrix **/
	DataDepElem * next_task = NULL;

	SimpleTimer * t_prep, * t_compt, * t_post; /** timers **/

	cout << "<<<< =========== Program Setup =========== >>>>" << endl; cout.flush ();
	cout << func_name << " is running on " << ProcInfo () << endl;
	cout << "Matrix size = " << matrix_size << " x " << matrix_size << endl;
	cout << "Data type = " << data_type << " (size = " << sizeof (ElemType) << ")" << endl;
	cout << "Data Volume per Matrix = " << PrintBytes (matrix_elems * sizeof (ElemType)) << endl;
	cout << "Data Dependency Matrix Filename = " << dep_filename << endl;
	cout << "Threads = " << threads << endl;
	cout.flush ();

	if (matrix_size % partition_size != 0) {
		cerr << "matrix_size (" << matrix_size << ") has to be a multiple of partition_size ("
			<< partition_size << ")" << endl;
		return 2;
	}

	cout << "<<<< =========== Data preparation =========== >>>>" << endl; cout.flush ();
	t_prep = new SimpleTimer ("Data preparation");

	ElemType * A = Strassen::scalloc <ElemType, SizeType> (matrix_elems);
	ElemType * B = Strassen::scalloc <ElemType, SizeType> (matrix_elems);
	ElemType * C = Strassen::scalloc <ElemType, SizeType> (matrix_elems);

	/** row or column major allocation (this is the assumption, no change necessary) **/
	bool A_row_major = true;
#if defined MAT_MULTI
	bool B_row_major = false;
#else
	bool B_row_major = true;
#endif
	bool C_row_major = true;

	genRandomArray <ElemType, SizeType> (A, matrix_elems, 0, 10);
	genRandomArray <ElemType, SizeType> (B, matrix_elems, 0, 10);
	genRandomArray <ElemType, SizeType> (C, matrix_elems, 0, 0);

	int A_matrix_start_id, A_matrix_end_id;
	int B_matrix_start_id, B_matrix_end_id;
	int result_matrix_start_id, result_matrix_end_id;

	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	/** Assume the task id is in a row from A, to B and then to C. **/

	A_matrix_start_id = depd_matrix.minTaskId (); /** start task id of matrix A **/
	Strassen::partition <ElemType, SizeType> (A, A_row_major, matrix_size, partition_size, A_matrix_start_id, buffer);
	A_matrix_end_id = buffer.size ();

	B_matrix_start_id = buffer.size () + 1;
	Strassen::partition <ElemType, SizeType> (B, B_row_major, matrix_size, partition_size, B_matrix_start_id, buffer);
	B_matrix_end_id = buffer.size ();

	result_matrix_start_id = buffer.size () + 1;
	Strassen::partition <ElemType, SizeType> (C, C_row_major, matrix_size, partition_size, result_matrix_start_id, buffer);
	result_matrix_end_id = buffer.size ();

	/** make resultant blocks all persistent **/
	for (int i = result_matrix_start_id; i <= result_matrix_end_id; i++)
		depd_matrix.setModuleVisits (i, -1);

#ifndef DO_VERIFY /** deallocate the raw matrices **/
	Strassen::sfree <ElemType, SizeType> (A, matrix_elems);
	Strassen::sfree <ElemType, SizeType> (B, matrix_elems);
	Strassen::sfree <ElemType, SizeType> (C, matrix_elems);
#endif

	/** make matrix A avaiable **/
	for (int id = A_matrix_start_id; id <= A_matrix_end_id; id++) {
		depd_matrix.setAvailModule (id);
		buffer [id]->setVisits (depd_matrix.getModuleVisits (id));
	}

	/** make matrix B avaiable **/
	for (int id = B_matrix_start_id; id <= B_matrix_end_id; id++) {
		depd_matrix.setAvailModule (id);
		buffer [id]->setVisits (depd_matrix.getModuleVisits (id));
	}

	/** make matrix C as a persistent matrix **/
	for (int id = result_matrix_start_id; id <= result_matrix_end_id; id++) {
		buffer [id]->setVisits (depd_matrix.getModuleVisits (id));
	}

	/** make other intermediate data modules **/
	for (int id = result_matrix_end_id + 1; id <= depd_matrix.maxTaskId (); id ++) {
		int size = matrix_size / partition_size;
		int vis = depd_matrix.getModuleVisits (id);
		DataModule <ElemType> * block = new DataModule <ElemType> (id, size, vis);
		buffer.add (block);
	}

	/** find all of the available tasks (initial) **/
	keep_running = true;
	do {
		next_task = depd_matrix.next_task ();
		if (next_task != NULL) {
			TaskMap <ElemType> :: add (scheduler, next_task->task_id, 
				buffer [next_task->input [0]], buffer [next_task->input [1]], buffer [next_task->output]);
			delete next_task;
		} else { keep_running = false; }
	} while (keep_running);

	delete t_prep; /** <================= **/

#if 0 /** for testing purpose **/
	TaskMap <ElemType> :: add (scheduler, MatrixAdd, buffer [1], buffer [5], buffer [9]);
	TaskMap <ElemType> :: add (scheduler, MatrixAdd, buffer [2], buffer [6], buffer [10]);
	TaskMap <ElemType> :: add (scheduler, MatrixAdd, buffer [3], buffer [7], buffer [11]);
	TaskMap <ElemType> :: add (scheduler, MatrixAdd, buffer [4], buffer [8], buffer [12]);
#endif

	cout << "(" << ProcInfo::GetThreadId () << ") pending tasks = " << scheduler.pending () << endl; 
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	cout << "<<<< =========== Matrix Computation (Begin) =========== >>>>" << endl; cout.flush ();
	writeLog ("started", memlog);
	writeLog (Strassen::_nowMemoryUse, memlog);

	t_compt = new SimpleTimer ("Matrix Computation");

	scheduler.run (); /** run the scheduler **/

	/** append the tasks based on the availability of data modules **/
	try {

	do {
		/** search for new available modules **/
		while (newAvailModules.pop (new_id)) {
			depd_matrix.setAvailModule (new_id);
		}

		/** search for new tasks **/
		next_task = depd_matrix.next_task ();
		if (next_task != NULL) {
			TaskMap <ElemType> :: add (scheduler, next_task->task_id, 
				buffer [next_task->input [0]], buffer [next_task->input [1]], buffer [next_task->output]);
			delete next_task;
		} /**  insert a new task based on availability of data modules **/

	} while (! depd_matrix.empty ());

	} catch (std::exception & e) {
		cerr << "Runtime error: " << e.what () << endl; cerr.flush ();
	}

	/****************************************************/

	scheduler.wait (); /** wait all jobs to all done **/

	writeLog (Strassen::_nowMemoryUse, memlog);
	writeLog ("ended", memlog);
	cout << "<<<< =========== Matrix Computation (End) =========== >>>>" << endl; cout.flush ();
	delete t_compt; /** <================= **/

	cout << "(" << ProcInfo::GetThreadId () << ") number of threads = " << threads << endl;
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	cout << "<<<< =========== Result Verification =========== >>>>" << endl; cout.flush ();
	t_post = new SimpleTimer ("Result Verification");

#ifdef DO_VERIFY
	/** copy submatrices to the resultant matrix for results verification. **/
	Strassen::copyback <ElemType, SizeType> (C, C_row_major, matrix_size, partition_size, result_matrix_start_id, buffer);

#if defined MAT_MULTI
	tradMatrixMatrixMulti<ElemType, SizeType> (A, B, C, matrix_size);
#else
	tradMatrixMatrixAdd<ElemType, SizeType> (A, B, C, matrix_size);
#endif
#endif

	delete t_post; /** <================= **/

	cout << "<<<< ========================= >>>>" << endl; cout.flush ();
	buffer.clear ();

	cout << "Peak Memory Usage = " << Strassen::_maxMemoryUse << " bytes" << endl;
	cout << "Peak Memory Usage = " << PrintBytes (Strassen::_maxMemoryUse) << endl;
	cout.flush ();

	cout << "<<<< =========== ALL DONE =========== >>>>" << endl; cout.flush ();

	if (A != NULL) Strassen::sfree <ElemType, SizeType> (A, matrix_elems);
	if (B != NULL) Strassen::sfree <ElemType, SizeType> (B, matrix_elems);
	if (C != NULL) Strassen::sfree <ElemType, SizeType> (C, matrix_elems);
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

#if defined WIN32
	cout << "Press any key to continue..." << endl;  cout.flush ();
	getchar ();
#endif

	return 0;
} /** main **/

void usage (const char * prog)
{
	cout << "Usage of " << prog << endl;
	cout << prog << " <matrix_size> <data_dep_map.filename> <data_visits.filename> [number_of_threads]" << endl;
	cout << "<matrix_size> must be a positive integer and a multiple of 7" << endl;
	cout << "<data_dep_map.filename> is the filename containing the data dependency matrix" << endl;
	cout << "<data_visits.filename> is the filename containing the visit counts of data modules" << endl;
	cout << "[number_of_threads] must be a positive integer if given" << endl;
	cout.flush ();
}

template <typename DataType, typename SizeType>
long double arrayRMSD (DataType * _array1, DataType * _array2, SizeType _size)
{
	long double ret = 0.0;
	for (SizeType elem = 0; elem < _size; elem++)
		ret += (_array1[elem] - _array2[elem]) * (_array1[elem] - _array2[elem]);
	ret = ret / _size;
	ret = static_cast<long double> (pow ((double) ret, 2.0));
	return ret;
} /** arrayRMSD **/

template <typename DataType, typename SizeType>
void tradMatrixMatrixMulti (DataType * A, DataType * B, DataType * C, SizeType matrix_size)
{
	SimpleTimer * timer = NULL;
	SizeType matrix_elems = matrix_size * matrix_size;
	DataType * ret = Strassen::scalloc <DataType, SizeType> (matrix_elems);
	if (ret == NULL) return;
	timer = new SimpleTimer ("Tranditional matrix-matrix multiplication");
#if 0
	for (int row = 0; row < matrix_size; row++) {
		for (int col = 0; col < matrix_size; col ++) {
			RowMajorElem (ret, row, col, matrix_size) = 0.0;
			for (int k = 0; k < matrix_size; k++)
				RowMajorElem (ret, row, col, matrix_size) += 
					RowMajorElem (A, row, k, matrix_size) * ColMajorElem (B, k, col, matrix_size);
		}
	}
#else
	Strassen::multi <DataType, SizeType> (A, B, ret, matrix_size);
#endif
	delete timer;
	double error = static_cast<double> (arrayRMSD<DataType, SizeType> (C, ret, matrix_elems));
	std::cout << "(" << ProcInfo::GetThreadId () 
		<< ") tolerance = " << std::fixed << error 
		<< " " << ((abs(error) <= epsilon) ? "Accepted" : "Failed");
	int i = 3;
	while ((i--) > 0) { Sleep (1000); std::cout << "."; std::cout.flush (); }
	std::cout << endl;
	Strassen::sfree <DataType, SizeType> (ret, matrix_elems);
} /** tradMatrixMatrixMulti **/

template <typename DataType, typename SizeType>
void tradMatrixMatrixAdd (DataType * A, DataType * B, DataType * C, SizeType matrix_size)
{
	SimpleTimer * timer = NULL;
	SizeType matrix_elems = matrix_size * matrix_size;
	DataType * ret = Strassen::scalloc <DataType, SizeType> (matrix_elems);
	if (ret == NULL) return;
	timer = new SimpleTimer ("Tranditional matrix-matrix addition");
#if 0
	for (SizeType row = 0; row < matrix_size; row++) {
		for (SizeType col = 0; col < matrix_size; col ++) {
			RowMajorElem (ret, row, col, matrix_size) = 
				RowMajorElem (A, row, col, matrix_size) + ColMajorElem (B, row, col, matrix_size);
		}
	}
#else
	Strassen::add <DataType, SizeType> (A, B, ret, matrix_size);
#endif
	delete timer;
	double error = static_cast<double> (arrayRMSD<DataType, SizeType> (C, ret, matrix_elems));
	cout << "(" << ProcInfo::GetThreadId () << ") tolerance = "
		<< std::fixed << error << " " << ((abs(error) <= epsilon) ? "Accepted" : "Failed");
	int i = 3;
	while ((i--) > 0) { Sleep (1000); std::cout << "."; std::cout.flush (); }
	std::cout << endl;
	cout.flush ();
	Strassen::sfree <DataType, SizeType> (ret, matrix_elems);
} /** tradMatrixMatrixAdd **/

template <typename DataType, typename SizeType>
bool isEqual (DataType * A, DataType * B, SizeType total_size, DataType tolerance)
{
	for (SizeType elem = 0; elem < total_size; elem++)
		if (abs (A [elem] - B [elem]) > tolerance) return false;
	return true;
} /** isEqual **/