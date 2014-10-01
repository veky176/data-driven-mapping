/**
Data Mapping Algrotihm based Cholesky Algorithm
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
#include "cholesky.h"
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
void genRandomSPDMatrix (DataType *, SizeType, DataType, DataType);

template <typename DataType, typename SizeType>
void checkCholeskyResult (DataType * S, DataType * R, SizeType matrix_size);

template <typename DataType, typename SizeType>
void checkCholeskyResult2 (DataType * S, DataType * L, SizeType matrix_size);

boost::lockfree::queue <int> newAvailModules (100);

LogFile logger ("runtime.log", 0);
LogFile errlog ("runtime.err", 0);
LogFile memlog ("runtime.mem", 0);

#define DO_TEST 0
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

	int partition_size = 4; /** partition method (config) **/
	int new_id;
	
	ThreadPool <ElemType> scheduler (threads); /** scheduler pool **/
	DataModuleContainer <ElemType> buffer; /** memory buffers **/
	DataMapping depd_matrix (dep_filename.c_str (), 3, 100); /** data depd matrix **/
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
		cout << "matrix_size (" << matrix_size << ") has to be a multiple of partition_size ("
			<< partition_size << ")" << endl;
		return 2;
	}

	press_key ();

#if DO_TEST  /** Tests for individual Functions **/
	cout << "<<<< =========== Function tests =========== >>>>" << endl; cout.flush ();
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();
	test_lu (1000);
	test_trtri (true, 500);
	test_trtri (false, 500);
	test_trsm (1000);
	test_inverse_up_triangular ();
	test_inverse_lo_triangular ();
	test_potrf (1000);
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();
	press_key ();
#endif

	cout << "<<<< =========== Data preparation =========== >>>>" << endl; cout.flush ();
	t_prep = new SimpleTimer ("Data preparation");

	ElemType * S = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);

	/** row-major allocation (this is the assumption, no change necessary) **/

	genRandomSPDMatrix <ElemType, SizeType> (S, matrix_size, -2, 2);

	int S_matrix_start_id, S_matrix_end_id;
	int result_matrix_start_id, result_matrix_end_id;

	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	S_matrix_start_id = depd_matrix.minTaskId (); /** 1 **/
	Cholesky::partition <ElemType, SizeType> (S, matrix_size, partition_size, S_matrix_start_id, buffer); /** create SPD matrix S **/
	S_matrix_end_id = buffer.size (); /** 10 **/

#if 0 /** check up the copyback function **/
	ElemType * S2 = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);
	Cholesky::copyback <ElemType, SizeType> (S2, matrix_size, partition_size, S_matrix_start_id, buffer);
	if (is2ArrayEqual <ElemType, SizeType> (S, S2, matrix_elems, (ElemType) epsilon)) {
		cout << "(" << ProcInfo::GetThreadId() << ") Cholesky::copyback passed" << endl;
	} else {
		cout << "(" << ProcInfo::GetThreadId() << ") Cholesky::copyback failed" << endl;
	}
	cout.flush ();

	Cholesky::sfree <ElemType, SizeType> (S2, matrix_elems);
	press_key ();
#endif

	result_matrix_end_id = depd_matrix.maxTaskId (); /** should be 70 **/
	result_matrix_start_id = result_matrix_end_id + 1 - ((partition_size + 1) * (partition_size)) / 2; /** should be 61 **/
	
	cout << "S matrix data id is in [" << S_matrix_start_id << ", " << S_matrix_end_id << "]" << endl;
	cout << "Result matrix data id is in [" << result_matrix_start_id << ", " << result_matrix_end_id << "]" << endl;
	cout.flush ();

	/** make matrix S avaiable **/
	for (int id = S_matrix_start_id; id <= S_matrix_end_id; id ++) {
		depd_matrix.setAvailModule (id);
		buffer [id]->setVisits (depd_matrix.getModuleVisits (id));
	}

	/** make other intermediate data modules **/
	for (int id = S_matrix_end_id + 1; id < result_matrix_start_id; id++) {
		int size = matrix_size / partition_size;
		int vis = depd_matrix.getModuleVisits (id);
		DataModule <ElemType> * block = new DataModule <ElemType> (id, size, vis);
		buffer.add (block);
	}

	/** make the resultant matrix as a persistent matrix **/
	for (int id = result_matrix_start_id; id <= result_matrix_end_id; id++) {
		int size = matrix_size / partition_size;
		int vis = -1;
		DataModule <ElemType> * block = new DataModule <ElemType> (id, size, vis);
		buffer.add (block);
	}

#ifdef DO_VERIFY
	for (int id = 1; id <= 10; id ++) buffer [id]->setVisits (-1);
	for (int id = 21; id <= 30; id ++)  buffer [id]->setVisits (-1);
	for (int id = 41; id <= 50; id ++)  buffer [id]->setVisits (-1);
#else
	/** deallocate the raw matrices **/
	Cholesky::sfree <ElemType, SizeType> (S, matrix_elems);
#endif

	/** find all of the available tasks (initial) **/
	keep_running = true;
	do {
		next_task = depd_matrix.next_task ();
		if (next_task != NULL) {
			TaskMap <ElemType> :: add (scheduler, next_task->task_id, 
				buffer [next_task->input [0]], buffer [next_task->input [1]], buffer [next_task->input [2]],
				buffer [next_task->output]);
			delete next_task;
		} else { keep_running = false; }
	} while (keep_running);

	delete t_prep; /** <================= **/

	cout << "(" << ProcInfo::GetThreadId () << ") pending tasks = " << scheduler.pending () << endl; 
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	// press_key ();
	cout << "<<<< =========== Cholesky Factorization (Begin) =========== >>>>" << endl; cout.flush ();
	writeLog ("started", memlog);
	writeLog (Cholesky::_nowMemoryUse, memlog);

	t_compt = new SimpleTimer ("Cholesky Factorization");

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
				buffer [next_task->input [0]], buffer [next_task->input [1]], buffer [next_task->input [2]], 
				buffer [next_task->output]);
			delete next_task;
		} /**  insert a new task based on availability of data modules **/

	} while (! depd_matrix.empty ());

	} catch (std::exception & e) {
		cout << "Runtime error: " << e.what () << endl; cout.flush ();
	}

	/****************************************************/

	scheduler.wait (); /** wait all jobs to all done **/

	writeLog (Cholesky::_nowMemoryUse, memlog);
	writeLog ("ended", memlog);
	cout << "<<<< =========== Cholesky Factorization (End) =========== >>>>" << endl; cout.flush ();
	delete t_compt; /** <================= **/

	cout << "(" << ProcInfo::GetThreadId () << ") number of threads = " << threads << endl;
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	press_key ();
	cout << "<<<< =========== Result Verification =========== >>>>" << endl; cout.flush ();
	t_post = new SimpleTimer ("Result Verification");

#ifdef DO_VERIFY
	ElemType * L = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);  /** S = L * trans (L) **/
	ElemType * invL = Cholesky::scalloc <ElemType, SizeType> (matrix_elems); /** invL = inv (L)  **/
	ElemType * invS = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);  /** invS = inv (S) **/

	/** copy back to a lower triangular matrix : L such that S = L * trans (L) **/
	Cholesky::copyback2 <ElemType, SizeType> (L, matrix_size, partition_size, 21, buffer); /** 21 ~ 30 **/
	SaveMatrix <ElemType, SizeType> (L, matrix_size, "L_matrix.dat").save ();

	/** copy back to a lower triangular matrix : invL such that invL = inv (L) **/
	Cholesky::copyback2 <ElemType, SizeType> (invL, matrix_size, partition_size, 41, buffer);  /** 41 ~ 50 **/
	SaveMatrix <ElemType, SizeType> (invL, matrix_size, "invL_matrix.dat").save ();

	/** copy back to a symmetric matrix : invS such that invS =  inv (S) **/
	Cholesky::copyback <ElemType, SizeType> (invS, matrix_size, partition_size, result_matrix_start_id, buffer); /** 61 ~ 70 **/
	SaveMatrix <ElemType, SizeType> (invS, matrix_size, "invS_matrix.dat").save ();

	cout << "(" << ProcInfo::GetThreadId() << ") << STEP 1 >> validating the Cholesky factorization of a SPD matrix" << endl; cout.flush ();
	checkCholeskyResult2 <ElemType, SizeType> (S, L, matrix_size); /** check if S = L * trans (L) **/

	cout << "(" << ProcInfo::GetThreadId() << ") << STEP 2 >> validating the inversion of a lower triangular matrix" << endl; cout.flush ();
	checkCholeskyResult <ElemType, SizeType> (L, invL, matrix_size);  /** check if invL = inv (L) **/

	cout << "(" << ProcInfo::GetThreadId() << ") << STEP 3 >> validating the inversion of a SPD matrix" << endl; cout.flush ();
	checkCholeskyResult <ElemType, SizeType> (S, invS, matrix_size);  /** check if S = inv (S) **/

	Cholesky::sfree <ElemType, SizeType> (S, matrix_elems);
	Cholesky::sfree <ElemType, SizeType> (L, matrix_elems);
	Cholesky::sfree <ElemType, SizeType> (invL, matrix_elems);
	Cholesky::sfree <ElemType, SizeType> (invS, matrix_elems);

#else

	ElemType * ResultMatrix = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);
	switch (result_matrix_start_id) 
	{
		/** copy back a lower triangular matrix **/
	case 21:
	case 41:
		Cholesky::copyback2 <ElemType, SizeType> (ResultMatrix, matrix_size, partition_size, result_matrix_start_id, buffer);
		break;
		/** copy back a symmetric matrix **/
	case 61:
	default:
		Cholesky::copyback <ElemType, SizeType> (ResultMatrix, matrix_size, partition_size, result_matrix_start_id, buffer);
	}
	// SaveMatrix <ElemType, SizeType> (ResultMatrix, matrix_size, "result_matrix.dat").save ();
	Cholesky::sfree <ElemType, SizeType> (ResultMatrix, matrix_elems);

#endif

	delete t_post; /** <================= **/

	cout << "<<<< ========================= >>>>" << endl; cout.flush ();
	buffer.clear ();

	cout << "Peak Memory Usage = " << Strassen::_maxMemoryUse << " bytes" << endl;
	cout << "Peak Memory Usage = " << PrintBytes (Strassen::_maxMemoryUse) << endl;
	cout.flush ();

	cout << "<<<< =========== ALL DONE =========== >>>>" << endl; cout.flush ();

	if (S != NULL) Strassen::sfree <ElemType, SizeType> (S, matrix_elems);
	cout << "(" << ProcInfo::GetThreadId () << ") current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
	cout.flush ();

	press_key ();

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
void genRandomSPDMatrix 
	(DataType * _array, SizeType _size, DataType min, DataType max)
{
	SimpleTimer timer (__FUNCTION__);

#if 0 /** The windows random number is not that random **/
	SizeType matrix_elems = _size * _size;
	DataType random, value;
	srand (time(NULL));

	DataType * A = Cholesky::scalloc <DataType, SizeType> (matrix_elems);
	DataType * B = Cholesky::scalloc <DataType, SizeType> (matrix_elems);

	genRandomArray <DataType, SizeType> (A, matrix_elems, min, max);
	Cholesky::transpose <DataType, SizeType> (A, B, _size, false);

	Cholesky::multi2 (A, B, _array, _size, true);

	SaveMatrix <DataType, SizeType> (A, _size, "matrix_inputA.txt").save ();
	SaveMatrix <DataType, SizeType> (B, _size, "matrix_inputB.txt").save ();
	SaveMatrix <DataType, SizeType> (_array, _size, "matrix_input.txt").save ();

	Cholesky::sfree <DataType, SizeType> (A, matrix_elems);
	Cholesky::sfree <DataType, SizeType> (B, matrix_elems);
#endif

	std::stringstream filename;

	/** load a SPD matrix from file **/
	filename << "input_" << _size << ".txt";

	/** Load a lower triangular matrix from file **/
	// filename << "matrix_LoTri_" << _size << ".txt";

	std::ifstream file;
	
	file.open (filename.str().c_str());

	if (! file.is_open ()) {
		cout << "File open error: cannot open " << filename.str() << endl;
		cout.flush ();
		return;
	}

	cout << "(" << ProcInfo::GetThreadId () << ") reading " << filename.str () << endl;
	cout.flush ();
	std::string line, subline;

	for (int row = 0; row < _size; row++) 
	{
		line.clear ();
		getline (file, line);
		if (line.empty ()) {
			cout << "Error: not enough lines" << endl; cout.flush ();
			press_key ();
			return;
		}
		for (int col = 0; col < _size; col ++) {
			subline = line.substr (0, line.find_first_of (' '));
			if (subline.empty ()) {
				cout << "Error: not enough element in Line " << row << endl;
				cout.flush ();
				press_key ();
				return;
			}
			_array [row * _size + col] = static_cast <DataType> (atof (subline.c_str()));
			line.erase (0, subline.size () + 1);
		}
	}

#ifdef CHECK_ASSUMPTIONS
	if (Cholesky::_is_spd_matrix <DataType, SizeType> (_array, _size, true)) {
		cout << "(" << ProcInfo::GetThreadId () << ") input matrix is a SPD matrix" << endl;
		cout.flush ();
	}
#endif

#if 0 /** print the input matrix for double-check only **/
	cout << PrintMatrix <DataType, SizeType> (_array, _size * _size) << endl;
	press_key ();
#endif
} /** genRandomSPDMatrix **/

/** assume: S and R are both row-major and R = inv (S) **/
template <typename DataType, typename SizeType>
void checkCholeskyResult (DataType * S, DataType * R, SizeType matrix_size)
{
	SizeType matrix_elems = matrix_size * matrix_size;
	ElemType * temp = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);

	if (temp == NULL) {
		cout << "Runtime error: insufficient memory!" << endl;
		cout.flush ();
		return;
	}

	bool S_LoTri = Cholesky::_is_lower_triangular <DataType, SizeType> (S, matrix_size, true);
	bool S_UpTri = Cholesky::_is_upper_triangular <DataType, SizeType> (S, matrix_size, true);
	bool S_Sym = Cholesky::_is_symmetric <DataType, SizeType> (S, matrix_size, true);

	bool R_LoTri = Cholesky::_is_lower_triangular <DataType, SizeType> (R, matrix_size, true);
	bool R_UpTri = Cholesky::_is_upper_triangular <DataType, SizeType> (R, matrix_size, true);
	bool R_Sym = Cholesky::_is_symmetric <DataType, SizeType> (R, matrix_size, true);

	if (S_LoTri && R_LoTri == false) {
		cout << "(" << ProcInfo::GetThreadId() 
			<< ") Error: inversion of a lower triangular matrix should be a lower triangular matrix" << endl;
	}
	if (S_UpTri && R_UpTri == false) {
		cout << "(" << ProcInfo::GetThreadId() 
			<< ") Error: inversion of a upper triangular matrix should be a upper triangular matrix" << endl;
	}
	if (S_Sym && R_Sym == false) {
		cout << "(" << ProcInfo::GetThreadId() 
			<< ") Error: inversion of a symmatric matrix should be a symmatric matrix" << endl;
	}
	cout << "(" << ProcInfo::GetThreadId() << ") input matrix S is a"
		<< ((S_LoTri) ? " lower triangular" : "") << ((S_UpTri) ? " upper triangular" : "")
		<< ((S_Sym) ? " symmatric" : "") << " matrix" << endl;
	cout << "(" << ProcInfo::GetThreadId() << ") input matrix R is a"
		<< ((R_LoTri) ? " lower triangular" : "") << ((R_UpTri) ? " upper triangular" : "")
		<< ((R_Sym) ? " symmatric" : "") << " matrix" << endl;
	cout.flush ();

	Cholesky::multi2 <DataType, SizeType> (S, R, temp, matrix_size, true);
	
	long double error = 0.0;
	DataType value;

	for (SizeType row = 0; row < matrix_size; row ++) {
		for (SizeType col = 0; col < matrix_size; col ++) {
			value = RowMajorElem (temp, row, col, matrix_size);
			if (row == col)  error += (value - 1.0) * (value - 1.0);
			else  error += value * value;
		}
	}
	//error = error / matrix_elems;
	error = static_cast<long double> (pow ((double) error, 2.0));

	std::cout << "(" << ProcInfo::GetThreadId () << ") tolerance = " << std::fixed << error
		<< " " << ((abs(error) <= epsilon) ? "Accepted" : "Failed");

	int i = 3;
	while ((i--) > 0) { Sleep (1000); std::cout << "."; std::cout.flush (); }
	std::cout << endl;

	Cholesky::sfree <ElemType, SizeType> (temp, matrix_elems);
	return;
}


/** assume: S and L are both row-major and S = L * trans (L) **/
template <typename DataType, typename SizeType>
void checkCholeskyResult2 (DataType * S, DataType * L, SizeType matrix_size)
{
	SizeType matrix_elems = matrix_size * matrix_size;
	ElemType * temp = Cholesky::scalloc <ElemType, SizeType> (matrix_elems);

	if (temp == NULL) {
		cout << "Runtime error: insufficient memory!" << endl;
		cout.flush ();
		return;
	}

	bool S_Sym = Cholesky::_is_symmetric <DataType, SizeType> (S, matrix_size, true);
	bool L_LoTri = Cholesky::_is_lower_triangular <DataType, SizeType> (L, matrix_size, true);

	if (! S_Sym) {
		cout << "(" << ProcInfo::GetThreadId() << ") Error: input matrix S should be a symmatric matrix" << endl;
	}
	if (! L_LoTri) {
		cout << "(" << ProcInfo::GetThreadId() << ") Error: input matrix L should be a lower triangular matrix" << endl;
	}

	cout << "(" << ProcInfo::GetThreadId() << ") input matrix S is"
		<< ((S_Sym) ? "" : " not") << " a symmetric matrix" << endl;
	cout << "(" << ProcInfo::GetThreadId() << ") input matrix L is"
		<< ((L_LoTri) ? "" : " not") << " a lower triangular matrix" << endl;
	cout.flush ();

	Cholesky::trmm4 <DataType, SizeType> (L, L, temp, matrix_size, true);
	
	long double error = 0.0;
	DataType value1, value2;

	for (SizeType row = 0; row < matrix_size; row ++) {
		for (SizeType col = 0; col < matrix_size; col ++) {
			value1 = RowMajorElem (S, row, col, matrix_size);
			value2 = RowMajorElem (temp, row, col, matrix_size);
			error += (value1 - value2) * (value1 - value2);
		}
	}
	// error = error / matrix_elems;
	error = static_cast<long double> (pow ((double) error, 2.0));

	std::cout << "(" << ProcInfo::GetThreadId () << ") tolerance = " << std::fixed << error
		<< " " << ((abs(error) <= epsilon) ? "Accepted" : "Failed");

	int i = 3;
	while ((i--) > 0) { Sleep (1000); std::cout << "."; std::cout.flush (); }
	std::cout << endl;

	Cholesky::sfree <ElemType, SizeType> (temp, matrix_elems);
	return;
}