#ifndef _TASKS_CHOLESKY_H
#define _TASKS_CHOLESKY_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <time.h>

#include "procinfo.h"
#include "../cholesky/cholesky.h"
#include "data_module.h"
#include "thread_pool.h"
/**
 Use boost to get a date string in _datestr ()
 **/
#include "boost\date_time\posix_time\posix_time.hpp"
#include "boost\date_time\posix_time\posix_time_io.hpp"
#include "boost\lockfree\queue.hpp"
#include "boost\lockfree\stack.hpp"

using namespace std;
using namespace boost::posix_time;

extern boost::lockfree::queue <int> newAvailModules;

/**************************/
template <typename T>
void CholeskyPOTRF (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 1;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: " 	<< "(data " << R->id() 
		<< ") = chol (data " << A->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::potrf <T, size_t> (A->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyPOTRF **/


/** TRSM : R = A / trans(B) **/
template <typename T>
void CholeskyTRSM (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") = (data " << A->id () 
		<< ") / trans (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trsm <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();
	
	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRSM **/


/** TRSM2 : R = A / B **/
template <typename T>
void CholeskyTRSM2 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ")  = (data " << A->id () 
		<< ") / (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trsm2 <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRSM2 **/


/** TRSM3 : R = - A / B **/
template <typename T>
void CholeskyTRSM3 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ")  = - (data " << A->id () 
		<< ") / (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trsm3 <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRSM3 **/


/** SYRK : R = A - B * trans (B) **/
template <typename T>
void CholeskySYRK (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: "
		<< "(data " << R->id() << ") = (data " << A->id () << ") - (data " << B->id () 
		<< ") * trans (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::syrk <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskySYRK **/


/** SYRK2 : R = A + trans (B) * B **/
template <typename T>
void CholeskySYRK2 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C, DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: "
		<< "(data " << R->id() << ") = (data " << A->id () << ") + trans (data " << B->id () 
		<< ") * (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::syrk2 <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskySYRK2 **/


/** GEMM : R = A - B * trans (C) **/
template <typename T>
void CholeskyGEMM (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 3;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: "
		<< "(data " << R->id() << ") = (data " << A->id () << ") - (data " << B->id () 
		<< ") * trans (data " << C->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::gemm <T, size_t> (A->data (), B->data (), C->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyGEMM **/


/** GEMM2 : R = A - B * C **/
template <typename T>
void CholeskyGEMM2 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C, DataModule<T> * R, bool refresh) 
{
	int inputs_n = 3;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: "
		<< "(data " << R->id() << ") = (data " << A->id () << ") - (data " << B->id () << ") * (data " << C->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::gemm2 <T, size_t> (A->data (), B->data (), C->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyGEMM2 **/


/** GEMM3 : R = A + trans (B) * C **/
template <typename T>
void CholeskyGEMM3 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C, DataModule<T> * R, bool refresh) 
{
	int inputs_n = 3;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: "
		<< "(data " << R->id() << ") = (data " << A->id () << ") + trans (data " << B->id () << ") * (data " << C->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::gemm3 <T, size_t> (A->data (), B->data (), C->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyGEMM3 **/


/** TRTRI : R = inv (A) **/
template <typename T>
void CholeskyTRTRI (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 1;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") = inv (data " << A->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trtri <T, size_t> (A->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRTRI **/


/** LAUUM : R = trans (A) * A **/
template <typename T>
void CholeskyLAUUM (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 1;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") = trans (data " 
		<< A->id () << ") * (data " << A->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::lauum <T, size_t> (A->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyLAUUM **/


/** triangular matrix-matrix multiply : TRMM : R = A * B **/
template <typename T>
void CholeskyTRMM (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id()
		<< ") = (data " << A->id () << ") * (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trmm <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRMM **/


/** triangular matrix-matrix multiply : TRMM : R = - A * B **/
template <typename T>
void CholeskyTRMM2 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id()
		<< ") = - (data " << A->id () << ") * (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trmm2 <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRMM2 **/


/** triangular matrix-matrix multiply : TRMM : R = trans (A) * B **/
template <typename T>
void CholeskyTRMM3 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id()
		<< ") = trans (data " << A->id () << ") * (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trmm3 <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRMM3 **/


/** triangular matrix-matrix multiply : TRMM : R = A * trans (B) **/
template <typename T>
void CholeskyTRMM4 (DataModule<T> * A, DataModule<T> * B, DataModule<T> * C,
	DataModule<T> * R, bool refresh) 
{
	int inputs_n = 2;
	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);

	/** check-up addresses **/
	if (inputs_n >= 1 && A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (inputs_n >= 2 && B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (inputs_n >= 3 && C->data () == NULL)  cout << "Runtime error: invalid address for (data " << C->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id()
		<< ") = (data " << A->id () << ") * trans (data " << B->id () << ")" << std::endl;
	logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Cholesky::trmm4 <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	if (inputs_n >= 1)  A->addVisits ();
	if (inputs_n >= 2)  B->addVisits ();
	if (inputs_n >= 3)  C->addVisits ();

	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] task: (data " << R->id() << ") is available" << std::endl;
	// logger.write (ss);
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** CholeskyTRMM4 **/

/**************************/

#endif