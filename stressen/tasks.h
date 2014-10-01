/** This file defined all of the tasks **/
#ifndef _TASKS_H
#define _TASKS_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <time.h>

#include "procinfo.h"
#include "strassen.h"
#include "data_module.h"
#include "thread_pool.h"
/**
 Use boost to get a date string in _datestr ()
 **/
#include "boost\date_time\posix_time\posix_time.hpp"
#include "boost\date_time\posix_time\posix_time_io.hpp"

using namespace std;
using namespace boost::posix_time;

enum TaskId {
	MatrixAdd = 1,
	MatrixMinus,
	MatrixMulti,
	ChoPOTRF = 4,				/** R = chol (A) **/
	ChoTRSM,						/** R = A / trans (B) **/
	ChoTRSM2,					/** R = A / B **/
	ChoTRSM3,					/** R = - A / B **/
	ChoSYRK,						/** R = A - B * trans (B) **/
	ChoSYRK2,					/** R = A + trans (B) * B **/
	ChoGEMM,					/** R = A - B * trans (C) **/
	ChoGEMM2,					/** R = A - B * C **/
	ChoGEMM3,					/** R = A + trans (B) * C **/
	ChoTRTRI,						/** R = inv (A) **/
	ChoLAUUM,					/** R = trans (A) * A **/
	ChoTRMM,					/** R = A * B **/
	ChoTRMM2,					/** R = -A * B **/
	ChoTRMM3,					/** R = trans (A) * B **/
	ChoTRMM4,					/** R = A * trans (B) **/
}; 

#include "../cholesky/tasks_cholesky.h"
#include "../cholesky/cholesky.h"

/**************************/
/** R = A + B **/
template <typename T>
void StrassenAdd (DataModule<T> * A, DataModule<T> * B, DataModule<T> * R, bool refresh) 
{
	/** check-up addresses **/
	if (A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt
		<< "] task: (data " << A->id() << ") + (data " << B->id ()
		<< ") = (data " << R->id () << ")" << std::endl;
	cout << ss.str ();  cout.flush ();  ss.str ("");
	/** computing **/
	Strassen::add <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	A->addVisits ();
	B->addVisits ();
	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt
		<< "] task: (data " << R->id() << ") is available" << std::endl;
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** StrassenAdd **/

/** R = A - B **/
template <typename T>
void StrassenMinus (DataModule<T> * A, DataModule<T> * B, DataModule<T> * R, bool refresh) 
{
	/** check-up addresses **/
	if (A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt
		<< "] task: (data " << A->id() << ") - (data " << B->id ()
		<< ") = (data " << R->id () << ")" << std::endl;
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Strassen::minus <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	A->addVisits ();
	B->addVisits ();
	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt
		<< "] task: (data " << R->id() << ") is available" << std::endl;
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** StrassenMinus **/

/** R = A * B **/
template <typename T>
void StrassenMulti (DataModule<T> * A, DataModule<T> * B, DataModule<T> * R, bool refresh)
{
	/** check-up addresses **/
	if (A->data () == NULL)  cout << "Runtime error: invalid address for (data " << A->id () << ")";
	if (B->data () == NULL)  cout << "Runtime error: invalid address for (data " << B->id () << ")";
	if (R->data () == NULL)  cout << "Runtime error: invalid address for (data " << R->id () << ")";
	cout.flush ();

	std::stringstream ss;
	std::string dt;
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt
		<< "] task: (data " << A->id() << ") * (data " << B->id ()
		<< ") = (data " << R->id () << ")" << std::endl;
	cout << ss.str ();  cout.flush ();  ss.str ("");

	/** computing **/
	Strassen::multi <T, size_t> (A->data (), B->data (), R->data (), A->elems (), refresh);
	/** post-processing **/
	A->addVisits ();
	B->addVisits ();
	while (! newAvailModules.push ((int) R->id ()));
	now = boost::posix_time::second_clock::local_time ();
	dt = boost::posix_time::to_simple_string (now);
	ss << "(" << ProcInfo::GetThreadId() << ") [" << dt
		<< "] task: (data " << R->id() << ") is available" << std::endl;
	cout << ss.str ();  cout.flush ();  ss.str ("");
} /** StrassenMulti **/

/**************************/

template <typename DataType>
class TaskMap {

public:
	TaskMap () {}
	~TaskMap () {}

	/** for Strassen algorithm **/
	static void add (
		ThreadPool <DataType> & pool, TaskId task_id,
		DataModule <DataType> * InPut1,
		DataModule <DataType> * InPut2,
		DataModule <DataType> * Result);

	/** for Cholesky algorithm **/
	static void add (
		ThreadPool <DataType> & pool, TaskId task_id,
		DataModule <DataType> * InPut1,
		DataModule <DataType> * InPut2, 
		DataModule <DataType> * InPut3,
		DataModule <DataType> * Result);

	static const char * taskId2Str (TaskId task_id);

}; /** class TaskMap **/

template <typename DataType>
const char * TaskMap <DataType> :: taskId2Str (TaskId task_id)
{
	switch (task_id) {
	case MatrixAdd : return "Matrix Addition";
	case MatrixMinus : return "Matrix Minus";
	case MatrixMulti : return "Matrix Multiplication";
	case ChoPOTRF : return "Cholesky PORTF";
	case ChoTRSM : return "Cholesky TRSM";
	case ChoTRSM2 : return "Cholesky TRSM2";
	case ChoTRSM3 : return "Cholesky TRSM3";
	case ChoSYRK : return "Cholesky SYRK";
	case ChoSYRK2 : return "Cholesky SYRK2";
	case ChoGEMM : return "Cholesky GEMM";
	case ChoGEMM2 : return "Cholesky GEMM2";
	case ChoGEMM3 : return "Cholesky GEMM3";
	case ChoTRTRI : return "Cholesky TRTRI";
	case ChoLAUUM : return "Cholesky LAUUM";
	case ChoTRMM : return "Cholesky TRMM";
	case ChoTRMM2 : return "Cholesky TRMM2";
	case ChoTRMM3 : return "Cholesky TRMM3";
	case ChoTRMM4 : return "Cholesky TRMM4";
	}
	return "not defined operation";
}


template <typename DataType>
void TaskMap <DataType> :: add (
	ThreadPool <DataType> & pool,	TaskId task_id,
	DataModule <DataType> * InPut1, DataModule <DataType> * InPut2,
	DataModule <DataType> * Result ) 
{
	size_t size;

	switch (task_id) {
	case MatrixAdd :
		pool.getPoolPtr ()->schedule (boost::bind (StrassenAdd <DataType>, InPut1, InPut2, Result, true));
		break;
		
	case MatrixMinus :
		pool.getPoolPtr ()->schedule (boost::bind (StrassenMinus <DataType>, InPut1, InPut2, Result, true));
		break;

	case MatrixMulti :
		pool.getPoolPtr ()->schedule (boost::bind (StrassenMulti <DataType>, InPut1, InPut2, Result, true));
		break;

	default:
		cout << ProcInfo () << " : unknown task_id = " << task_id << endl;
		cout.flush ();
	}
} /** TaskMap <DataType> :: add **/


template <typename DataType>
void TaskMap <DataType> :: add (
	ThreadPool <DataType> & pool,	TaskId task_id,
	DataModule <DataType> * InPut1,
	DataModule <DataType> * InPut2,
	DataModule <DataType> * InPut3,
	DataModule <DataType> * Result ) 
{
	size_t size;

	switch (task_id) {
	case ChoPOTRF :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyPOTRF <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoTRSM :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRSM <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoTRSM2 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRSM2 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;
		
	case ChoTRSM3 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRSM3 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoSYRK :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskySYRK <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoSYRK2 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskySYRK2 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoGEMM :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyGEMM <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoGEMM2 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyGEMM2 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoGEMM3 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyGEMM3 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoTRTRI :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRTRI <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoLAUUM :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyLAUUM <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoTRMM :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRMM <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;
		
	case ChoTRMM2 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRMM2 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	case ChoTRMM3 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRMM3 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	
	case ChoTRMM4 :
		pool.getPoolPtr ()->schedule (boost::bind (CholeskyTRMM4 <DataType>, InPut1, InPut2, InPut3, Result, true));
		break;

	default:
		cout << ProcInfo () << " : unknown task_id = " << task_id << endl;
		cout.flush ();
	}
} /** TaskMap <DataType> :: add **/


/**************************/

#endif