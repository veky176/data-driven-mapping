#ifndef _DATA_MAPPING_H
#define _DATA_MAPPING_H

#include "tasks.h"
#include <fstream>

using namespace std;

/** one item for the data dependence matrix **/
struct DataDepElem {
	TaskId task_id;
	int * input; /** pointer to the input argument array **/
	int _argc; /** number of input arguments **/
	int output;
	DataDepElem (TaskId id, int argc, int * argv, int ret)
		: task_id (id), output (ret), _argc (argc), input (NULL) {
			input = (int *) calloc (argc, sizeof  (int));
			if (argv != NULL) 
				for (size_t i = 0; i < _argc; i++) 
					input [i] = argv [i];
	}
	~DataDepElem () { if (input != NULL) free (input); }

	/** return a copy of itself **/
	DataDepElem * copy () {
		return (new DataDepElem (this->task_id, this->_argc, this->input, this->output));
	} /** copy **/

	DataDepElem & operator= (const DataDepElem & other) {
		if (this != &other) {
			task_id = other.task_id;
			_argc = other._argc;
			output = other.output;
			for (size_t i = 0; i < _argc; i++) input [i] = other.input [i];
		}
	} /** instance copy **/

	bool operator== (const DataDepElem & other) {
		if (this == &other) return true;
		if (task_id != other.task_id) return false;
		if (_argc != other._argc) return false;
		if (output != other.output)  return false;
		for (size_t i = 0; i < _argc; i++)
			if (input [i] != other.input [i])  return false;
		return true;
	} /** instance compare **/

	friend std::ostream & operator<< (std::ostream & os, DataDepElem & obj)
	{
		stringstream ss;
		ss << "TASK (" << TaskMap<int>::taskId2Str (obj.task_id) << ") inputs = [";
		if (obj._argc > 0)  ss << obj.input [0];
		for (int i = 1; i < obj._argc - 1; i ++) {
			if (obj.input [i] <= 0)  continue; /** task id > 0 **/
			ss << " " << obj.input [i];
		}
		if (obj.input [obj._argc - 1] > 0)  ss << " " << obj.input [obj._argc - 1];		
		ss << "] output = [" << obj.output << "]" << endl;
		os << ss.str ();
		os.flush ();
		ss.str ("");
		return os;
	} /** operator<< **/
}; 


class DataMapping : public SimpleMutexLock
{
private:
	/** data dependency map **/
	map<int, DataDepElem *> _dep;
	
	int _min_task_id;
	int _max_task_id;

	map<int, int> _visits;
	
	bool * avail_modules;
	int _max_size; /** max number of modules **/

	void _read (const char * filename, int argc = 2);
	void _visits_info (int default_value = -1); /** build _visits **/
	
	bool _task_avail (DataDepElem *); /** test if a task is available to go **/

public:
	DataMapping (const char * filename, int input_count, int max_size)
		: _max_size (max_size), avail_modules (NULL),
		_min_task_id (1), _max_task_id (1)
	{
		avail_modules = (bool *) calloc (_max_size, sizeof (bool));
		for (size_t each = 0; each < _max_size; each ++)  avail_modules [each] = false;
		_read (filename, input_count);
		_visits_info ();
	}
	~DataMapping () {
		if (avail_modules != NULL) free (avail_modules);
		_dep.clear ();
	}

	int maxTaskId () { return this->_max_task_id; }
	int minTaskId () { return this->_min_task_id; }

	size_t size () { 
		this->lock ();
		size_t sz = _dep.size (); 
		this->unlock ();
		return sz;
	}
	bool empty () {
		this->lock ();
		bool ret = _dep.empty (); 
		this->unlock ();
		return ret;
	}

	void setAvailModule (int module_id) {
		this->lock ();
		if (module_id >= _max_size) {
			cout << "module id " << module_id << " is out of scope, max_size = "
				<< _max_size << endl; cout.flush ();
			this->unlock ();
			return;
		}
		avail_modules [module_id] = true;
		this->unlock ();
		std::stringstream ss;
		std::string dt;
		boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
		dt = boost::posix_time::to_simple_string (now);
		ss << "(" << ProcInfo::GetThreadId () << ") [" << dt << "] (data " << module_id
			<< ") is available" << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
	} /** setAvailModule **/

	bool checkAvailModule (int module_id, bool thread_safe = true) {
		if (thread_safe) this->lock ();
		if (module_id >= _max_size) {
			cout << "module id " << module_id << " is out of scope, max_size = "
				<< _max_size << endl; cout.flush ();
			if (thread_safe) this->unlock ();
			return false;
		}
		bool ret = avail_modules [module_id];
		if (thread_safe) this->unlock ();
		return ret;
	} /** checkAvailModule **/

	void setModuleVisits (int id, int default_value = -1) {
		this->lock ();
		_visits [id] = default_value;
		this->unlock ();
	} /** set a visit number of a data module **/

	int getModuleVisits (int id, int default_value = -1) {
		this->lock ();
		int vis = default_value;
		if (_visits.find (id) != _visits.end ())  vis = _visits [id];
		this->unlock ();
		return vis;
	} /** return visits of a data module **/

	/** search for the next task **/
	DataDepElem *  next_task ();

}; /** class DataMapping **/

#endif /** _DATA_MAPPING_H **/