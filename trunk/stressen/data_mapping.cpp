#include "data_mapping.h"

using namespace std;

#define READ_VISITS_FILE 0
// #define INS_BARRIER

void DataMapping :: _read 
	(const char * filename, int argc)
{
	cout << "(" << ProcInfo::GetThreadId () 
		<< ") reading " << filename << endl;
	cout.flush ();
	std::ifstream file (filename);

	while (file.is_open () && file.good ()) {
		string line, msg;
		getline (file, line);
		while (line.c_str() [0] == ' ') line.erase (0, 1);

		msg = line.substr (0, line.find_first_of (' '));
		if (msg.empty () || msg.c_str() [0] == '#') {
			/** this is an empty or comment line **/
			continue;
		}
		TaskId task_id = static_cast <TaskId> (atoi (msg.c_str()));
		line.erase (0, msg.size ());
		while (line.c_str() [0] == ' ') line.erase (0, 1);

		msg = line.substr (0, line.find_first_of (' '));
		int * inputs = (int *) calloc (argc, sizeof (int));

		for (int i = 0; i < argc; i++) {
			msg = line.substr (0, line.find_first_of (' '));
			if (msg.empty ()) {
				cout << "Error: not enough parameters [" << msg << "]" << endl;
				cout.flush ();
				return;
			}
			inputs [i] = static_cast <int> (atoi (msg.c_str ()));
			line.erase (0, msg.size ());
			while (line.c_str() [0] == ' ') line.erase (0, 1);
		}

		msg = line.substr (0, line.find_first_of (' ')); 
		if (msg.empty ()) {
			cout << "Error: not enough parameters [" << msg << "]" << endl;
			cout.flush ();
			return;
		}
		int output = static_cast <int> (atoi (msg.c_str ()));
		
		int id = _dep.size ();
		DataDepElem * elem = new DataDepElem (task_id, argc, inputs, output);

		_dep.insert (std::pair <int, DataDepElem *> (id, elem));

		std::cout << (* _dep [id]);

		if (inputs != NULL)  free (inputs);
	}

	if (file.is_open ()) file.close ();
	else {
		cout << "Error: cannot open " << filename << endl;
		cout.flush ();
	}
	
#if READ_VISITS_FILE
	/** read the visits file **/
	cout << "(" << ProcInfo::GetThreadId () 
		<< ") reading " << filename_visits << endl;
	cout.flush ();

	file.open (filename_visits);
	
	while (file.is_open () && file.good ()) {
		string line, msg;
		getline (file, line);
		while (line.c_str() [0] == ' ') line.erase (0, 1);

		int id, vis = -1;

		msg = line.substr (0, line.find_first_of (' '));
		if (msg.empty ()) {
			cout << "Error: not enough parameters" << endl;
			cout.flush ();
			return;
		}
		id = static_cast <int> (atoi (msg.c_str()));

		line.erase (0, msg.size ());
		while (line.c_str() [0] == ' ') line.erase (0, 1);
		msg = line.substr (0, line.find_first_of (' '));
		if (msg.empty ()) {
			cout << "Error: not enough parameters" << endl;
			cout.flush ();
			return;
		}
		vis = static_cast <int> (atoi (msg.c_str()));

		this->_visits [id] = vis;
		cout << "Max. visit for module " << id << " = " << _visits [id] << endl;
		cout.flush ();
	}
	
	if (file.is_open ()) file.close ();
	else {
		cout << "Error: cannot open " << filename_visits << endl;
		cout.flush ();
	}
#endif
}


void DataMapping :: _visits_info (int default_value)
{
	map<int, DataDepElem *>::iterator iter;

	for (iter = _dep.begin (); iter != _dep.end (); iter ++) {
		_min_task_id = min (iter->second->output, _min_task_id);
		_max_task_id = max (iter->second->output, _max_task_id);
		for (int i = 0; i < iter->second->_argc; i++) {
			if (iter->second->input [i] <= 0) continue;
			_min_task_id = min (iter->second->input [i], _min_task_id);
			_max_task_id = max (iter->second->input [i], _max_task_id);
		}
	}

	cout << "min. task id = " << _min_task_id << endl;
	cout << "max. task id = " << _max_task_id << endl;
	cout.flush ();

	for (int id = _min_task_id; id <= _max_task_id; id++)
		_visits [id] = default_value;

	for (iter = _dep.begin (); iter != _dep.end (); iter ++) {
		for (int i = 0; i < iter->second->_argc; i++) {
			if (iter->second->input [i] <= 0) continue; /** task id > 0 **/
			int id = iter->second->input [i];
			if (_visits [id] == default_value)
				_visits [id] = 1;
			else
				_visits [id] += 1;
		}
	}

	for (map<int, int>::iterator it = _visits.begin (); it != _visits.end (); it ++) {
		cout << "max visit for (data " << it->first << ") = " << _visits [it->first] << endl;
	}
	cout.flush ();
}


/** search for the next task **/
DataDepElem * DataMapping :: next_task ()
{
	int key = -1;
	DataDepElem * task = NULL;
	this->lock ();
	if (_dep.empty ()) {
		this->unlock ();
		return NULL;
	}

	map<int, DataDepElem *>::iterator iter;

	for (iter = _dep.begin (); iter != _dep.end () && task == NULL; iter ++) 
	{
		task = iter->second;
		if (_task_avail (task)) {
			key = iter->first; 
		} else {
			task = NULL; /** go to the next **/ 
		}
	} // for-loop

	if (task != NULL) {
		_dep.erase (key);
	}

	this->unlock ();
	return task;
} /** DataMapping :: next_task **/

/** verify if a task is available **/
bool DataMapping :: _task_avail (DataDepElem * task)
{
	if (task == NULL) return false;
	for (int i = 0; i < task->_argc; i++) {
		if (task->input [i] <= 0)  continue; /** task id > 0 **/
		if (! checkAvailModule (task->input [i], false))
			return false;
	}

#ifdef INS_BARRIER
	if (task->output >= 31 && task->output <= 50) 
	{  /** make sure modules [1, 30] all done **/
		for (int id = 1; id <= 30; id ++)
			if (! checkAvailModule (id, false)) return false;
	}
	if (task->output >= 51 && task->output <= 70) 
	{  /** make sure modules [31, 50] all done **/
		for (int id = 31; id <= 50; id ++)
			if (! checkAvailModule (id, false)) return false;
	}
#endif

	std::stringstream ss;
	ss << "(" << ProcInfo::GetThreadId () << ") available task: " << (* task);
	cout << ss.str ();
	cout.flush ();
	ss.str ("");
	return true;
}