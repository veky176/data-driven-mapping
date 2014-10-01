#ifndef _DATA_MODULE_H
#define _DATA_MODULE_H
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>

#include "procinfo.h"
#include "timer.h"
#include "mutexlock.h"


using namespace std;

template <typename DataType>
class DataModule : public SimpleMutexLock
{
private:
	int _maxVisits; /** maximum number of use, once the use times reached to the threshold, the instance will be removed **/
	int _visits; /** how many times the module was used **/
	int module_id;
	size_t _elems; /** numer of elems in the module, not memory size **/
	DataType * _data;

	void _allocate ();

public:
	DataModule (int myid, size_t myElems, int maxVisits);
	~DataModule () {
		this->clear ();
	} /** ~DataModule **/

	size_t size () { return _elems * _elems * sizeof (DataType); } /** memory size **/
	size_t elems () { return _elems; } /** number of elements in one dimension **/

	DataType * data () { 
		if (_data == NULL)  _allocate ();
		return _data; 
	}

	int id () { return module_id; }
	DataModule * getPointer () { return this; }

	void clear (bool thread_safe = true) {
		if (thread_safe)  this->lock ();
		if (_data != NULL) {
			std::stringstream ss;
			std::string dt;
			boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
			dt = boost::posix_time::to_simple_string (now);
			ss << "(" << ProcInfo::GetThreadId() << ") [" << dt << "] erase (data " << module_id << ") ";
			Strassen::sfree <DataType, size_t> (_data, _elems * _elems);
			_data = NULL;
			ss << "current memory usage = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
			std::cout << ss.str ();
			std::cout.flush ();
			ss.str ("");
		}
		if (thread_safe)  this->unlock ();
	} /** clear **/

	void setVisits (int vis) {
		std::stringstream ss;
		this->lock ();
		_maxVisits = vis;
		ss << "(" << ProcInfo::GetThreadId () << ") (data " << module_id
			<< ") set max_visits = " << _maxVisits << endl;
		this->unlock ();
		cout << ss.str (); cout.flush (); ss.str ("");
	}
	void addVisits () { 
		this->lock (); 
		_visits += 1;
		if (_maxVisits > 0 && _visits >= _maxVisits)
			this->clear (false);
		this->unlock (); 
	}
	int getVisits () { return _visits; }
	int getMaxVisits () { return _maxVisits; }


}; /** DataModule **/

template <typename DataType>
DataModule<DataType> :: DataModule 
	(int myid, size_t myElems, int maxVisits)
	: module_id (myid), _elems (myElems), _data (NULL), _maxVisits (maxVisits), _visits (0)
{
} /** DataModule :: DataModule **/


template <typename DataType>
void DataModule<DataType> :: _allocate ()
{
	this->lock ();
	if (_data != NULL) { this->unlock (); return; }
	_data = Strassen::scalloc <DataType, size_t> (_elems * _elems);
	if (_data == NULL)  {
		cout << ProcInfo() << ": failed to allocate the memory" << endl;
		cout.flush ();
	} else {
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId () << ") allocated (data "
			<< this->id () << ") elems = " << this->elems ()
			<< ", current memory size = " << PrintBytes (Strassen::_nowMemoryUse) << endl;
		cout << ss.str (); cout.flush (); ss.str ("");
	}
	this->unlock ();
} /** DataModule :: DataModule **/

/********************************/

template <typename DataType>
class DataModuleContainer : public SimpleMutexLock 
{
	typedef DataModule <DataType> ElemType;
	typedef ElemType * ElemPtr;
	typedef std::map <int, ElemPtr> Mapper;

private:
	Mapper _buffer;

	void _clear_all (); /** Thread Safe **/

public:
	DataModuleContainer () {}
	~DataModuleContainer () { _clear_all (); }

	bool add (ElemPtr ptr, bool thread_safe = true) {
		if (ptr == NULL) return false;
		if (thread_safe)  this->lock ();
		_buffer [ptr->id ()] = ptr->getPointer ();
		cout << "(" << ProcInfo::GetThreadId () << ") added data module (id = "
			<< ptr->id () << ", elems = " << ptr->elems ()
			<< ", visits = " << ptr->getMaxVisits ()
			<< ", ptr = " << ptr << ")" << endl;
		cout.flush ();
		if (thread_safe)  this->unlock ();
		return find (ptr->id (), thread_safe);
	} /** add **/

	bool erase (int id, bool thread_safe = true) {
		if (thread_safe)  this->lock ();
		Mapper::iterator iter = _buffer.find (id);
		if (iter == _buffer.end ()) {
			if (thread_safe) this->unlock ();
			return false; /** not found the element **/
		}
		try { 
			_buffer.erase (iter); 
		} catch (std::exception & e) {
			cout << "Runtime error at erasing a map element (code = " << e.what () << ")" << endl;
			cout.flush (); 
			if (thread_safe) this->unlock ();
			return false; /** deletion error **/
		}
		if (thread_safe) this->unlock ();
		return (! find (id, thread_safe)); /** job well done **/
	} /** erase **/

	bool find (int id, bool thread_safe = false) {
		if (thread_safe) {
			this->lock ();
			bool ret = (_buffer.find (id) != _buffer.end ()); 
			this->unlock ();
			return ret;
		}
		return (_buffer.find (id) != _buffer.end ()); 
	}
	bool find (int id, ElemPtr &ptr, bool thread_safe = false) {
		if (thread_safe)  this->lock ();
		Mapper::iterator iter = _buffer.find (id);
		if (iter == _buffer.end ()) {
			if (thread_safe)  this->unlock ();
			return false;
		}
		ptr = iter->second;
		if (thread_safe)  this->unlock ();
		return true;
	} /** find **/

	size_t size (bool thread_safe = false) { return _buffer.size ();  }
	bool empty () { return _buffer.empty (); }
	void clear () { _clear_all (); }

	/** add an element **/
	void operator+= (ElemPtr elem) { return this->add (elem); } /** operator+= **/
	/** delete an element **/
	void operator -= (int elemId) { this->erase (elemId); } /** operator-= **/

	/** access one element **/
	ElemPtr operator[] (int elemId) {
		ElemPtr ptr = NULL;
		if (find (elemId, ptr, false))  return ptr;
		return NULL;
	}

}; /** DataModuleContainer **/

template <typename DataType>
void DataModuleContainer <DataType> :: _clear_all ()
{
	if (empty ()) return;
	this->lock ();
	while (! _buffer.empty ()) {
		ElemPtr elem = _buffer.begin()->second;
#if defined DEBUG && DEBUG > 2
		cout << "(" << ProcInfo::GetThreadId () << ") delete (data "
			<< elem->id () << ") size = " << elem->size () << ", ptr = " << elem << endl;
		cout.flush ();
#endif
		if (elem != NULL) delete elem;
		_buffer.erase (_buffer.begin ());
	}
	this->unlock ();
} /** DataModuleContainer :: _clear_all **/


#endif /** #ifndef _DATA_MODULE_H **/