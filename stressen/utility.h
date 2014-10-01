#ifndef _UTILITY_H
#define _UTILITY_H

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "mutexlock.h"

using namespace std;

/** return the time in milliseconds **/
long getCurrentTime ();

class LogFile : public SimpleMutexLock
{
	int _buf_size;
	std::string _buffer;
	std::string _filename;
	std::ofstream _file;

	void _open () {
		this->lock ();
		_file.open (_filename.c_str());
		if (! _file.is_open ()) {
			cout << "Runtime Error: cannot open the file " << _filename.c_str() << endl;
			cout.flush ();
		} else {
			cout << "INFO: logfile = " << _filename.c_str () << endl;
			cout.flush ();
		}
		this->unlock ();
	}

	void _close () { 
		this->lock ();
		if (_file.is_open ()) {
			_dump ();
			_file.close (); 
			cout << "INFO: closed " << _filename.c_str() << endl;
			cout.flush ();
		}
		this->unlock ();
	}

	void _dump () {
		if (_file.is_open ()) { _file << _buffer; _file.flush (); }
		_buffer.clear ();
	}

	void _write (std::string &s) {
		this->lock ();
		_buffer += s;
		if (_buffer.size () > _buf_size) _dump ();
		this->unlock ();
	}

public:
	LogFile (const char * filename, int bufsize = 1000000)
		: _filename (filename), _buf_size (bufsize)
	{
			_open (); 
	}

	~LogFile () 
	{ 
		_close (); 
	}

	void write (std::stringstream & ss) {
		_write (ss.str ());
	}
	void write (std::string &s) {
		_write (s);
	}
	void write (const char * str) {
		_write (std::string (str));
	}

}; 


template <typename DataType, typename SizeType> 
void genRandomArray (DataType *, SizeType, DataType, DataType);

template <typename DataType, typename SizeType> 
void genRandomArray (DataType * _array, SizeType _size, DataType min, DataType max)
{
	SimpleTimer timer (__FUNCTION__);

	DataType random;
	srand (time(NULL));
	for (SizeType elem = 0; elem < _size; elem ++) {
		if (elem % 10 == 0) srand (time(NULL));
		random = ((DataType) rand()) / (DataType) RAND_MAX;
		_array [elem] = min + random * (max - min);
	}
} /** genRandomArray **/

template <typename DataType>
void readMatrixFromFile (DataType * _array, int _size, const char * filename)
{
	SimpleTimer timer (__FUNCTION__);
	std::ifstream file;

	file.open (filename);

	if (! file.is_open ()) {
		cout << "File open error: cannot open " << filename << endl;
		cout.flush ();
		return;
	}

	cout << "(" << ProcInfo::GetThreadId () << ") reading " << filename << endl;
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
} /** readMatrixFromFile **/


template <typename DataType, typename SizeType>
bool is2ArrayEqual (DataType * ptr1, DataType * ptr2, SizeType elems, DataType eps)
{
	DataType diff = 0;
	for (SizeType i = 0; i < elems; i ++) {
		diff = static_cast <DataType> (abs (ptr1 [i] - ptr2 [i]));
		if (diff > eps) return false;
	}
	return true;
} /** is2ArrayEqual **/

void press_key ();

void writeLog (std::stringstream & ss, LogFile & log);
void writeLog (std::string & ss, LogFile & log);
void writeLog (const char * ss, LogFile & log);
void writeLog (size_t size, LogFile & log);
void writeLog (int size, LogFile & log);
void writeLog (double number, LogFile & log);


extern LogFile logger;
extern LogFile errlog;
extern LogFile memlog;

#endif