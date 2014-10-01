#ifndef _THREAD_POOL_H
#define _THREAD_POOL_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <deque>

#include "boost\version.hpp"
/** change TIME_UTC to TIME_UTC_ for new versions... **/
#if BOOST_VERSION > 105000
#define TIME_UTC TIME_UTC_
#endif

#include "threadpool\boost\threadpool.hpp"


using namespace std;
using namespace boost::threadpool;

template <typename DataModulePtr>
class ThreadPool {

private:
	int _maxThreads;
	pool * _taskpool;

	/** init the pool **/
	void _init (int val = 0) {
		if (this->_taskpool == NULL)  this->_taskpool = new pool (val);
	}
	/** delete pool **/
	void _del_pool () {
		if (this->_taskpool != NULL)  delete _taskpool;
		_taskpool = NULL;
	}
	/** wait for all tasks to finish **/
	void _wait_all () {
		this->_taskpool->wait ();
	}
	/** if there is any task still in the pool **/
	bool _empty_pool () {
		return this->_taskpool->empty ();
	}

public:
	ThreadPool (int max_threads) : _maxThreads (max_threads), _taskpool (NULL)
	{
		this->_init (0);
	}
	virtual ~ThreadPool () 
	{
		this->_wait_all ();
		this->_empty_pool ();
		this->_del_pool ();
	}

	void run () {
		cout << ProcInfo () << " running " << _maxThreads << " threads" << endl;
		cout.flush ();
		this->_taskpool->size_controller().resize (_maxThreads);
	}

	/** Returns the number of tasks which are ready for execution. **/
	int pending () { return this->_taskpool->pending (); }

	void wait () { this->_wait_all (); }
	bool empty () { return this->_empty_pool () }

	pool * getPoolPtr () { return this->_taskpool; }

}; /** class ThreadPool **/

#endif