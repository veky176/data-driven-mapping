/**
 This file defined a mutex lock class.
 **/
#ifndef __MUTEXLOCK_MACRO
#define __MUTEXLOCK_MACRO

#include "boost\thread\mutex.hpp"
#include "boost\thread\thread.hpp"

using namespace std;

/*********************************************************/

class SimpleMutexLock {
private:
	boost::mutex mtx_;
	
public:
	SimpleMutexLock () {}
	virtual ~SimpleMutexLock () {}
	
	bool try_lock () { return mtx_.try_lock (); }
	bool try_lock (unsigned int timedout_millisecond) {
		unsigned int elapsed_time = 0;
		bool result = try_lock ();
		while (!result && (elapsed_time++ < timedout_millisecond)) {
			boost::this_thread::sleep (boost::posix_time::milliseconds (1));
			result = try_lock ();
		}
		return result;
	} /// try_lock ()			
	
	void lock () { mtx_.lock (); }
	void unlock () { mtx_.unlock (); }
	
}; // class SimpleMutexLock

/*********************************************************/
#endif