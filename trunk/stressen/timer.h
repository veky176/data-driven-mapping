#ifndef _TIMER_H
#define _TIMER_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>

#include "procinfo.h"
#include "boost\date_time\posix_time\posix_time.hpp"

using namespace std;


class SimpleTimer {
private:
	// time_t prev_time;
	boost::posix_time::ptime prev_time;
	std::string _name;

public:
	SimpleTimer (const char * name) : _name (name) 
	{
		// std::time (& prev_time);
		prev_time = boost::posix_time::microsec_clock::local_time ();;
	}
	~SimpleTimer () {
		std::cout << (* this);
		_name.clear ();
	}
	
	double elapsed_seconds () {
		boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time ();
		boost::posix_time::time_duration diff = now - prev_time;
		double millisec = diff.total_milliseconds ();
		return (millisec / 1000.0);
		// return std::difftime (std::time (NULL), prev_time);
	}

	friend std::ostream& operator<< (std::ostream & os, SimpleTimer & object)
	{
		std::stringstream ss;
		ss << "(" << ProcInfo::GetThreadId () << ") "
			<< object._name << " consumed " << std::fixed << std::setprecision (3)
			<< object.elapsed_seconds ()
			<< " seconds" << std::endl;
		os << ss.str ();
		os.flush ();
		ss.str ("");
		return os;
	}
}; /** SimpleTimer **/


#endif