#include "utility.h"
#include "procinfo.h"

#include "boost\date_time\posix_time\posix_time.hpp"

using namespace std;

/** return the time in milliseconds **/
long getCurrentTime ()
{
	static boost::posix_time::ptime starter = boost::posix_time::microsec_clock::local_time ();
	boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time ();
//	static boost::posix_time::ptime starter = boost::posix_time::second_clock::local_time ();
//	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time ();
	boost::posix_time::time_duration diff = now - starter;
	return static_cast <long> (diff.total_milliseconds ());
}

void press_key ()
{
#if defined WIN32
	cout << "press any key to continue... " << endl;
	cout.flush ();
	getchar ();
	cout << "keep running" << endl;
	cout.flush ();
#endif
} /** press_key **/


void writeLog (std::stringstream & ss, LogFile & log)
{
	if (ss.str ().empty ()) return;
	std::stringstream _ss;
	_ss << "(" << ProcInfo::GetThreadId() << ") [" << getCurrentTime ()
		<< "] " << ss.str () << endl;
	log.write (_ss.str ());
	_ss.str ("");
}

void writeLog (std::string & ss, LogFile & log)
{
	if (ss.empty ()) return;
	std::stringstream _ss;
	_ss << "(" << ProcInfo::GetThreadId() << ") [" << getCurrentTime ()
		<< "] " << ss << endl;
	log.write (_ss.str ());
	_ss.str ("");
}

void writeLog (const char * ss, LogFile & log)
{
	if (ss == NULL || strlen (ss) == 0) return;
	std::stringstream _ss;
	_ss << "(" << ProcInfo::GetThreadId() << ") [" << getCurrentTime ()
		<< "] " << ss << endl;
	log.write (_ss.str ());
	_ss.str ("");
}

void writeLog (size_t size, LogFile & log)
{
	std::stringstream _ss;
	_ss << "(" << ProcInfo::GetThreadId() << ") [" << getCurrentTime ()
		<< "] " << size << endl;
	log.write (_ss.str ());
	_ss.str ("");
}

void writeLog (int size, LogFile & log)
{
	std::stringstream _ss;
	_ss << "(" << ProcInfo::GetThreadId() << ") [" << getCurrentTime ()
		<< "] " << size << endl;
	log.write (_ss.str ());
	_ss.str ("");
}

void writeLog (double number, LogFile & log)
{
	std::stringstream _ss;
	_ss << "(" << ProcInfo::GetThreadId() << ") [" << getCurrentTime ()
		<< "] " << std::fixed << number << endl;
	log.write (_ss.str ());
	_ss.str ("");
}

