#ifndef _PROCINFO_H
#define _PROCINFO_H
#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
/**
 Boost.Process is not (yet) an official Boost C++ library.
 It must be downloaded and installed separately.
 As Boost.Process is header-only the two directories boost and libs in the zip archive
 only need to be copied to your Boost directory.
 As Boost.Process has been tested with Boost 1.36.0 this is the minimum version supported.
 **/
#include "boost\process.hpp"
#include "boost\log\attributes\current_process_id.hpp"
#include "boost\log\attributes\current_process_name.hpp"
#include "boost\log\attributes\current_thread_id.hpp"
#include "boost\thread.hpp"

using namespace std;
using namespace boost::process;
using namespace boost::this_thread;

/** get the process and thread info **/
class ProcInfo {
	
public:
	unsigned long proc_id;  /** process ID **/
	unsigned long pproc_id;  /** parent process ID **/
	unsigned long thread_id;  /** thread ID **/
	std::string proc_name;  /** processor name **/
	
	ProcInfo ();
	virtual ~ProcInfo ();
	
	/** clear up the memory **/
	virtual void clear ();
	
	/** get all attributes at once **/
	void GetAttrs ();
	
	/** get process ID **/
	static unsigned long GetProcId ();
	
	/** get parent process ID **/
	static unsigned long GetPProcId ();
	
	/** get processor/CPU name as a string **/
	static const char * GetProcName (std::string& name);
	static const char * GetCPUName (std::string& name);
	
	/** get the current thread ID **/
	static unsigned long GetThreadId (int base = 16);
	
	/** output the CPU, process and thread infor... **/
	friend std::ostream& operator<< (std::ostream&, ProcInfo& object);
	
}; //  ProcInfo


#endif