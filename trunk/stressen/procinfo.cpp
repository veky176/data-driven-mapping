#include "procinfo.h"

using namespace std;
using namespace boost::process;
using namespace boost::this_thread;

#ifdef WIN32
#include <Windows.h>
#include <TlHelp32.h>
namespace win32 {
	
	/** get the processor name **/
	void GetProcessorName (std::string& name, int& namelen, int MaxNameLen = 1024) 
	{
		int ret_code = SOCKET_ERROR;
		char * _ProcName = (char *) calloc (MaxNameLen, sizeof (char));
#pragma comment(lib, "ws2_32.lib")
		WSAData wsa_data;
		WSAStartup(MAKEWORD(2, 2), &wsa_data);
		ret_code = gethostname(_ProcName, MaxNameLen);
		name = (ret_code != SOCKET_ERROR) ? (_ProcName) : ("unknown");
		namelen = name.size ();
		WSACleanup();
		if (_ProcName != NULL)  delete _ProcName;
	} // GetProcessorName

	
	/** get the process ID of its parent **/
	unsigned long GetParentProcessID (unsigned long pid)
	{
		unsigned long ppid = 0;
		HANDLE h = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);
		PROCESSENTRY32 pe = { 0 };
		pe.dwSize = sizeof(PROCESSENTRY32);
		
		if (pid == 0)
			pid = GetCurrentProcessId ();
		
		if( Process32First(h, &pe)) {
			do {
				if (pe.th32ProcessID == pid) {
					char idstr [128];
					sprintf (idstr, "%i", pe.th32ParentProcessID);
					ppid = strtoul (idstr, NULL, 0);
				}
			} while(Process32Next(h, &pe));
		}
		
		CloseHandle(h);
		
		return ppid;
	} // GetParentProcessID
} /** namespace win32 **/
#endif


ProcInfo :: ProcInfo () {
	this->GetAttrs ();
} // ProcInfo :: ProcInfo


ProcInfo :: ~ProcInfo () {
	this->clear ();
} // ProcInfo :: ~ProcInfo


void ProcInfo :: clear () {
	this->proc_name.clear ();
}


void ProcInfo :: GetAttrs () {
	proc_id = GetProcId ();
	pproc_id = GetPProcId ();
	thread_id = GetThreadId ();
	GetProcName (proc_name);
} // ProcInfo :: GetAttrs


unsigned long ProcInfo :: GetProcId () 
{
	boost::process::self & s = self::get_instance ();
	boost::process::process::id_type id = s.get_id ();
	return (unsigned long) id;
} // ProcInfo :: GetProcId


unsigned long ProcInfo :: GetPProcId ()
{
#if defined WIN32 || defined _WIN32
	unsigned long pid = ProcInfo::GetProcId ();
	unsigned long id = win32::GetParentProcessID (pid);
#else
		pid_t id = getppid ();
#endif
		return (unsigned long) id;
	} // ProcInfo :: GetPProcId


const char * ProcInfo :: GetProcName (std::string & name) 
{
#ifdef WIN32 || _WIN32
	int namelen;
	win32::GetProcessorName (name, namelen);
#else
	char * _ProcName = (char *) calloc (namelen, sizeof (char));
	int ret_code = gethostname (_ProcName, namelen);
	if (ret_code == 0) 
		name.assign (_ProcName);
	if (_ProcName != NULL)  free (_ProcName);
#endif

		return name.c_str();
	} // ProcInfo :: GetProcName


const char * ProcInfo :: GetCPUName (std::string & name) 
{
	return ProcInfo::GetProcName (name);
} // ProcInfo :: GetCPUName


unsigned long ProcInfo :: GetThreadId (int base) 
{
	/**
	Get this_thread id from boost, then convert to a string (hex),
	last, convert it to an integer.
	**/
	boost::thread::id & myid = boost::this_thread::get_id ();
	std::stringstream idstr;
	idstr << myid;
	return std::stoul (idstr.str(), nullptr, base); 
} // ProcInfo :: GetThreadId


std::ostream& operator<< (std::ostream& os, ProcInfo& object)
{
	object.GetAttrs ();
	os << object.proc_name << " (pid " << object.proc_id
		<< ", ppid " << object.pproc_id
		<< ", thread " << object.thread_id
		<< ")";
	return os;
} // operator<<