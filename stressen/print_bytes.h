#ifndef _PRINT_BYTES_H
#define _PRINT_BYTES_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <ctime>

using namespace std;

class PrintBytes {
	size_t _base;
	size_t _bytes;

public:
	PrintBytes (size_t bytes, size_t base = 1000) :
	  _bytes (bytes), _base (base) {}
	virtual ~PrintBytes () {}

	friend ostream& operator<< (ostream& os, PrintBytes & obj)
	{
		if (obj._bytes < obj._base) {
			os << obj._bytes << " bytes"; return os;
		}
		if (obj._bytes < obj._base * obj._base) {
			double base = pow ((double) obj._base, 1.0);
			os << std::fixed << std::setprecision (3) << (obj._bytes / base) << " KB"; return os;
		}
		if (obj._bytes < obj._base * obj._base * obj._base) {
			double base = pow ((double) obj._base, 2.0);
			os << std::fixed << std::setprecision (3) << (obj._bytes / base) << " MB"; return os;
		}
		double base = pow ((double) obj._base, 3.0);
		os << std::fixed << std::setprecision (3) << (obj._bytes / base) << " GB";
		return os;
	}

}; /** class PrintBtyes **/

template <typename DataType, typename SizeType>
class PrintMatrix {
	DataType * _ptr;
	SizeType _size;

public:
	PrintMatrix (DataType * ptr, SizeType size)
		: _ptr (ptr), _size (size) {
	}

	friend std::ostream & operator<< (std::ostream & os, PrintMatrix & obj)
	{
		os << "(ptr = " << obj._ptr << ") = ";
		for (SizeType id = 0; id < obj._size; id ++)
			os << " " << std::fixed << std::setprecision (2) << obj._ptr [id];
		os << endl;
		os.flush ();
		return os;
	}

}; /** PrintMatrix **/


template <typename DataType, typename SizeType>
class SaveMatrix {
private:
	DataType * ptr;
	SizeType row;
	SizeType col;
	std::string filename;

public:

	SaveMatrix (DataType * pointer, SizeType mat_size, const char * fname) 
		: ptr (pointer), row (mat_size), col (mat_size), filename (fname) {}

	SaveMatrix (DataType * pointer, SizeType mat_row, SizeType mat_col, const char * fname) 
		: ptr (pointer), row (mat_row), col (mat_col), filename (fname) {}

	void save (int decimal = 10) 
	{
		SimpleTimer timer (__FUNCTION__);

		std::ofstream file;
		file.open (filename.c_str());
		if (! file.is_open ()) return;
		for (SizeType r = 0; r < row; r++) {
			for (SizeType c = 0; c < col; c++)
				file << std::fixed << std::setprecision (decimal) << ptr [r * col + c] << " ";
			file << endl;
		}
		file.close ();
	}

}; /** class SaveMatrix **/

#endif