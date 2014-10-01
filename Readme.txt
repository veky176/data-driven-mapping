	In this folder, you can find the source codes and
input data-dependence matricies.

To run the Strassen algorithm, you need:

[1] You need to install Boost and Boost.Process and ThreadPool.
	[1.a] Boost 1.55.0.
	[1.b] Boost.Process (http://www.highscore.de/boost/process/)
	[1.c] ThreadPool 0.2.5 (http://threadpool.sourceforge.net/index.html)

[2] You can find the source codes for Strassen in ./strassen.
	The main funcion is in ./strassen/strassen_main.cpp.
	The data-dependency matrix files are:
	[2.a] Strassen method  : ./strassen-dmp/mat_strassen.dep
	[2.b] Native MM method : ./strassen-dmp/mat_multi.dep

[3] You can use MSVC (Microsoft Visual Studio 2010) to run it.
	Command Arguments: 5000 mat_strassen.dep 10
	It means matrix size = 5000, using 10 processes, and
	using the data-dependence matrix in mat_strassen.dep.
	
	
To run the Cholesky inversion algorithm, you need:

[1] Do the same procedure as mention in above [1].

[2] You need the source codes in ./strassen and ./cholesky.
	The main funcion is in ./cholesky/cholesky_main.cpp.
	The data-dependency matrix file is:
	./cholesky-dmp/mat_cholesky.dep
	
[3] You need to generate a SPD matrix as input for Cholesky.
	Here is a program that can help create a SPD matrix.
	You can find the code: ./cholesky-dmp/generateSPDmatrix.m
	Usage:
	[3.a] A = generateSPDmatrix(5000); % create the matrix
	[3.b] dlmwrite ('input_5000.txt', A, ' '); % save the matrix
	Then, make sure you have the file - input_5000.txt in your
	working directory.	
	
[4] You can use MSVC to run the program.
	Command Arguments: 5000 mat_cholesky.dep 7
	It means matrix size = 5000, using 7 processes, and
	using the data-dependence matrix in mat_cholesky.dep.
	
	6/5/2014