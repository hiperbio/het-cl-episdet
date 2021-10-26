
// Libraries
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <sstream>
#include <omp.h>
#include <immintrin.h>
using namespace std;

// Macros
#define TABLE_MAX_SIZE 748
#define TABLE_ERROR -0.0810

// Classes
class SNP
{
	public:

	int samplesize;			// #patients (rows in the file data)
	int locisize; 			// #SNPs
	int data_col;			// #SNPs + class (columns in the file data)
	unsigned long ***data;	// SNPs data matrix
	unsigned long *states;	// Disease states vector
	int nrows;				// #rows in data matrix
	int ncols;				// #columns in data matrix
	char** SNPnames;		// names of the SNPs (first row of file)
	int classvalues[2];		// #patients in class '0' and class '1'
	
	void input_data(const char* path);	// input data reader
    void generate_data(int num_pac, int num_snp);
	void destroy();						// destructor
};

// Destructor
void SNP::destroy()
{
	int i, j;

	// delete matrix "SNPnames"
	for(i = 0; i < locisize; i++)
		delete []SNPnames[i];
	delete []SNPnames;
	
	// delete matrix "data"
	for(i = 0; i < nrows - 1; i++)
	{
		for(j = 0; j < 3; j++)
			delete []data[i][j];
		delete []data[i];
	}
	delete []data;

	// delete vector "states"
	delete []states;
}

// Input data reader
void SNP::generate_data(int num_pac, int num_snp)
{
	int i, j, x, new_j, temp;
	string line;
	unsigned char **tdata;

	data_col = num_snp + 1;
	locisize = num_snp;
	samplesize = num_pac;
	nrows = data_col;
	ncols = samplesize;

	// create matrixes "SNPnames" and "tdata"
	SNPnames = new char*[data_col];
	for(i = 0; i < data_col; i++)
		SNPnames[i] = new char[50];
	tdata = new unsigned char*[nrows];
	for(i = 0; i < nrows; i++)
		tdata[i] = new unsigned char[ncols];
	
	// initialize "classvalues"
	classvalues[0] = 0;
	classvalues[1] = 0;

	// fill matrix "tdata"
    srand(100);
    for(j = 0; j < ncols; j++)
        for(i = 0; i < nrows - 1; i++)
            tdata[i][j] = rand() % 3;
    i = data_col - 1;
    for(j = 0; j < ncols; j++)
    {
        temp = rand() % 2;
        tdata[i][j] = temp;
        classvalues[temp]++;
    }

	// convert data matrix
	ncols = ceil(1.0 * ncols / 64);
	data = new unsigned long**[nrows - 1];
	for(i = 0; i < nrows - 1; i++)
	{
		data[i] = new unsigned long*[3];
		for(j = 0; j < 3; j++)
			data[i][j] = new unsigned long[ncols];
	}
	states = new unsigned long[ncols];

	for(i = 0; i < nrows - 1; i++)
	{
		new_j = 0;
		for(j = 0; j + 63 < samplesize; j += 64)
		{
			// set to zero
			data[i][0][new_j] = 0;
			data[i][1][new_j] = 0;
			data[i][2][new_j] = 0;
			// loop through the 64 columns
			for(x = 0; x < 64; x++)
			{
				// left shift by 1
				data[i][0][new_j] <<= 1;
				data[i][1][new_j] <<= 1;
				data[i][2][new_j] <<= 1;
				// set appropriate position to 1
				data[i][tdata[i][j + x]][new_j] |= 1;
			}
			// get disease state
			states[new_j] = 0;
			for(x = 0; x < 64; x++)
			{
				states[new_j] <<= 1;
				states[new_j] |= tdata[nrows - 1][j + x];
			}
			// update index
			new_j++;
		}

		// repeat for remainder
		if(j != samplesize)
		{
			data[i][0][new_j] = 0;
			data[i][1][new_j] = 0;
			data[i][2][new_j] = 0;
			for(x = 0; j + x < samplesize; x++)
			{
				data[i][0][new_j] <<= 1;
				data[i][1][new_j] <<= 1;
				data[i][2][new_j] <<= 1;
				data[i][tdata[i][j + x]][new_j] |= 1;
			}
			states[new_j] = 0;
			for(x = 0; j + x < samplesize; x++)
			{
				states[new_j] <<= 1;
				states[new_j] |= tdata[nrows - 1][j + x];
			}
			new_j++;
		}
	}
	
	// delete matrix "tdata"
	for(i = 0; i < nrows; i++)
		delete []tdata[i];
	delete []tdata;

	// end
}

// Input data reader
void SNP::input_data(const char* path)
{
	int i, j, x, new_j, temp;
	string line;
	unsigned char **tdata;

	ifstream in(path);
	getline(in, line);
	istringstream test(line);
	string word;

	// get "data_col" and "locisize"
	i = 0;
	while(!test.eof())
	{
		getline(test, word, ',');
		i++;
	}
	data_col = i;
	locisize = i - 1;

	// get "samplesize"
	j = 0;
	while(!in.eof())
	{   
		getline(in, line);
		j++;
	}
	samplesize = j - 1;
	in.close();
	
	nrows = data_col;
	ncols = samplesize;

	// create matrixes "SNPnames" and "tdata"
	SNPnames = new char*[data_col];
	for(i = 0; i < data_col; i++)
		SNPnames[i] = new char[50];
	tdata = new unsigned char*[nrows];
	for(i = 0; i < nrows; i++)
		tdata[i] = new unsigned char[ncols];
	
	// initialize "classvalues"
	classvalues[0] = 0;
	classvalues[1] = 0;

	ifstream in1(path);
	getline(in1, line);
	istringstream test1(line);

	// fill matrix "SNPnames"
	i = 0;
	while(!test1.eof())
	{
		if(i == data_col)
			break;
		getline(test1, word, ',');
		strcpy(SNPnames[i], word.c_str());
		i++;
	}

	// fill matrix "tdata"
	i = 0;
	while(!in1.eof())
	{
		if(i == samplesize)
			break;
		getline(in1, line);
		istringstream values(line);
		j = 0;
		while(!values.eof())
		{
			getline(values, word, ',');
			istringstream int_iss(word);
			int_iss >> temp;
			tdata[j++][i] = temp;
			if(j == data_col)
				classvalues[temp]++;
		}
		i++;
	}
	in1.close();

	// convert data matrix
	ncols = ceil(1.0 * ncols / 64);
	data = new unsigned long**[nrows - 1];
	for(i = 0; i < nrows - 1; i++)
	{
		data[i] = new unsigned long*[3];
		for(j = 0; j < 3; j++)
			data[i][j] = new unsigned long[ncols];
	}
	states = new unsigned long[ncols];

	for(i = 0; i < nrows - 1; i++)
	{
		new_j = 0;
		for(j = 0; j + 63 < samplesize; j += 64)
		{
			// set to zero
			data[i][0][new_j] = 0;
			data[i][1][new_j] = 0;
			data[i][2][new_j] = 0;
			// loop through the 64 columns
			for(x = 0; x < 64; x++)
			{
				// left shift by 1
				data[i][0][new_j] <<= 1;
				data[i][1][new_j] <<= 1;
				data[i][2][new_j] <<= 1;
				// set appropriate position to 1
				data[i][tdata[i][j + x]][new_j] |= 1;
			}
			// get disease state
			states[new_j] = 0;
			for(x = 0; x < 64; x++)
			{
				states[new_j] <<= 1;
				states[new_j] |= tdata[nrows - 1][j + x];
			}
			// update index
			new_j++;
		}

		// repeat for remainder
		if(j != samplesize)
		{
			data[i][0][new_j] = 0;
			data[i][1][new_j] = 0;
			data[i][2][new_j] = 0;
			for(x = 0; j + x < samplesize; x++)
			{
				data[i][0][new_j] <<= 1;
				data[i][1][new_j] <<= 1;
				data[i][2][new_j] <<= 1;
				data[i][tdata[i][j + x]][new_j] |= 1;
			}
			states[new_j] = 0;
			for(x = 0; j + x < samplesize; x++)
			{
				states[new_j] <<= 1;
				states[new_j] |= tdata[nrows - 1][j + x];
			}
			new_j++;
		}
	}
	
	// delete matrix "tdata"
	for(i = 0; i < nrows; i++)
		delete []tdata[i];
	delete []tdata;

	// end
}

// Global variables
SNP SNPs;
double *addlogtable;
double time_begin, time_end;

// Functions

// Returns value from addlog table or approximation, depending on number
double addlog(int n)
{
	if(n < TABLE_MAX_SIZE)
		return addlogtable[n];
	else
	{
		double x = (n + 0.5)*log(n) - (n - 1)*log(exp(1)) + TABLE_ERROR;
		return x;
	}
}

// Computes addlog table
double my_factorial(int n)
{
	double i;
	double z = 0;

	if(n < 0)
	{
		printf("Error: n should be a non-negative number.\n");
		return 0;
	}
	if(n == 0)
		return 0;
	if(n == 1)
		return 0;

	z = addlogtable[n - 1] + log(n);
	return z;
}

// Fills frequency table for a given
void fill_table(int *array, int prev_i, int pos, int k, uint r[][2], int *r_index, __m256i snps[][3], __m256i state)
{
    int curr_i, i;
    __m256i result, result0, result1;
	unsigned long aux0[4], aux1[4];

    if(prev_i <= 2)
    {
        for(curr_i = 0; curr_i <= 2; curr_i++)
        {
            array[pos] = curr_i;
            if(pos + 1 < k)
            {
                // recall function
                fill_table(array, curr_i, pos + 1, k, r, r_index, snps, state);
            }
            else
            {
				// fill frequency table
                result = snps[0][array[0]];
                for(i = 1; i < k; i++)
                    result = _mm256_and_si256(result, snps[i][array[i]]);
                result0 = _mm256_andnot_si256(state, result);
                result1 = _mm256_and_si256(state, result);
				_mm256_storeu_si256((__m256i *) aux0, result0);
				_mm256_storeu_si256((__m256i *) aux1, result1);
				for(i = 0; i < 4; i++)
				{
					r[*r_index][0] += _popcnt64(aux0[i]);
					r[*r_index][1] += _popcnt64(aux1[i]);
				}
                (*r_index)++;
            }
        }
    }
    // end
}

// Computes Bayesian K2 score for a given set of k SNPs
double bayesian(int r_size, int *set, int k)
{
	double score;
	int i, j, x, r_index, array[k];

	int m = SNPs.nrows;
	int n = SNPs.ncols;
	int old_n = SNPs.samplesize;
	__m256i snps[k][3], state, mask_v;
	unsigned long mask[4];
	
	// create frequency vector r
	uint r[r_size][2];
	for(i = 0; i < r_size; i++)
	{
		r[i][0] = 0;
		r[i][1] = 0;
	}

	// loop data columns
	for(j = 0; j + 3 < n; j += 4)
	{
		for(x = 0; x < k; x++)
		{
			snps[x][0] = _mm256_loadu_si256((__m256i const *) &(SNPs.data[set[x]][0][j]));
			snps[x][1] = _mm256_loadu_si256((__m256i const *) &(SNPs.data[set[x]][1][j]));
			snps[x][2] = _mm256_loadu_si256((__m256i const *) &(SNPs.data[set[x]][2][j]));
		}
		state = _mm256_loadu_si256((__m256i const *) &(SNPs.states[j]));
		
		// fill frequency table
		r_index = 0;
		for(i = 0; i <= 2; i++)
		{
			array[0] = i;
			fill_table(array, i, 1, k, r, &r_index, snps, state);
		}
	}

	// repeat for remainder
	if(j < n)
	{
		// set mask
		i = 0;
		for(x = j; x < n; x++)
		{
			mask[i] = 0x8000000000000000;
			i++;
		}
		while(i < 4)
		{
			mask[i] = 0;
			i++;
		}
		mask_v = _mm256_loadu_si256((__m256i const *) &(mask[0]));

		for(x = 0; x < k; x++)
		{
			snps[x][0] = _mm256_maskload_epi64((long const *) &(SNPs.data[set[x]][0][j]), mask_v);
			snps[x][1] = _mm256_maskload_epi64((long const *) &(SNPs.data[set[x]][1][j]), mask_v);
			snps[x][2] = _mm256_maskload_epi64((long const *) &(SNPs.data[set[x]][2][j]), mask_v);
		}
		state = _mm256_maskload_epi64((long const *) &(SNPs.states[j]), mask_v);

		// fill frequency table
		r_index = 0;
		for(i = 0; i <= 2; i++)
		{
			array[0] = i;
			fill_table(array, i, 1, k, r, &r_index, snps, state);
		}
	}

	// compute the K2 score
	score = 0.0;
	for(i = 0; i < r_size; i++)
		score += addlog(r[i][0] + r[i][1] + 1) - addlog(r[i][0]) - addlog(r[i][1]);
	score = fabs(score);

	return score;
}

// Recursive function to go through all combinations of size choose k
void recursive(int* comb, int* best_sol, double* best_score, int prev_i, int position, int size, int k, int* t_count, int r_size)
{
	int curr_i, t_id, n_ts, x;
	double curr_score;

	// get thread information and table size
	t_id = omp_get_thread_num();
	n_ts = omp_get_num_threads();

	// run if previous element wasn't the final on the data set
	if(prev_i + 1 < size)
	{
		// iterate through every remaining value for this position
		for(curr_i = prev_i + 1; curr_i < size; curr_i++)
		{
			comb[position] = curr_i;
			// advance to next position if we haven't reached the final combination
			if(position + 1 < k)
			{
				recursive(comb, best_sol, best_score, curr_i, position + 1, size, k, t_count, r_size);
			}
			// if we have a final combination but it's not to be processed by this thread
			else if(*t_count != 0)
			{
				// decrease counter and ignore combination
				(*t_count)--;
			}
			// if we have a final combination and it's adequate for this thread
			else
			{
				// get score
				curr_score = bayesian(r_size, comb, k);
				// compare score
				if(curr_score < *best_score)
				{
					// update score and solution
					*best_score = curr_score;
					for(x = 0; x < k; x++)
						best_sol[x] = comb[x];
				}
				// reset counter
				*t_count += n_ts - 1;
			}
		}
	}
}

// Exhaustively searches for the best Bayesian K2 score using the recursive method
double exhaustive(int k, int* sol)
{
	double score, scores[omp_get_max_threads()];
	int solutions[omp_get_max_threads()][k], x, j, r_size;

	r_size = (int) pow(3.0, k);

	time_begin = omp_get_wtime();
	#pragma omp parallel
	{
		int i, t_id, comb[k], best_sol[k], position, size, t_count;
		double best_score;

		// get thread information
		t_id = omp_get_thread_num();

		// set parameters for first recursive function call
		best_score = 100000000.0;
		position = 1;
		size = SNPs.locisize;
		t_count = t_id;

		// call recursive function
		for(i = 0; i <= size - k; i++)
		{
			comb[0] = i;
			recursive(comb, best_sol, &best_score, i, position, size, k, &t_count, r_size);
		}

		// register score and solution
		scores[t_id] = best_score;
		for(i = 0; i < k; i++)
			solutions[t_id][i] = best_sol[i];
		
		// synchronize threads
		#pragma omp barrier
		// end of parallel region
	}

	// get best score and solution from all threads
	score = 100000000.0;
	for(x = 0; x < omp_get_max_threads(); x++)
	{
		if(scores[x] < score)
		{
			score = scores[x];
			for(j = 0; j < k; j++)
				sol[j] = solutions[x][j];
		}
	}
	time_end = omp_get_wtime();

	return score;
}

// Main
int main(int argc, char **argv)
{
	int i, k, addlogsize;
	double score, begin, end;

	if(argc == 4)
    {
        int num_pac = atoi(argv[1]);
        int num_snp = atoi(argv[2]);
        k = atoi(argv[3]);
        SNPs.generate_data(num_pac, num_snp);
    }
    else if(argc == 3)
    {
        k = atoi(argv[2]);
        SNPs.input_data(argv[1]);
    }
    int sol[k];

	// create addlog table (up to TABLE_MAX_SIZE positions at max)
	if(SNPs.samplesize <= TABLE_MAX_SIZE)
		addlogsize = SNPs.samplesize;
	else
		addlogsize = TABLE_MAX_SIZE;
	addlogtable = new double[addlogsize];
	for(i = 0; i < addlogsize; i++)
		addlogtable[i] = my_factorial(i);

	// call exhaustive search function
	cout << "Starting computation of K2 Score..." << endl;
	score = exhaustive(k, sol);
	cout << "... done!" << endl << "K2 Score: " << score << endl;
	cout << "Solution: ";
	for(i = 0; i < k; i++)
		cout << sol[i] << " ";
	cout << endl;
	double interval = double(time_end - time_begin);

	// free dynamic memory
	SNPs.destroy();
	delete addlogtable;

	// display elapsed time
	cout << "Total time: " << interval << " s" << endl;
	return 0;
}
