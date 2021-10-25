
// Libraries
#include <CL/opencl.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <immintrin.h>

// Macros
#define TABLE_MAX_SIZE 748
#define TABLE_ERROR -0.0810
#define MAX_K 2
#define MAX_R_SIZE 9 // = 3^MAX_K
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define MY_DEVICE_TYPE CL_DEVICE_TYPE_GPU
#define N_WORK_ITEMS 76800
#define N_WORK_ITEMS_PER_GROUP 256
#define INITIAL_CHUNK_SIZE 8000


// SNP class
class SNP
{
	public:

	int samplesize;			// #patients (rows in the file data)
	int locisize; 			// #SNPs
	int data_col;			// #SNPs + class (columns in the file data)
	unsigned int ***data;	// SNPs data matrix
	unsigned int *states;	// Disease states vector
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
	std::string line;
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
	ncols = ceil(1.0 * ncols / 32);
	while(ncols % 8 != 0)
		ncols++;
	data = new unsigned int**[nrows - 1];
	for(i = 0; i < nrows - 1; i++)
	{
		data[i] = new unsigned int*[3];
		for(j = 0; j < 3; j++)
			data[i][j] = new unsigned int[ncols];
	}
	states = new unsigned int[ncols];

	for(i = 0; i < nrows - 1; i++)
	{
		new_j = 0;
		for(j = 0; j + 31 < samplesize; j += 32)
		{
			// set to zero
			data[i][0][new_j] = 0;
			data[i][1][new_j] = 0;
			data[i][2][new_j] = 0;
			// loop through the 32 columns
			for(x = 0; x < 32; x++)
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
			for(x = 0; x < 32; x++)
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
	std::string line;
	unsigned char **tdata;

	std::ifstream in(path);
	getline(in, line);
	std::istringstream test(line);
	std::string word;

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

	std::ifstream in1(path);
	getline(in1, line);
	std::istringstream test1(line);

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
		std::istringstream values(line);
		j = 0;
		while(!values.eof())
		{
			getline(values, word, ',');
			std::istringstream int_iss(word);
			int_iss >> temp;
			tdata[j++][i] = temp;
			if(j == data_col)
				classvalues[temp]++;
		}
		i++;
	}
	in1.close();

	// convert data matrix
	ncols = ceil(1.0 * ncols / 32);
	while(ncols % 8 != 0)
		ncols++;
	data = new unsigned int**[nrows - 1];
	for(i = 0; i < nrows - 1; i++)
	{
		data[i] = new unsigned int*[3];
		for(j = 0; j < 3; j++)
			data[i][j] = new unsigned int[ncols];
	}
	states = new unsigned int[ncols];

	for(i = 0; i < nrows - 1; i++)
	{
		new_j = 0;
		for(j = 0; j + 31 < samplesize; j += 32)
		{
			// set to zero
			data[i][0][new_j] = 0;
			data[i][1][new_j] = 0;
			data[i][2][new_j] = 0;
			// loop through the 32 columns
			for(x = 0; x < 32; x++)
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
			for(x = 0; x < 32; x++)
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

// Computes n choose r
unsigned int choose(unsigned int n, unsigned int k)
{
    if(k > n)
		return 0;
    if(k * 2 > n)
		k = n - k;
    if(k == 0)
		return 1;

    unsigned int result = n;
    for(int i = 2; i <= k; i++)
	{
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

// Computes the mth combination in the set of n choose r combinations
void getcombination(int* comb, int n, int r, int m)
{
	int a, b, x, i;

	a = n;
	b = r;
	x = choose(n, r) - 1 - m;

	for(i = 0; i < r; i++)
	{
		a--;
		while(choose(a, b) > x)
			a--;
		comb[i] = n - 1 - a;
		x = x - choose(a, b);
		b--;
	}
}

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
double bayesian(int *set, int k)
{
	double score;
	int i, j, x, r_index, array[k];

	int m = SNPs.nrows;
	int n = SNPs.ncols;
	int old_n = SNPs.samplesize;
	__m256i snps[k][3], state, mask_v;
	unsigned int mask[8];
	
	// create frequency vector r
	uint r[MAX_R_SIZE][2];
	for(i = 0; i < MAX_R_SIZE; i++)
	{
		r[i][0] = 0;
		r[i][1] = 0;
	}

	// loop data columns
	for(j = 0; j + 7 < n; j += 8)
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
	// compute the K2 score
	score = 0.0;
	for(i = 0; i < MAX_R_SIZE; i++)
		score += addlog(r[i][0] + r[i][1] + 1) - addlog(r[i][0]) - addlog(r[i][1]);
	score = fabs(score);

	return score;
}


// Main
int main(int argc, char** argv)
{
    int i, j, k, addlogsize, offset;
	
	k = 2;
    if(argc == 3)
    {
        int num_pac = atoi(argv[1]);
        int num_snp = atoi(argv[2]);
        SNPs.generate_data(num_pac, num_snp);
    }
    else if(argc == 2)
    {
        SNPs.input_data(argv[1]);
    }

    // get OpenCL platform and device
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    if(all_platforms.size() == 0)
    {
        std::cout << "No platforms found." << std::endl;
        exit(1);
    }
    cl::Platform default_platform;
    for(i = 0; i < all_platforms.size(); i++)
    {
        if(std::string(&(all_platforms[i].getInfo<CL_PLATFORM_NAME>())[0]).find("Graphics") != std::string::npos)
            default_platform = all_platforms[i];
    }
    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
    if(all_devices.size() == 0)
    {
        std::cout << "No devices found." << std::endl;
        exit(1);
    }
    cl::Device default_device = all_devices[0];
    std::cout << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << std::endl;

    // create OpenCL context
    cl::Context context(default_device);

    // create OpenCL program
    std::ifstream program_file("kernel.cl");
    std::string program_string(std::istreambuf_iterator<char>(program_file), (std::istreambuf_iterator<char>()));
    cl::Program::Sources program_source;
    program_source.push_back({program_string.c_str(), program_string.length() + 1});
    cl::Program program(context, program_source);
    
	char build_options[100];
	sprintf(build_options, "-D NCOLS=%d", SNPs.ncols);
    if(program.build(build_options) != CL_SUCCESS)
    {
        std::cout << "Error building:" << std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    // create OpenCL command queue
    cl::CommandQueue queue(context, default_device, CL_QUEUE_PROFILING_ENABLE);

	// create addlog table (up to TABLE_MAX_SIZE positions at max)
	addlogsize = TABLE_MAX_SIZE;
	addlogtable = new double[addlogsize];
	for(i = 0; i < addlogsize; i++)
		addlogtable[i] = my_factorial(i);

    // create OpenCL buffers on device
    cl::Buffer dev_scores(context, CL_MEM_READ_WRITE, N_WORK_ITEMS * sizeof(double));
    cl::Buffer dev_solutions(context, CL_MEM_READ_WRITE, N_WORK_ITEMS * MAX_K * sizeof(int));
	cl::Buffer dev_combinations(context, CL_MEM_READ_WRITE, N_WORK_ITEMS * MAX_K * sizeof(int));
    cl::Buffer dev_data(context, CL_MEM_READ_ONLY, (SNPs.nrows - 1) * SNPs.ncols * 3 * sizeof(unsigned int));
    cl::Buffer dev_states(context, CL_MEM_READ_ONLY, SNPs.ncols * sizeof(unsigned int));
    cl::Buffer dev_addlogtable(context, CL_MEM_READ_ONLY, addlogsize * sizeof(double));

    // copy buffers from host to device
	offset = 0;
	for(i = 0; i < SNPs.nrows - 1; i++)
	{
		for(j = 0; j < 3; j++)
		{
			queue.enqueueWriteBuffer(dev_data, CL_TRUE, offset * sizeof(unsigned int), SNPs.ncols * sizeof(unsigned int), SNPs.data[i][j]);
			offset += SNPs.ncols;
		}
	}
	queue.enqueueWriteBuffer(dev_states, CL_TRUE, 0, SNPs.ncols * sizeof(unsigned int), SNPs.states);
    queue.enqueueWriteBuffer(dev_addlogtable, CL_TRUE, 0, addlogsize * sizeof(double), addlogtable);

	// set buffer containing scores in device
	queue.enqueueFillBuffer(dev_scores, 100000000.0, 0, N_WORK_ITEMS * sizeof(double));

	// setup OpenCL kernel call
    cl::Kernel kernel_bayesian(program, "kernel_bayesian");
    kernel_bayesian.setArg(0, dev_scores);
    kernel_bayesian.setArg(1, dev_solutions);
	kernel_bayesian.setArg(2, dev_combinations);
    kernel_bayesian.setArg(3, dev_data);
    kernel_bayesian.setArg(4, dev_states);
    kernel_bayesian.setArg(5, dev_addlogtable);

	// create OpenCL event to monitor execution status
	cl::Event execution_event;

	// create buffers to receive results
    int n_ts = omp_get_max_threads();
    double score, scores_gpu[N_WORK_ITEMS], scores_cpu[n_ts];
	int sol[k], sols_gpu[N_WORK_ITEMS][MAX_K], sols_cpu[n_ts][MAX_K];
	for(i = 0; i < N_WORK_ITEMS; i++)
		scores_gpu[i] = 100000000.0;
	for(i = 0; i < n_ts; i++)
		scores_cpu[i] = 100000000.0;
	score = 100000000.0;

	// create buffer with information to transfer
	int combinations_gpu[N_WORK_ITEMS][k];
	int total_combs = choose(SNPs.locisize, k);
	// initialize event
	queue.enqueueWriteBuffer(dev_combinations, CL_TRUE, 0, 1 * MAX_K * sizeof(int), combinations_gpu, nullptr, &execution_event);

	// create array for thread synchronization
	volatile bool sync[n_ts + 1];
	for(i = 0; i < n_ts; i++)
		sync[i] = false;
	sync[n_ts] = true;
	volatile int chunk = INITIAL_CHUNK_SIZE;
	volatile int combinations_cpu;
	volatile double p_cpu;

	// start counting time
	time_begin = omp_get_wtime();

	// spawn OpenMP threads
	#pragma omp parallel
	{
		// create thread-private variables
		int comb[k], n, a, b, c, i0, i1;
		int t_id;
		double best_score, curr_score;
		int best_sol[k], initial_comb;
		int counter = 0;

		best_score = 100000000.0;
		// get thread ID
		t_id = omp_get_thread_num();

		// Thread 0 works as a scheduler
		if(t_id == 0)
		{
			int e, total;
			bool d;
			double t_gpu;
			cl_ulong t_gpu_start;
			cl_ulong t_gpu_end;
			// start counting time
			std::cout << "Starting computation of K2 Score..." << std::endl;

			if(total_combs > N_WORK_ITEMS)
			{
				// fill buffer with combinations
				total = 0;
				a = 0;
				for(i0 = 0; i0 <= SNPs.locisize - k; i0++)
				{
					for(i1 = i0 + 1; i1 < SNPs.locisize; i1++)
					{
						combinations_gpu[a][0] = i0;
						combinations_gpu[a][1] = i1;
						a++;
						// check if buffer is filled
						if(a == N_WORK_ITEMS)
						{
							total += a;
							while(a != 0)
							{
								// check if GPU is free
								if(execution_event.getInfo<CL_EVENT_COMMAND_EXECUTION_STATUS>() == 0)
								{
									// get GPU performance
									t_gpu_start = execution_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
									t_gpu_end = execution_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
									t_gpu = double(t_gpu_end - t_gpu_start) / 1E9; // OCL times are in nanoseconds
									// send buffer to device
									queue.enqueueWriteBuffer(dev_combinations, CL_TRUE, 0, N_WORK_ITEMS * MAX_K * sizeof(int), combinations_gpu);
									// enqueue OpenCL buffer
									queue.enqueueNDRangeKernel(kernel_bayesian, cl::NullRange, cl::NDRange(N_WORK_ITEMS), cl::NDRange(N_WORK_ITEMS_PER_GROUP), nullptr, &execution_event);
									// reset GPU buffer iterator
									a = 0;
								}

								// check if CPU is free
								d = false;
								for(b = 1; b < n_ts; b++)
									d = d | sync[b];
								if(d == false)
								{
									// get initial combination for CPU
									combinations_cpu = total;
									
									// get new CPU chunk size
									chunk = t_gpu * p_cpu;
									if(chunk <= 0)
										chunk = INITIAL_CHUNK_SIZE;

									// signal other threads to start processing
									for(b = 1; b < n_ts; b++)
										sync[b] = true;
									
									// get last combination to be processed
									if(total + (n_ts - 1) * chunk < total_combs)
										total = total + (n_ts - 1) * chunk;
									else
										total = total_combs;
								}

								if(chunk > INITIAL_CHUNK_SIZE)
								{
									counter++;
									getcombination(comb, SNPs.locisize, k, total);
									curr_score = bayesian(comb, k);
									// update best score and solution
									if(curr_score < best_score)
									{
										best_score = curr_score;
										for(e = 0; e < k; e++)
											best_sol[e] = comb[e];
									}
									total++;
								}
							}
							// update current last combination
							getcombination(comb, SNPs.locisize, k, total);
							i0 = comb[0];
							i1 = comb[1];
						}
					}
				}
				// signal other threads to end
				sync[n_ts] = false;
				// register scores and solutions
				scores_cpu[t_id] = best_score;
				for(e = 0; e < k; e++)
					sols_cpu[t_id][e] = best_sol[e];
				// send buffer to device
				queue.enqueueWriteBuffer(dev_combinations, CL_TRUE, 0, a * MAX_K * sizeof(int), combinations_gpu);
				// enqueue OpenCL kernel
				queue.enqueueNDRangeKernel(kernel_bayesian, cl::NullRange, cl::NDRange(a), cl::NDRange(1));
				// get scores and solutions
				queue.enqueueReadBuffer(dev_scores, CL_TRUE, 0, N_WORK_ITEMS * sizeof(double), scores_gpu);
				queue.enqueueReadBuffer(dev_solutions, CL_TRUE, 0, N_WORK_ITEMS * MAX_K * sizeof(int), (*sols_gpu));
			}
			else
			{
				// if the total number of combinations is low, let the CPU cores do all the work
				combinations_cpu = 0;
				chunk = total_combs;
				// signal other threads to start processing
				for(b = 1; b < n_ts; b++)
					sync[b] = true;
				// signal other threads to end
				sync[n_ts] = false;
				// register scores and solutions
				scores_cpu[t_id] = best_score;
				for(a = 0; a < k; a++)
					sols_cpu[t_id][a] = best_sol[a];
			}
		}
		else // Remaining threads work on score calculation
		{
			double t_cpu;
			int mult;
			// infinite loop
			while(1)
			{
				// check if thread has pending work
				if(sync[t_id] == true)
				{
					t_cpu = omp_get_wtime();
					// set up parameters
					initial_comb = combinations_cpu;
					n = chunk * (n_ts - 1) + 1;
					getcombination(comb, SNPs.locisize, k, initial_comb);
					i0 = comb[0];
					i1 = comb[1];
					mult = t_id;
					// process combinations
					for(; i0 <= SNPs.locisize - k && n != 0; i0++)
					{
						for(; i1 < SNPs.locisize && n != 0; i1++)
						{
							n--;
							mult--;
							if(mult == 0)
							{
								counter++;
								// get combination
								comb[0] = i0;
								comb[1] = i1;
								// compute score
								curr_score = bayesian(comb, k);
								// update best score and solution
								if(curr_score < best_score)
								{
									best_score = curr_score;
									for(c = 0; c < k; c++)
										best_sol[c] = comb[c];
								}
								mult = n_ts;
							}
						}
						i1 = i0 + 2;
					}
					// get scheduler performance
					p_cpu = chunk / (omp_get_wtime() - t_cpu);
					// signal thread as free
					sync[t_id] = false;
				}
				// check for termination signal
				if(sync[n_ts] == false && sync[t_id] == false)
				{
					// register scores and solutions
					scores_cpu[t_id] = best_score;
					for(a = 0; a < k; a++)
						sols_cpu[t_id][a] = best_sol[a];
					// leave infinite loop
					break;
				}
			}
		}
		// end of parallel region
	}
    // compute final score and solution
    for(i = 0; i < N_WORK_ITEMS; i++)
    {
        if(scores_gpu[i] < score)
        {
            score = scores_gpu[i];
            for(j = 0; j < k; j++)
                sol[j] = sols_gpu[i][j];
        }
    }
	for(i = 0; i < n_ts; i++)
	{
		if(scores_cpu[i] < score)
		{
			score = scores_cpu[i];
			for(j = 0; j < k; j++)
				sol[j] = sols_cpu[i][j];
		}
	}

	time_end = omp_get_wtime();
    std::cout << "... done!" << std::endl << "K2 Score: " << score << std::endl;
	std::cout << "Solution: ";
	for(i = 0; i < k; i++)
		std::cout << sol[i] << " ";
	std::cout << std::endl;
	double interval = double(time_end - time_begin);

	// free dynamic memory
	SNPs.destroy();
	delete addlogtable;

	// display elapsed time
	std::cout << "Total time: " << interval << " s" << std::endl;

	// end
	return 0;
}
