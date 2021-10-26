
// Libraries
#include <CL/opencl.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <omp.h>

// Macros
#define TWOLOG_MAX_SIZE 1500
#define MAX_K 2
#define MAX_R_SIZE 9 // = 3^MAX_K
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define MY_DEVICE_TYPE CL_DEVICE_TYPE_GPU
#define N_WORK_ITEMS 76800
#define N_WORK_ITEMS_PER_GROUP 256


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

// SNP destructor
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
		for(j = 0; j < ncols; j++)
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
	while(ncols % 4 != 0)
		ncols++;
	data = new unsigned int**[nrows - 1];
	for(i = 0; i < nrows - 1; i++)
	{
		data[i] = new unsigned int*[ncols];
		for(j = 0; j < ncols; j++)
			data[i][j] = new unsigned int[3]();
	}
	states = new unsigned int[ncols]();

	for(i = 0; i < nrows - 1; i++)
	{
		new_j = 0;
		for(j = 0; j + 31 < samplesize; j += 32)
		{
			// set to zero
			data[i][new_j][0] = 0;
			data[i][new_j][1] = 0;
			data[i][new_j][2] = 0;
			// loop through the 64 columns
			for(x = 0; x < 32; x++)
			{
				// left shift by 1
				data[i][new_j][0] <<= 1;
				data[i][new_j][1] <<= 1;
				data[i][new_j][2] <<= 1;
				// set appropriate position to 1
				data[i][new_j][tdata[i][j + x]] |= 1;
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
			data[i][new_j][0] = 0;
			data[i][new_j][1] = 0;
			data[i][new_j][2] = 0;
			for(x = 0; j + x < samplesize; x++)
			{
				data[i][new_j][0] <<= 1;
				data[i][new_j][1] <<= 1;
				data[i][new_j][2] <<= 1;
				data[i][new_j][tdata[i][j + x]] |= 1;
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
	while(ncols % 4 != 0)
		ncols++;
	data = new unsigned int**[nrows - 1];
	for(i = 0; i < nrows - 1; i++)
	{
		data[i] = new unsigned int*[ncols];
		for(j = 0; j < ncols; j++)
			data[i][j] = new unsigned int[3]();
	}
	states = new unsigned int[ncols]();

	for(i = 0; i < nrows - 1; i++)
	{
		new_j = 0;
		for(j = 0; j + 31 < samplesize; j += 32)
		{
			// set to zero
			data[i][new_j][0] = 0;
			data[i][new_j][1] = 0;
			data[i][new_j][2] = 0;
			// loop through the 64 columns
			for(x = 0; x < 32; x++)
			{
				// left shift by 1
				data[i][new_j][0] <<= 1;
				data[i][new_j][1] <<= 1;
				data[i][new_j][2] <<= 1;
				// set appropriate position to 1
				data[i][new_j][tdata[i][j + x]] |= 1;
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
			data[i][new_j][0] = 0;
			data[i][new_j][1] = 0;
			data[i][new_j][2] = 0;
			for(x = 0; j + x < samplesize; x++)
			{
				data[i][new_j][0] <<= 1;
				data[i][new_j][1] <<= 1;
				data[i][new_j][2] <<= 1;
				data[i][new_j][tdata[i][j + x]] |= 1;
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
double *twologtable, hy_me, log2total_me;
int total_me;
double time_begin, time_end;

// Functions

// Returns value from twolog table or log2() result, depending on number
double twolog(int n)
{
	if(n < TWOLOG_MAX_SIZE)
		return twologtable[n];
	else
	{
		double x = log2(n);
		return x;
	}
}

// Main
int main(int argc, char** argv)
{
    int i, i0, i1, j, k, offset;

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
    int sol[k];

	// compute global variables
	total_me = SNPs.samplesize;
	log2total_me = log2(total_me);
	hy_me = - (1.0*SNPs.classvalues[0]/total_me) * log2((1.0*SNPs.classvalues[0]/total_me)) - (1.0*SNPs.classvalues[1]/total_me) * log2((1.0*SNPs.classvalues[1]/total_me));

	// create twolog table (up to TWOLOG_MAX_SIZE positions at max)
	int twologsize = TWOLOG_MAX_SIZE;
	twologtable = new double[twologsize];
	for(i = 1; i < twologsize; i++)
		twologtable[i] = log2(i);
    twologtable[0] = 0.0f;

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
	sprintf(build_options, "-D NCOLS=%d -D HY_ME=%f -D LOG2TOTAL_ME=%f -D TOTAL_ME=%d", SNPs.ncols, hy_me, log2total_me, total_me);
    if(program.build(build_options) != CL_SUCCESS)
    {
        std::cout << "Error building:" << std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    // create OpenCL command queue
    cl::CommandQueue queue(context, default_device);

    // create OpenCL buffers on device
    cl::Buffer dev_scores(context, CL_MEM_READ_WRITE, N_WORK_ITEMS * sizeof(double));
    cl::Buffer dev_solutions(context, CL_MEM_READ_WRITE, N_WORK_ITEMS * MAX_K * sizeof(int));
	cl::Buffer dev_combinations(context, CL_MEM_READ_WRITE, N_WORK_ITEMS * MAX_K * sizeof(int));
    cl::Buffer dev_data(context, CL_MEM_READ_ONLY, (SNPs.nrows - 1) * SNPs.ncols * 3 * sizeof(unsigned int));
    cl::Buffer dev_states(context, CL_MEM_READ_ONLY, SNPs.ncols * sizeof(unsigned int));
    cl::Buffer dev_addlogtable(context, CL_MEM_READ_ONLY, twologsize * sizeof(double));

    // create buffer with information to transfer
	int combinations[N_WORK_ITEMS][k];

    // copy buffers from host to device
	offset = 0;
	for(i = 0; i < SNPs.nrows - 1; i++)
	{
		for(j = 0; j < SNPs.ncols; j++)
		{
			queue.enqueueWriteBuffer(dev_data, CL_TRUE, offset * sizeof(unsigned int), 3 * sizeof(unsigned int), SNPs.data[i][j]);
			offset += 3;
		}
	}
	queue.enqueueWriteBuffer(dev_states, CL_TRUE, 0, SNPs.ncols * sizeof(unsigned int), SNPs.states);
    queue.enqueueWriteBuffer(dev_addlogtable, CL_TRUE, 0, twologsize * sizeof(double), twologtable);

	// set buffer containing scores in device
	queue.enqueueFillBuffer(dev_scores, 100000000.0, 0, N_WORK_ITEMS * sizeof(double));

	// setup OpenCL kernel call
    std::cout << "Starting computation of K2 Score..." << std::endl;
    cl::Kernel kernel_bayesian(program, "kernel_bayesian");
    kernel_bayesian.setArg(0, dev_scores);
    kernel_bayesian.setArg(1, dev_solutions);
	kernel_bayesian.setArg(2, dev_combinations);
    kernel_bayesian.setArg(3, dev_data);
    kernel_bayesian.setArg(4, dev_states);
    kernel_bayesian.setArg(5, dev_addlogtable);

	// create buffers to receive results
    double score, scores[N_WORK_ITEMS];
    int sols[N_WORK_ITEMS][MAX_K];
	for(i = 0; i < N_WORK_ITEMS; i++)
		scores[i] = 100000000.0;

	// start counting time
	time_begin = omp_get_wtime();

	// compute combinations
	i = 0;
    for(i0 = 0; i0 <= SNPs.locisize - k; i0++)
	{
		for(i1 = i0 + 1; i1 < SNPs.locisize; i1++)
		{
			combinations[i][0] = i0;
			combinations[i][1] = i1;
			i++;
			if(i == N_WORK_ITEMS)
			{
				// make sure queue is clear
				queue.finish();
				// update combinations to be processed
				queue.enqueueWriteBuffer(dev_combinations, CL_TRUE, 0, N_WORK_ITEMS * MAX_K * sizeof(int), combinations);
				// run OpenCL kernel
				queue.enqueueNDRangeKernel(kernel_bayesian, cl::NullRange, cl::NDRange(N_WORK_ITEMS), cl::NDRange(N_WORK_ITEMS_PER_GROUP));
				i = 0;
			}
		}
	}
	// run once again for remainder
	queue.finish();
	queue.enqueueWriteBuffer(dev_combinations, CL_TRUE, 0, i * MAX_K * sizeof(int), combinations);
	queue.enqueueNDRangeKernel(kernel_bayesian, cl::NullRange, cl::NDRange(i), cl::NDRange(1));

    // copy results from device to host
	queue.finish();
    queue.enqueueReadBuffer(dev_scores, CL_TRUE, 0, N_WORK_ITEMS * sizeof(double), scores);
    queue.enqueueReadBuffer(dev_solutions, CL_TRUE, 0, N_WORK_ITEMS * MAX_K * sizeof(int), (*sols));
	
    // compute final score and solution
	score = 100000000.0;
    for(i = 0; i < N_WORK_ITEMS; i++)
    {
        if(scores[i] < score)
        {
            score = scores[i];
            for(j = 0; j < k; j++)
                sol[j] = sols[i][j];
        }
    }
	
	time_end = omp_get_wtime();
    std::cout << "... done!" << std::endl << "ME Score: " << score << std::endl;
	std::cout << "Solution: ";
	for(i = 0; i < k; i++)
		std::cout << sol[i] << " ";
	std::cout << std::endl;
	double interval = double(time_end - time_begin);
	SNPs.destroy();
	delete twologtable;

	// display elapsed time
	std::cout << "Total time: " << interval << " s" << std::endl;

	// end
	return 0;
}
