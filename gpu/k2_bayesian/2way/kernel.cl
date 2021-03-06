
#define TABLE_MAX_SIZE 748 // 5060
#define TABLE_ERROR -0.0810 //-0.08104
#define MAX_K 2
#define MAX_R_SIZE 9 // = 3^MAX_K
#pragma OPENCL EXTENSION cl_khr_fp64 : enable


// Returns value from addlog table or approximation, depending on number
double addlog(int n, __constant double* addlogtable)
{
	if(n < TABLE_MAX_SIZE)
		return addlogtable[n];
	else
	{
		double x = (n + 0.5)*log(convert_double(n)) - (n - 1)*log(convert_double(exp(convert_double(1)))) + TABLE_ERROR;
		return x;
	}
}

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

// Computes Bayesian K2 score for a given set of k SNPs
double bayesian(int* set, __constant unsigned int* data, __constant unsigned int* states, __constant double* addlogtable)
{
	double score;
	int i, j, x, r_index;
	int i0, i1, array[MAX_K];
	uint4 snps[MAX_K][3], state, result, result0, result1, count0, count1;

	// create frequency vector r
	uint r[MAX_R_SIZE][2];
	for(i = 0; i < MAX_R_SIZE; i++)
	{
		r[i][0] = 0;
		r[i][1] = 0;
	}

	// loop data columns
	for(j = 0; j + 3 < NCOLS; j += 4)
	{
		for(x = 0; x < MAX_K; x++)
		{
			snps[x][0] = (uint4) (data[(set[x] * NCOLS + j + 0) * 3 + 0], data[(set[x] * NCOLS + j + 1) * 3 + 0], data[(set[x] * NCOLS + j + 2) * 3 + 0], data[(set[x] * NCOLS + j + 3) * 3 + 0]);
			snps[x][1] = (uint4) (data[(set[x] * NCOLS + j + 0) * 3 + 1], data[(set[x] * NCOLS + j + 1) * 3 + 1], data[(set[x] * NCOLS + j + 2) * 3 + 1], data[(set[x] * NCOLS + j + 3) * 3 + 1]);
			snps[x][2] = (uint4) (data[(set[x] * NCOLS + j + 0) * 3 + 2], data[(set[x] * NCOLS + j + 1) * 3 + 2], data[(set[x] * NCOLS + j + 2) * 3 + 2], data[(set[x] * NCOLS + j + 3) * 3 + 2]);
		}
		state = (uint4) (states[j + 0], states[j + 1], states[j + 2], states[j + 3]);

		// fill frequency table
		r_index = 0;
		for(i0 = 0; i0 <= 2; i0++)
		{
			for(i1 = 0; i1 <= 2; i1++)
			{
				// define array
				array[0] = i0;
				array[1] = i1;

				// fill frequency table
				result = snps[0][array[0]];
				for(i = 1; i < MAX_K; i++)
					result = result & snps[i][array[i]];
				result0 = result & ~state;
				result1 = result & state;
				count0 = popcount(result0);
				count1 = popcount(result1);
				r[r_index][0] += count0.x + count0.y + count0.z + count0.w;
				r[r_index][1] += count1.x + count1.y + count1.z + count1.w;
				r_index++;
			}
		}
		// end for loop
	}
	
	// compute the K2 score
	score = 0.0;
	for(i = 0; i < MAX_R_SIZE; i++)
		score += addlog(r[i][0] + r[i][1] + 1, addlogtable) - addlog(r[i][0], addlogtable) - addlog(r[i][1], addlogtable);
	score = fabs(convert_double(score));

	return score;
}


// Computes K2 score for a combination of SNPs
__kernel void kernel_bayesian(__global double* scores, __global int* solutions, __global int* combs, __global int* ncombs, __constant unsigned int* data, __constant unsigned int* states, __constant double* addlogtable)
{	
	// get thread ID and number of threads
	size_t t_id = get_global_id(0);
	size_t n_ts = get_global_size(0);

	// get combination
	int comb[MAX_K];
	comb[0] = combs[t_id * MAX_K + 0];
	comb[1] = combs[t_id * MAX_K + 1];

	// update combination index
	ncombs[t_id] += n_ts;

	// compute the K2 score
	double curr_score = bayesian(comb, data, states, addlogtable);
	if(curr_score < scores[t_id])
	{
		// update score and solution
		scores[t_id] = curr_score;
		solutions[t_id * MAX_K + 0] = comb[0];
		solutions[t_id * MAX_K + 1] = comb[1];
	}
	// end kernel	
}

// Computes next combination of SNPs
__kernel void kernel_combo(__global int* combs, __global int* ncombs)
{
	// get thread ID
	size_t t_id = get_global_id(0);

	// compute next combination
	int comb[MAX_K];
	if(ncombs[t_id] < NCOMB)
		getcombination(comb, SIZE, MAX_K, ncombs[t_id]);
	else
	{
		comb[0] = 0;
		comb[1] = 1;
	}

	// update combination
	combs[t_id * MAX_K + 0] = comb[0];
	combs[t_id * MAX_K + 1] = comb[1];
	// end kernel
}