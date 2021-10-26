
#define TABLE_MAX_SIZE 748
#define TABLE_ERROR -0.0810
#define TWOLOG_MAX_SIZE 1500
#define MAX_K 3
#define MAX_R_SIZE 27 // = 3^MAX_K
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

// Returns value from twolog table or log2() result, depending on number
double twolog(int n, __constant double* twologtable)
{
	if(n < TWOLOG_MAX_SIZE)
		return twologtable[n];
	else
	{
		double x = log2(convert_double(n));
		return x;
	}
}

// Computes Bayesian K2 score for a given set of k SNPs
void bayesian(int* set, double *scores, __constant unsigned int* data, __constant unsigned int* states, __constant double* addlogtable, __constant double* twologtable)
{
	double k2_score, me_score, hx, hxy;
	int i, j, x, r_index;
	int i0, i1, i2, array[MAX_K];
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
				for(i2 = 0; i2 <= 2; i2++)
				{
					// define array
					array[0] = i0;
					array[1] = i1;
					array[2] = i2;

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
		}
		// end for loop
	}

	// compute the K2 and ME scores
	k2_score = 0.0;
	hx = 0.0;
	hxy = 0.0;
	for(i = 0; i < MAX_R_SIZE; i++)
	{
		k2_score += addlog(r[i][0] + r[i][1] + 1, addlogtable) - addlog(r[i][0], addlogtable) - addlog(r[i][1], addlogtable);
		hx = hx - (r[i][0] + r[i][1]) * (twolog(r[i][0] + r[i][1], twologtable) - LOG2TOTAL_ME);
		hxy = hxy - r[i][0] * (twolog(r[i][0], twologtable) - LOG2TOTAL_ME);
		hxy = hxy - r[i][1] * (twolog(r[i][1], twologtable) - LOG2TOTAL_ME);
	}
	k2_score = fabs(convert_double(k2_score));
	me_score = HY_ME + (hx - hxy)/TOTAL_ME;
	me_score = 1/me_score;

	// end
	scores[0] = k2_score;
	scores[1] = me_score;
}


// Kernel
__kernel void kernel_bayesian(__global double* scores, __global int* combinations, __constant unsigned int* data, __constant unsigned int* states, __constant double* addlogtable, __constant double* twologtable)
{
	// get current work item id
	int t_id = get_global_id(0);

	// create variables to use in score calculation
	int comb[MAX_K];
	double curr_scores[2];

	// get combination
	comb[0] = combinations[t_id * MAX_K + 0];
	comb[1] = combinations[t_id * MAX_K + 1];
	comb[2] = combinations[t_id * MAX_K + 2];

	// compute new scores
	bayesian(comb, curr_scores, data, states, addlogtable, twologtable);
	
	// register scores
	scores[t_id * 2 + 0] = curr_scores[0];
	scores[t_id * 2 + 1] = curr_scores[1];

	// end of kernel
}