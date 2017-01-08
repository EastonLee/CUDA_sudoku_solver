/*
 ============================================================================
 Name        : c_sudoku.c
 Author      : Easton Lee
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA_sudoku_solver
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#define ROWS 25
#define COLS ROWS
#define ROOT 5
#define NO_SOLUTION 0
#define MAY_HAVE_SOLUTION 1
#define GOT_SOLUTION 2
#define NO_CHANGE_SO_PAUSE -1
//typedef (long long) EASTON_TYPE;

#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }

#define Runtime(...) do{\
	struct timeval start, end;\
	gettimeofday(&start, NULL);\
	__VA_ARGS__;\
	gettimeofday(&end, NULL);\
	printf(": Used time %f\n",(double) (end.tv_usec - start.tv_usec) / 1000000 + (double)(end.tv_sec-start.tv_sec));\
}while(0)

__device__ bool *d_B_change_occur , *d_B_no_solution , *d_B_got_solution ;
 bool h_B_change_occur = 1, h_B_no_solution = 0, h_B_got_solution = 0;

long long original_matrix[ROWS][COLS] = { 0 };
//{ROOT,0,7,9,0,0,0,2,0,4,0,0,0,0,0,0,0,0,0,0,0,3,6,0,1,0,5,0,8,0,0,0,0,2,1,7,0,0,0,4,0,0,0,0,0,0,0,0,0,0,9,0,5,0,0,0,8,0,0,3,0,0,9,0,0,0,0,7,2,0,0,0,0,0,3,0,0,0,6,0,0};
long long matrix[ROWS][COLS] ={0};
//{0,0,0,0,0,0,0,8,0,0,9,0,0,0,0,0,0,0,7,0,8,9,0,0,0,4,0,0,0,4,0,0,0,0,7,0,0,0,0,0,6,0,2,0,0,8,0,0,0,0,2,3,9,0,2,0,0,3,0,0,0,0,0,3,0,0,5,0,0,6,0,7,0,4,0,6,0,1,5,0,0};
//{ 0, 0, 4, 9, 5, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3,	0, 0, 1, 0, 0, 0, 0, 9, 0, 8, 1, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	0, 5, 7, 0, 0, 6, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 8,	0, 5, 0, 0, 1, 0, 9, 5, 0, 9, 0, 0, 0, 4, 8, 0 };
//{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//{0,0,4,9,5,0,0,0,0,0,0,0,2,0,0,0,3,0,0,1,0,0,0,0,9,0,8,1,9,0,0,0,0,0,0,0,0,0,0,0,0,0,5,7,0,0,6,3,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,8,0,5,0,0,1,0,9,5,0,9,0,0,0,4,8,0};
//{3,0,7,9,0,0,0,2,0,4,0,0,0,0,0,0,0,0,0,0,0,3,6,0,1,0,5,0,8,0,0,0,0,2,1,7,0,0,0,4,0,0,0,0,0,0,0,0,0,0,9,0,5,0,0,0,8,0,0,3,0,0,9,0,0,0,0,7,2,0,0,0,0,0,3,0,0,0,6,0,0};
__device__ long long (*d_matrix)[ROWS][COLS];//donnot read this shit-->//donot treat as a pointer, but as a matrix

int tree_length = 0, divide_times = 0, tree_change_times=0;
__device__ int *d_conflict_pos , *d_exhaust_pos;
int h_conflict_pos = -1, h_exhaust_pos = -1;


int simple_node_link[ROWS * COLS * 2]  ;

typedef struct node {
	long long (*current_matrix)[ROWS][COLS];
	int divide_from_x, divide_from_y;
	long long current_candicate;
	struct node * p_prev_node;
} node;

//node * head_node;
//node * tree;

void print_simple_node_link(int simple_node_link[ROWS * COLS * 2])
{
	int i;
	printf("\nnow the simple node link is like this\n");
	for(i=0;simple_node_link[i]!=-1;i++)
	{
		if(i%2 == 0)
			printf("(%d ",simple_node_link[i]);
		else
			printf("%d) ",simple_node_link[i]);
	}
	printf("\n\n");
}
__host__ __device__ int bit_count(long long foo) {
	long long i, sum = 0, tmp;
	for (i = 0; i < ROWS; i++) {
		tmp = 1 << i;
		sum += ((foo & tmp) >> i);
	}
	return sum;
}

__host__ __device__ int highest_bit(long long foo) {
	int i, tmp;
	for (i = ROWS-1; i >= 0; i--) {
		tmp = 1 << i;
		if ((foo & tmp) >> i)
			return i + 1;
	}
	return 0;
}

int calc_least_candicate(long long (*foo)[ROWS][COLS], int *least_x, int *least_y,
		int *least_candicate) {
	int i, tmp_x = -1, tmp_y = -1, tmp_least = 0xFF;
	for (i = 0; i < ROWS * COLS; i++) {
		if (((*foo)[i / ROWS][i % COLS]) & 0x8000000000000000) {
			long long tmp = (*foo)[i / ROWS][i % COLS];
			int bit_num = bit_count(tmp);
			if (bit_num < tmp_least) {
				tmp_x = i / COLS;
				tmp_y = i % COLS;
				tmp_least = bit_count((*foo)[tmp_x][tmp_y]);
			}
		}

	}

	*least_x = tmp_x;
	*least_y = tmp_y;
	*least_candicate = tmp_least;
	return tmp_least;
}

void print_matrix(long long foo[ROWS][COLS]) {
	int i;

	printf("\nThe current matrix is like this\n");
	for (i = 0; i < ROWS * COLS; i++) {
		if (!(i % COLS)) {
			printf("\n");
		}
		if ((foo[i / COLS][i % COLS]) & 0x8000000000000000)
			printf("_ ");
		else
			printf("%d ", foo[i / COLS][i % COLS]);

	}
	printf("\n");
}

void print_candicate_num(long long foo[ROWS][COLS]) {
	int i;

	printf("\nThe number of potential solutions is like this");
	for (i = 0; i < ROWS * COLS; i++) {
		if (!(i % COLS)) {
			printf("\n");
		}
		if ((foo[i / COLS][i % COLS]) & 0x8000000000000000)
			printf("%d ", bit_count(foo[i / COLS][i % COLS]));
		else
			printf("_ ");

	}
	printf("\n");
}

__global__ void VecAdd(const float* A, const float* B, float* C, int N) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < N)
		C[i] = A[i] + B[i];
}

__global__ void kernel_row_check(long long (*d_matrix)[ROWS][COLS],bool *d_B_change_occur,bool *d_B_no_solution,bool *d_B_got_solution,int *d_conflict_pos,int *d_exhaust_pos) {
	int i, bit_map = 0;

	*d_B_change_occur=	*d_B_no_solution=	*d_B_got_solution =0;
	int row = threadIdx.x;
	for (i = 0; i < COLS; i++) //calc bit_map
			{
		if (!(((*d_matrix)[row][i]) & 0x8000000000000000)) { //int debug = bit_map|(1<<((*d_matrix)[row][i]-1));int debug2 = 4|64;
			if ((bit_map | (1 << ((*d_matrix)[row][i] - 1))) == bit_map) {
				//printf("\nelement conflict, wrong branch\n");//TODO: i cannot use printf in kernel
				atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution |= 1;
				*d_conflict_pos = row * COLS + i;
				return;
			}
			bit_map |= 1 << ((*d_matrix)[row][i] - 1);
		}
	}

	for (i = 0; i < COLS; i++) {
		if (((*d_matrix)[row][i]) & 0x8000000000000000) {
			if (((~bit_map) & ((*d_matrix)[row][i])) != (*d_matrix)[row][i]) //some possibility should be cut
					{
				((*d_matrix)[row][i]) &= (~bit_map);
				//d_B_change_occur=1;
			}
			int tmp = bit_count((*d_matrix)[row][i]);
//			if (tmp == 0) {
//				*d_B_no_solution = 1;
//				d_exhaust_pos = row * COLS + i;
//				return;
//			}
			if (tmp == 1) //only one possible is left, we consider it's the right one
					{

				(*d_matrix)[row][i] = highest_bit((*d_matrix)[row][i]);
				atomicOr((int*)d_B_change_occur,1);//*d_B_change_occur |= 1;

				bit_map = 0; //update bit_map
				for (i = 0; i < COLS; i++) //recalc bit_map
						{
					if (!(((*d_matrix)[row][i]) & 0x8000000000000000)) {
						if ((bit_map | (1 << ((*d_matrix)[row][i] - 1)))
								== bit_map) {
							//printf("\nelement conflict, wrong branch\n");//TODO: i cannot use printf in kernel
							atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution |= 1;
							return;
						}
						bit_map |= 1 << ((*d_matrix)[row][i] - 1);
					}
				}

				//else return NO_SOLUTION;

			}

			if ((*d_matrix)[row][i] == 0x8000000000000000) //no possible value left, there is no solution, the puzzle is wrong
					{
				//printf("\nno value for matrix[%d][%d]\n", row, i);//TODO: i cannot use printf in kernel
				atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution |= 1;
				*d_exhaust_pos = row * COLS + i;
				return;
			}
		}
	}
	return;
}

__global__ void kernel_col_check(long long (*d_matrix)[ROWS][COLS],bool *d_B_change_occur,bool *d_B_no_solution,bool *d_B_got_solution,int *d_conflict_pos,int *d_exhaust_pos) {
	int i, bit_map = 0;

	int col = threadIdx.x;
	for (i = 0; i < ROWS; i++) //calc bit_map
			{
		if (!(((*d_matrix)[i][col]) & 0x8000000000000000)) {
			//int debug = (1<<((*d_matrix)[i][col]-1));int debug2 = bit_map|debug;
			if ((bit_map | (1 << ((*d_matrix)[i][col] - 1))) == bit_map) {
				//printf("\nelement conflict, wrong branch\n");//TODO
				atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
				*d_conflict_pos = i * COLS + col;
				return;
			}
			bit_map |= 1 << ((*d_matrix)[i][col] - 1);
		}
	}

	for (i = 0; i < ROWS; i++) {
		if (((*d_matrix)[i][col]) & 0x8000000000000000) {
			if (((~bit_map) & ((*d_matrix)[i][col])) != (*d_matrix)[i][col]) {
				((*d_matrix)[i][col]) &= (~bit_map);
				//d_B_change_occur=1;
			}

			int tmp = bit_count((*d_matrix)[i][col]);
			if (tmp == 0) {
				atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
				*d_exhaust_pos = i * COLS + col;
				return;
			}
			if (tmp == 1) //only one possible is left, we consider it's the right one
					{
				(*d_matrix)[i][col] = highest_bit((*d_matrix)[i][col]);
				atomicOr((int*)d_B_change_occur,1);//*d_B_change_occur = true;

				bit_map = 0; //update bit_map
				for (i = 0; i < ROWS; i++) //recalc bit_map
						{
					if (!(((*d_matrix)[i][col]) & 0x8000000000000000)) {
						if ((bit_map | (1 << ((*d_matrix)[i][col] - 1)))
								== bit_map) {
							//printf("\nelement conflict, wrong branch\n");//TODO
							atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
							return;
						}
						bit_map |= 1 << ((*d_matrix)[i][col] - 1);
					}
				}
			}

			if ((*d_matrix)[i][col] == 0x8000000000000000) {
				//printf("\nno value for matrix[%d][%d]\n", i, col);//TODO
				atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
				*d_exhaust_pos = i * COLS + col;
				return;
			}
		}
	}
	return;
}

__global__ void kernel_block_check(long long (*d_matrix)[ROWS][COLS],bool *d_B_change_occur,bool *d_B_no_solution,bool *d_B_got_solution,int *d_conflict_pos,int *d_exhaust_pos) {
	int i, j, bit_map = 0;

	int nth_block = threadIdx.x;
	int block_row = nth_block / ROOT;
	int block_cul = nth_block % ROOT;

	for (i = 0; i < ROOT; i++) //calc bit_map
			{
		for (j = 0; j < ROOT; j++) {
			if (!(((*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) & 0x8000000000000000)) {
				if ((bit_map
						| (1
								<< ((*d_matrix)[block_row * ROOT + i][block_cul * ROOT
										+ j] - 1))) == bit_map) {
					//printf("\nelement conflict, wrong branch\n");//TODO
					atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
					*d_conflict_pos = i * COLS + j;
					return;
				}
				bit_map |= 1
						<< ((*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]
								- 1);
			}

		}

	}

	for (i = 0; i < ROOT; i++) {
		for (j = 0; j < ROOT; j++) {

			if (((*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) & 0x8000000000000000) {
				if (((~bit_map)
						& ((*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]))
						!= (*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) {
					((*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) &=
							(~bit_map);
					//d_B_change_occur=1;
				}

				int tmp = bit_count(
						(*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]);
				if (tmp == 0) {
					atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
					*d_exhaust_pos = i * COLS + j;
					return;
				}
				if (tmp == 1) //only one possible is left, we consider it's the right one
						{
					(*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j] =
							highest_bit(
									(*d_matrix)[block_row * ROOT + i][block_cul * ROOT
											+ j]);
					atomicOr((int*)d_B_change_occur,1);//*d_B_change_occur = 1;


					bit_map = 0; //update bit_map
					for (i = 0; i < ROOT; i++) //recalc bit_map
							{
						for (j = 0; j < ROOT; j++) {
							if (!(((*d_matrix)[block_row * ROOT + i][block_cul * ROOT
									+ j]) & 0x8000000000000000)) {
								if ((bit_map
										| (1
												<< ((*d_matrix)[block_row * ROOT
														+ i][block_cul * ROOT + j]
														- 1))) == bit_map) {
									//printf("\nelement conflict, wrong branch\n");//TODO
									atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
									return;
								}
								bit_map |=
										1
												<< ((*d_matrix)[block_row * ROOT
														+ i][block_cul * ROOT + j]
														- 1);
							}

						}

					}
				}

				if ((*d_matrix)[block_row * ROOT + i][block_cul * ROOT + j]
						== 0x8000000000000000) {
					//printf("\nno value for matrix[%d][%d]\n", block_row * ROOT + i,block_cul * ROOT + j);//TODO
					atomicOr((int*)d_B_no_solution,1);//*d_B_no_solution = 1;
					*d_exhaust_pos = i * COLS + j;
					return;
				}
			}

		}

	}
	return;
}

void divide_to_new(node** foo) //donot know why, reference parameter cannot be used(like &parameter), so replace with its pointer
		{
	node *p_new_node = (node*) malloc(sizeof(struct node));
	p_new_node->current_matrix = (long long (*)[ROWS][COLS]) malloc(
			ROWS * COLS * sizeof(long long));
	memcpy(p_new_node->current_matrix, (*foo)->current_matrix,
			ROWS * COLS * sizeof(long long));
	(*(p_new_node->current_matrix))[(*foo)->divide_from_x][(*foo)->divide_from_y] =
			(*foo)->current_candicate;

	p_new_node->divide_from_x = p_new_node->divide_from_y =
			p_new_node->current_candicate = 0;
	p_new_node->p_prev_node = *foo;
	simple_node_link[tree_length * 2] = (*foo)->divide_from_x * COLS
			+ (*foo)->divide_from_y;
	simple_node_link[tree_length * 2 + 1] = (*foo)->current_candicate;
	*foo = p_new_node;

	tree_length++;
	tree_change_times++;
	divide_times++;
	printf("\ntree changes %d times\n",tree_change_times);

	printf("\ntry a new branch&add a new node\nit has %d nodes now!!!\n",
			tree_length);

}

void choose_best_candicate(node ** foo) {
	int least_candicate, least_x, least_y;

	calc_least_candicate((*foo)->current_matrix, &least_x, &least_y,
			&least_candicate);

	if (least_candicate == 0xFF) {//branch exhaust
		h_B_no_solution = 1; //send a branch-over signal
		h_exhaust_pos = 0xFF;
		return;
	}
	(*foo)->divide_from_x = least_x;
	(*foo)->divide_from_y = least_y;
	(*foo)->current_candicate = highest_bit(
			(*((*foo)->current_matrix))[least_x][least_y]);
	divide_to_new(foo);

	print_simple_node_link(simple_node_link);

}

int search_solution(long long (*p_matrix)[ROWS][COLS]) {
	int i;

	h_B_change_occur = 1;
	h_B_no_solution = 0;
	int debug_count = 0;
	h_conflict_pos = h_exhaust_pos = -1;

	CUDA_CHECK_RETURN(cudaMemcpy(&(*d_matrix[0][0]),&((*p_matrix)[0][0]), sizeof(long long)*ROWS*COLS,cudaMemcpyHostToDevice));
//	CUDA_CHECK_RETURN(cudaMemcpyToSymbol(d_matrix, *p_matrix, sizeof(long long)*ROWS*COLS, cudaMemcpyHostToDevice));
	while (h_B_change_occur && !h_B_no_solution) //everytime we search a matrix, we must loop the three phase until no change happen
	{
		debug_count++;
		h_B_change_occur = 0; //*d_B_no_solution = 0;
		//for (i = 0; (i < ROWS) && !*d_B_no_solution; i++) {
		kernel_row_check<<<1, ROWS>>>(d_matrix,d_B_change_occur,d_B_no_solution,d_B_got_solution,d_conflict_pos,d_exhaust_pos);
		CUDA_CHECK_RETURN(cudaMemcpy(&h_B_change_occur,d_B_change_occur,sizeof(bool),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_B_no_solution,d_B_no_solution,sizeof(bool),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_conflict_pos,d_conflict_pos,sizeof(int),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_exhaust_pos,d_exhaust_pos,sizeof(int),cudaMemcpyDeviceToHost));
		//}
		if(h_B_no_solution) break;
		//for (i = 0; (i < COLS) && !d_B_no_solution; i++) {
		kernel_col_check<<<1,ROWS>>>(d_matrix,d_B_change_occur,d_B_no_solution,d_B_got_solution,d_conflict_pos,d_exhaust_pos);
		CUDA_CHECK_RETURN(cudaMemcpy(&h_B_change_occur,d_B_change_occur,sizeof(bool),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_B_no_solution,d_B_no_solution,sizeof(bool),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_conflict_pos,d_conflict_pos,sizeof(int),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_exhaust_pos,d_exhaust_pos,sizeof(int),cudaMemcpyDeviceToHost));
		//}
		if(h_B_no_solution) break;
		//for (i = 0; (i < ROWS) && !d_B_no_solution; i++) {
		kernel_block_check<<<1,ROWS>>>(d_matrix,d_B_change_occur,d_B_no_solution,d_B_got_solution,d_conflict_pos,d_exhaust_pos);
		CUDA_CHECK_RETURN(cudaMemcpy(&h_B_change_occur,d_B_change_occur,sizeof(bool),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_B_no_solution,d_B_no_solution,sizeof(bool),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_conflict_pos,d_conflict_pos,sizeof(int),cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(&h_exhaust_pos,d_exhaust_pos,sizeof(int),cudaMemcpyDeviceToHost));
		//}
		if(h_B_no_solution) break;
	}

	if (h_B_no_solution) {

		return NO_SOLUTION;
	}

	else {
		CUDA_CHECK_RETURN(cudaMemcpy(*p_matrix, *d_matrix, sizeof(long long)*ROWS*COLS, cudaMemcpyDeviceToHost));//only this result deserve copied out
		int tmp = calc_least_candicate(p_matrix, &i, &i, &i);

		if (tmp == 0xFF) {
			h_B_got_solution = 1;
			return GOT_SOLUTION;
		}
		if (tmp > 1)
			return NO_CHANGE_SO_PAUSE;
	}

}

void backward_on_tree(node **p_tree) {
	int biggest_candicate;
	free((*p_tree)->current_matrix);
	free((*p_tree));
	(*p_tree) = (*p_tree)->p_prev_node; //backward
	tree_length--;
	tree_change_times++;
	printf("\ntree changes %d times\n",tree_change_times);
	simple_node_link[tree_length * 2] = -1;
	simple_node_link[tree_length * 2 + 1] = -1;
	printf("\ndelete a wrong node,\nit has %d nodes now!!!\n", tree_length);
	print_simple_node_link(simple_node_link);

	long long (*p_matrix)[ROWS][COLS] = (*p_tree)->current_matrix;
	biggest_candicate = highest_bit(
			(*p_matrix)[(*p_tree)->divide_from_x][(*p_tree)->divide_from_y]); //cut the highest bit
	(*p_matrix)[(*p_tree)->divide_from_x][(*p_tree)->divide_from_y] &= (~(1
			<< (biggest_candicate - 1)));

	int debug_bit_count = bit_count(
			(*p_matrix)[(*p_tree)->divide_from_x][(*p_tree)->divide_from_y]);
	if (bit_count(
			(*p_matrix)[(*p_tree)->divide_from_x][(*p_tree)->divide_from_y])
			== 1) {
		(*p_matrix)[(*p_tree)->divide_from_x][(*p_tree)->divide_from_y] =
				highest_bit(
						(*p_matrix)[(*p_tree)->divide_from_x][(*p_tree)->divide_from_y]);
		search_solution((*p_tree)->current_matrix);

	}

}

int main(void) {
	int i, tmp;
	tree_length = 0;
	tree_change_times  =0;
	memset(simple_node_link,-1,ROWS * COLS * 2*sizeof(int));

	int *debug_tree_length = &tree_length;
	int *debug_tree_change_times =&tree_change_times;
	int *debug_divide_times = &divide_times;
	int *debug_conflict_pos = &h_conflict_pos;
	int *debug_exhaust_pos = &h_exhaust_pos;
	int (*debug_simple_node_link)[ROWS * COLS * 2] = &simple_node_link;

	for (i = 0; i < ROWS * COLS; i++) {
		//scanf("%d",&tmp);
		tmp = matrix[i / COLS][i % COLS];
		matrix[i / COLS][i % COLS] = tmp ? tmp : 0xFFFFFFFFFFFFFFFF; //first bit means this is a solution, last 9 bits is the solution map;
		original_matrix[i / COLS][i % COLS] = tmp ? tmp : 0xFFFFFFFFFFFFFFFF;
	}



	//above is the c version
	//---------------------------------
	//below is the cuda version

	print_matrix(matrix);

	node head_node;
	head_node.current_matrix = &matrix;
	head_node.current_candicate = head_node.divide_from_x =
			head_node.divide_from_y = 0;

	node * tree;
	tree = &head_node;

	//choose_best_candicate(tree);
	//d_matrix = (long long(*)[ROWS][COLS])malloc(sizeof(long long*));//it is a global matrix initiated at the start, no need to alloc
	CUDA_CHECK_RETURN(cudaMalloc((void **)&d_matrix, (size_t)sizeof(long long)*ROWS*COLS));//TODO:remember to free the address
	CUDA_CHECK_RETURN(cudaMalloc((void **)&d_B_change_occur, (size_t)sizeof(bool)));
	CUDA_CHECK_RETURN(cudaMalloc((void **)&d_B_no_solution, (size_t)sizeof(bool)));
	CUDA_CHECK_RETURN(cudaMalloc((void **)&d_B_got_solution, (size_t)sizeof(bool)));
	CUDA_CHECK_RETURN(cudaMalloc((void **)&d_conflict_pos, (size_t)sizeof(int)));
	CUDA_CHECK_RETURN(cudaMalloc((void **)&d_exhaust_pos, (size_t)sizeof(int)));

	Runtime(
		int result = search_solution(tree->current_matrix); //result only comes from search_solution
		printf("\nafter search, ");
		print_matrix(*(tree->current_matrix));
		while ((result == NO_CHANGE_SO_PAUSE) || (result == NO_SOLUTION)) {
			if (result == NO_SOLUTION) {
				if (tree == &head_node) {
					break; //game is over, totally no solution!
				}
				backward_on_tree(&tree); //call choose_best_candicate, so tree may be changed again!!!!!!
			}

			if (h_B_got_solution == 1) {
				break;
			}

			if (h_B_no_solution == 1) //branch is over, no need to search
					{
				result = NO_SOLUTION;
				continue;
			}

			choose_best_candicate(&tree); //!!!NOTE:tree has been changed to a new one!

			if (h_B_no_solution == 1) //branch is over, no need to search
					{
				result = NO_SOLUTION;
				continue;
			}

			result = search_solution(tree->current_matrix); //result only comes from search_solution
			printf("\nafter search, ");
			print_matrix(*(tree->current_matrix));
		}

	);
	if (h_B_got_solution == 1) {
		printf("\nGot a solution!\n");
		print_matrix(*(tree->current_matrix));

	}

	else {
		printf("No solution!\nThe last found matrix is like this\n");
		print_matrix(*(tree->current_matrix));
		print_candicate_num(*(tree->current_matrix));

	}

	//if(*d_matrix)//XXX
		//cudaFree(*d_matrix);
	if(d_matrix)//no need to free
		cudaFree(d_matrix);

	return EXIT_SUCCESS;
}
