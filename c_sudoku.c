/*
 ============================================================================
 Name        : c_sudoku.c
 Author      : Easton Lee
 Version     :
 Copyright   : Your copyright notice
 Description : C_sudoku_solver
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#define ROWS 9
#define COLS ROWS
#define ROOT 3
#define NO_SOLUTION 0
#define MAY_HAVE_SOLUTION 1
#define GOT_SOLUTION 2
#define NO_CHANGE_SO_PAUSE -1
#define Runtime(...) do{\
	struct timeval start, end;\
	gettimeofday(&start, NULL);\
	__VA_ARGS__;\
	gettimeofday(&end, NULL);\
	printf("\n: Used time %f\n",(double) (end.tv_usec - start.tv_usec) / 1000000 + (double)(end.tv_sec-start.tv_sec));\
}while(0)

_Bool B_change_occur = 1, B_no_solution = 0, B_got_solution = 0;

unsigned long original_matrix[ROWS][COLS] = { 0 };
//{3,0,7,9,0,0,0,2,0,4,0,0,0,0,0,0,0,0,0,0,0,3,6,0,1,0,5,0,8,0,0,0,0,2,1,7,0,0,0,4,0,0,0,0,0,0,0,0,0,0,9,0,5,0,0,0,8,0,0,3,0,0,9,0,0,0,0,7,2,0,0,0,0,0,3,0,0,0,6,0,0};
unsigned long matrix[ROWS][COLS] ={0,0,0,0,0,0,0,8,0,0,9,0,0,0,0,0,0,0,7,0,8,9,0,0,0,4,0,0,0,4,0,0,0,0,7,0,0,0,0,0,6,0,2,0,0,8,0,0,0,0,2,3,9,0,2,0,0,3,0,0,0,0,0,3,0,0,5,0,0,6,0,7,0,4,0,6,0,1,5,0,0};
//{3,0,7,9,0,0,0,2,0,4,0,0,0,0,0,0,0,0,0,0,0,3,6,0,1,0,5,0,8,0,0,0,0,2,1,7,0,0,0,4,0,0,0,0,0,0,0,0,0,0,9,0,5,0,0,0,8,0,0,3,0,0,9,0,0,0,0,7,2,0,0,0,0,0,3,0,0,0,6,0,0};
//{0,0,0,0,0,0,0,8,0,0,9,0,0,0,0,0,0,0,7,0,8,9,0,0,0,4,0,0,0,4,0,0,0,0,7,0,0,0,0,0,6,0,2,0,0,8,0,0,0,0,2,3,9,0,2,0,0,3,0,0,0,0,0,3,0,0,5,0,0,6,0,7,0,4,0,6,0,1,5,0,0};
//{0,0,4,9,5,0,0,0,0,0,0,0,2,0,0,0,3,0,0,1,0,0,0,0,9,0,8,1,9,0,0,0,0,0,0,0,0,0,0,0,0,0,5,7,0,0,6,3,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,8,0,5,0,0,1,0,9,5,0,9,0,0,0,4,8,0};
//{3,0,7,9,0,0,0,2,0,4,0,0,0,0,0,0,0,0,0,0,0,3,6,0,1,0,5,0,8,0,0,0,0,2,1,7,0,0,0,4,0,0,0,0,0,0,0,0,0,0,9,0,5,0,0,0,8,0,0,3,0,0,9,0,0,0,0,7,2,0,0,0,0,0,3,0,0,0,6,0,0};

int tree_length = 0, divide_times = 0,tree_change_times=0;
int conflict_pos = -1, exhaust_pos = -1;
int simple_node_link[ROWS * COLS * 2] = { 0 };

typedef struct node {
	unsigned long (*current_matrix)[ROWS][COLS];
	int divide_from_x, divide_from_y;
	unsigned long current_candicate;
	struct node * p_prev_node;
} node;

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
//node * head_node;
//node * tree;

int bit_count(unsigned long foo) {
	unsigned long i, sum = 0, tmp;
	for (i = 0; i < ROWS; i++) {
		tmp = 1 << i;
		sum += ((foo & tmp) >> i);
	}
	return sum;
}

int highest_bit(unsigned long foo) {
	unsigned long i, tmp;
	for (i = ROWS-1; i >= 0; i--) {
		tmp = 1 << i;
		if ((foo & tmp) >> i)
			return i + 1;
	}
	return 0;
}

int calc_least_candicate(unsigned long (*foo)[ROWS][COLS], int *least_x, int *least_y,
		int *least_candicate) {
	int i, tmp_x = -1, tmp_y = -1, tmp_least = 0xFF;
	for (i = 0; i < ROWS*COLS; i++) {
		if (((*foo)[i / ROWS][i % ROWS]) & 0x8000000000000000) {
			unsigned long tmp = (*foo)[i / COLS][i % ROWS];
			int bit_num = bit_count(tmp);
			if (bit_num < tmp_least) {
				tmp_x = i / COLS;
				tmp_y = i % ROWS;
				tmp_least = bit_count((*foo)[tmp_x][tmp_y]);
			}
		}

	}

	*least_x = tmp_x;
	*least_y = tmp_y;
	*least_candicate = tmp_least;
	return tmp_least;
}

void print_matrix(unsigned long foo[ROWS][COLS]) {
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

void print_solution(unsigned long foo[ROWS][COLS]) {
	int i;

	printf("\nThe potential solution is like this");
	for (i = 0; i < ROWS * COLS; i++) {
		if (!(i % COLS)) {
			printf("\n");
		}
		if ((foo[i / COLS][i % COLS]) & 0x8000000000000000)
			printf("%d ", highest_bit(foo[i / COLS][i % COLS]));
		else
			printf("_ ");

	}
	printf("\n");
}
void row_check(unsigned long (*p_matrix)[ROWS][COLS], int row) {
	int i, bit_map = 0;

	for (i = 0; i < COLS; i++) //calc bit_map
			{
		if (!(((*p_matrix)[row][i]) & 0x8000000000000000)) { //int debug = bit_map|(1<<((*p_matrix)[row][i]-1));int debug2 = 4|64;
			if ((bit_map | (1 << ((*p_matrix)[row][i] - 1))) == bit_map) {
				printf("\nelement conflict, wrong branch\n");
				B_no_solution = 1;
				conflict_pos = row * ROWS + i;
				return;
			}
			bit_map |= 1 << ((*p_matrix)[row][i] - 1);
		}
	}

	for (i = 0; i < COLS; i++) {
		if (((*p_matrix)[row][i]) & 0x8000000000000000) {
			if (((~bit_map) & ((*p_matrix)[row][i])) != (*p_matrix)[row][i]) //some possibility should be cut
					{
				((*p_matrix)[row][i]) &= (~bit_map);
				//B_change_occur=1;
			}
			int tmp = bit_count((*p_matrix)[row][i]);
			if (tmp == 0) {
				B_no_solution = 1;
				exhaust_pos = row * ROWS + i;
				return;
			}
			if (tmp == 1) //only one possible is left, we consider it's the right one
					{

				(*p_matrix)[row][i] = highest_bit((*p_matrix)[row][i]);
				B_change_occur = 1;

				bit_map = 0; //update bit_map
				for (i = 0; i < COLS; i++) //recalc bit_map
						{
					if (!(((*p_matrix)[row][i]) & 0x8000000000000000)) {
						if ((bit_map | (1 << ((*p_matrix)[row][i] - 1)))
								== bit_map) {
							printf("\nelement conflict, wrong branch\n");
							B_no_solution = 1;
							return;
						}
						bit_map |= 1 << ((*p_matrix)[row][i] - 1);
					}
				}

				//else return NO_SOLUTION;

			}

			if ((*p_matrix)[row][i] == 0x8000000000000000) //no possible value left, there is no solution, the puzzle is wrong
					{
				printf("\nno value for matrix[%d][%d]\n", row, i);
				B_no_solution = 1;
				exhaust_pos = row * ROWS + i;
				return;
			}
		}
	}
	return;
}

void col_check(unsigned long (*p_matrix)[ROWS][COLS], int col) {
	int i, bit_map = 0;
	for (i = 0; i < ROWS; i++) //calc bit_map
			{
		if (!(((*p_matrix)[i][col]) & 0x8000000000000000)) {
			//int debug = (1<<((*p_matrix)[i][col]-1));int debug2 = bit_map|debug;
			if ((bit_map | (1 << ((*p_matrix)[i][col] - 1))) == bit_map) {
				printf("\nelement conflict, wrong branch\n");
				B_no_solution = 1;
				conflict_pos = i * ROWS + col;
				return;
			}
			bit_map |= 1 << ((*p_matrix)[i][col] - 1);
		}
	}

	for (i = 0; i < ROWS; i++) {
		if (((*p_matrix)[i][col]) & 0x8000000000000000) {
			if (((~bit_map) & ((*p_matrix)[i][col])) != (*p_matrix)[i][col]) {
				((*p_matrix)[i][col]) &= (~bit_map);
				//B_change_occur=1;
			}

			int tmp = bit_count((*p_matrix)[i][col]);
			if (tmp == 0) {
				B_no_solution = 1;
				exhaust_pos = i * ROWS + col;
				return;
			}
			if (tmp == 1) //only one possible is left, we consider it's the right one
					{
				(*p_matrix)[i][col] = highest_bit((*p_matrix)[i][col]);
				B_change_occur = 1;

				bit_map = 0; //update bit_map
				for (i = 0; i < ROWS; i++) //recalc bit_map
						{
					if (!(((*p_matrix)[i][col]) & 0x8000000000000000)) {
						if ((bit_map | (1 << ((*p_matrix)[i][col] - 1)))
								== bit_map) {
							printf("\nelement conflict, wrong branch\n");
							B_no_solution = 1;
							return;
						}
						bit_map |= 1 << ((*p_matrix)[i][col] - 1);
					}
				}
			}

			if ((*p_matrix)[i][col] == 0x8000000000000000) {
				printf("\nno value for matrix[%d][%d]\n", i, col);
				B_no_solution = 1;
				exhaust_pos = i * ROWS + col;
				return;
			}
		}
	}
	return;
}

void block_check(unsigned long (*p_matrix)[ROWS][COLS], int nth_block) {
	int i, j, bit_map = 0;
	int block_row = nth_block / ROOT;
	int block_cul = nth_block % ROOT;

	for (i = 0; i < ROOT; i++) //calc bit_map
			{
		for (j = 0; j < ROOT; j++) {
			if (!(((*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) & 0x8000000000000000)) {
				if ((bit_map
						| (1
								<< ((*p_matrix)[block_row * ROOT + i][block_cul * ROOT
										+ j] - 1))) == bit_map) {
					printf("\nelement conflict, wrong branch\n");
					B_no_solution = 1;
					conflict_pos = i * ROWS + j;
					return;
				}
				bit_map |= 1
						<< ((*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]
								- 1);
			}

		}

	}

	for (i = 0; i < ROOT; i++) {
		for (j = 0; j < ROOT; j++) {

			if (((*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) & 0x8000000000000000) {
				if (((~bit_map)
						& ((*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]))
						!= (*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) {
					((*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]) &=
							(~bit_map);
					//B_change_occur=1;
				}

				int tmp = bit_count(
						(*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]);
				if (tmp == 0) {
					B_no_solution = 1;
					exhaust_pos = i * ROWS + j;
					return;
				}
				if (tmp == 1) //only one possible is left, we consider it's the right one
						{
					(*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j] =
							highest_bit(
									(*p_matrix)[block_row * ROOT + i][block_cul * ROOT
											+ j]);
					B_change_occur = 1;

					bit_map = 0; //update bit_map
					for (i = 0; i < ROOT; i++) //recalc bit_map
							{
						for (j = 0; j < ROOT; j++) {
							if (!(((*p_matrix)[block_row * ROOT + i][block_cul * ROOT
									+ j]) & 0x8000000000000000)) {
								if ((bit_map
										| (1
												<< ((*p_matrix)[block_row * ROOT
														+ i][block_cul * ROOT + j]
														- 1))) == bit_map) {
									printf(
											"\nelement conflict, wrong branch\n");
									B_no_solution = 1;
									return;
								}
								bit_map |=
										1
												<< ((*p_matrix)[block_row * ROOT
														+ i][block_cul * ROOT + j]
														- 1);
							}

						}

					}
				}

				if ((*p_matrix)[block_row * ROOT + i][block_cul * ROOT + j]
						== 0x8000000000000000) {
					printf("\nno value for matrix[%d][%d]\n", block_row * ROOT + i,
							block_cul * ROOT + j);
					B_no_solution = 1;
					exhaust_pos = i * ROWS + j;
					return;
				}
			}

		}

	}
	return;
}

void divide_to_new(node** foo) {
	node *p_new_node = malloc(sizeof(struct node));
	p_new_node->current_matrix = malloc(ROWS*COLS * sizeof(unsigned long));
	memcpy(p_new_node->current_matrix, (*foo)->current_matrix,
			ROWS*COLS * sizeof(unsigned long));
	(*(p_new_node->current_matrix))[(*foo)->divide_from_x][(*foo)->divide_from_y] =
			(*foo)->current_candicate;

	p_new_node->divide_from_x = p_new_node->divide_from_y =
			p_new_node->current_candicate = 0;
	p_new_node->p_prev_node = *foo;
	simple_node_link[tree_length * 2] = (*foo)->divide_from_x * ROWS
			+ (*foo)->divide_from_y;
	simple_node_link[tree_length * 2 + 1] = (*foo)->current_candicate;
	*foo = p_new_node;

	tree_length++;
	tree_change_times++;
	printf("\ntree changes %d times\n",tree_change_times);
	divide_times++;

	printf("\ntry a new branch&add a new node\nit has %d nodes now!!!\n",
			tree_length);

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

	unsigned long (*p_matrix)[ROWS][COLS] = (*p_tree)->current_matrix;
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

void next_least_candicate_on_matrix(node * tree) {
	unsigned long * p_matrix = tree->current_matrix;

}

void choose_best_candicate(node ** foo) {
	int i, least_candicate, least_x, least_y;

	calc_least_candicate((*foo)->current_matrix, &least_x, &least_y,
			&least_candicate);

	if (least_candicate == 0xFF) {
		B_no_solution = 1; //send a branch-over signal
		exhaust_pos = 0xFF;
		return;
	}
	(*foo)->divide_from_x = least_x;
	(*foo)->divide_from_y = least_y;
	(*foo)->current_candicate = highest_bit(
			(*((*foo)->current_matrix))[least_x][least_y]);
	divide_to_new(foo);
	print_simple_node_link(simple_node_link);

}

int search_solution(unsigned long (*p_matrix)[ROWS][COLS]) {
	int i;

	B_change_occur = 1;
	B_no_solution = 0;
	int debug_count = 0;
	conflict_pos = exhaust_pos = -1;

	while (B_change_occur && !B_no_solution) //everytime we search a matrix, we must loop the three phase until no change happen
	{
		debug_count++;
		B_change_occur = 0; //B_no_solution = 0;
		for (i = 0; (i < ROWS) && !B_no_solution; i++) {
			row_check(p_matrix, i);
		}

		for (i = 0; (i < COLS) && !B_no_solution; i++) {
			col_check(p_matrix, i);
		}

		for (i = 0; (i < ROWS) && !B_no_solution; i++) {
			block_check(p_matrix, i);
		}
	}

	if (B_no_solution) {

		return NO_SOLUTION;
	}

	else {
		int tmp = calc_least_candicate(p_matrix, &i, &i, &i);
		if (tmp == 0xFF) {
			B_got_solution = 1;
			return GOT_SOLUTION;
		}
		if (tmp > 1)
			return NO_CHANGE_SO_PAUSE;
	}

}
int main(void) {
	//char * str = (char*)malloc(200);
	int i, tmp;
	tree_length = 0;
	tree_change_times=0;
	memset(simple_node_link,-1,ROWS * COLS * 2*sizeof(int));

	int *debug_tree_length = &tree_length;
	int *debug_tree_change_times = &tree_change_times;
	int *debug_divide_times = &divide_times;
	int *debug_conflict_pos = &conflict_pos;
	int *debug_exhaust_pos = &exhaust_pos;
	unsigned long (*debug_simple_node_link)[ROWS * COLS] = &simple_node_link;

	for (i = 0; i < ROWS * COLS; i++) {
		//scanf("%d",&tmp);
		tmp = matrix[i / COLS][i % COLS];
		matrix[i / COLS][i % COLS] = tmp ? tmp : 0xFFFFFFFFFFFFFFFF; //first bit means this is a solution, last 9 bits is the solution map;
		original_matrix[i / COLS][i % COLS] = tmp ? tmp : 0xFFFFFFFFFFFFFFFF;
	}

	print_matrix(matrix);

	node head_node;
	head_node.current_matrix = &matrix;
	head_node.current_candicate = head_node.divide_from_x =
			head_node.divide_from_y = 0;

	node * tree;
	tree = &head_node;

	//choose_best_candicate(tree);


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

			if (B_got_solution == 1) {
				break;
			}

			if (B_no_solution == 1) //branch is over, no need to search
					{
				result = NO_SOLUTION;
				continue;
			}

			choose_best_candicate(&tree); //!!!NOTE:tree has been changed to a new one!

			if (B_no_solution == 1) //branch is over, no need to search
					{
				result = NO_SOLUTION;
				continue;
			}

			result = search_solution(tree->current_matrix); //result only comes from search_solution
			printf("\nafter search, ");
			print_matrix(*(tree->current_matrix));
		}
	);


	if (B_got_solution == 1) {
		printf("\nGot a solution!\n");
		print_matrix(*(tree->current_matrix));

	}

	else {
		printf("No solution!\nThe last found matrix is like this\n");
		print_matrix(*(tree->current_matrix));
		print_solution(*(tree->current_matrix));

	}
	return EXIT_SUCCESS;
}
