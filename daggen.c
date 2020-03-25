/******************************************************************************
 * Copyright (c) 2007-2013. F. Suter, S. Hunold.
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL 2.1) which comes with this package.
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include "daggen_commons.h"
#include "mcts.h"
#include "heft.h"
#include "peft.h"
#include "cpop.h"
/*************************/
/** Static declarations **/
/*************************/
static DAG generateDAG(void);
static void generateTasks(DAG dag);
static void generateDependencies(DAG dag);
static void generateTransfers(DAG dag);
static double computeAvebw(DAG dag);
static double MaxSucc(DAG dag, Task task);
static void computeUpRank(DAG dag);
static double MaxPred(DAG dag, Task task);
static void computeDownRank(DAG dag);
static void computeSumRank(DAG dag);
static void computeCP(DAG dag);
static void generateCP(DAG dag);
static void freeDAG(DAG dag);

/************/
/** main() **/
/************/

int main(int argc, char **argv)
{
  time_t now = time(NULL);
  DAG dag, dag_read;
  Tree tree;
  Queue queue;
  int i,j,k,mcts_time;
  double mcts_makespan, heft_makespan, reduce, d_value, BB_mcts_makespan, peft_makespan,cpop_makespan,
         vio_makespan, SUM_CP, min_sum_W, temp_sum_W;

  srand((unsigned int)getpid() + (unsigned int)time(NULL));
  /* parse command line options */
  if (parseOptions(argc, argv) == -1)
  {
    printUsage();
    exit(1);
  }

  /* putting header */
  fprintf(OUTPUT, "// DAG automatically generated by daggen at %s// ",
          ctime(&now));
  for (i = 0; i < argc; i++)
  {
    fprintf(OUTPUT, "%s ", argv[i]);
  }
  fprintf(OUTPUT, "\n");

  FILE *fp_am;
  fp_am = fopen("nohup.out", "a+");
  init_processor();
  init_bandwidth();
  assign_processor_bandwidth();

  /* generating the dag */
  dag = generateDAG();

  /* generate output */
  if (dag)
  {
    if (global.dot_output)
      outputDOT(dag);
    else
      outputDAG(dag);  
  }

  SUM_CP = 0;
  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      SUM_CP = SUM_CP + dag->levels[i][j]->cp;
    }
  }
 // printf("sum_cp = %lf\n",SUM_CP);
  min_sum_W = 999999;
  for (k = 0; k < global.n_processors; k++)
  {
    temp_sum_W = 0;
    for (i = 0; i < dag->nb_levels; i++)
    {
      for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
      {
        temp_sum_W = temp_sum_W + dag->levels[i][j]->comp_costs[k];
      }
    }
    if (temp_sum_W < min_sum_W)
    {
      min_sum_W = temp_sum_W;
    }
  }
 //printf("min_sum_w=%lf\n",min_sum_W);

  //printf("main 1\n");
  queue = init_queue();
  //printf("main 2\n");
  tree = init_tree(dag, queue);
  clock_t start, end;
  start = clock();
  mcts_makespan = UCTsearch(dag, tree, queue);
  end = clock();
  mcts_time = (end - start);
  printf("mcts_makespan = %lf\n", mcts_makespan);
  // fprintf(fp_am, "%d %d %.1lf %lf %lf %lf %ld  ", global.n, global.n_processors, global.ccr,
  //         mcts_makespan, mcts_makespan / SUM_CP, min_sum_W / mcts_makespan, (end - start));
  //printf("UCTsearch time:%ld\n", (end - start));
  //alarm(1);
  BB_assign_processor_bandwidth();
  queue = BB_init_queue();
  tree = BB_init_tree(dag, queue);
  start = clock();
  BB_mcts_makespan = BB_UCTsearch(dag, tree, queue);
  end = clock();  
  printf("BB_makespan = %lf\n", BB_mcts_makespan);
  // fprintf(fp_am, "%lf %lf %lf %ld  ", BB_mcts_makespan, BB_mcts_makespan / SUM_CP, 
  //         min_sum_W / BB_mcts_makespan, (end - start));
  // printf("BB_makespan time:%ld\n", (end - start));
  //mcts_time = 300000;
  vio_assign_processor_bandwidth();
  start=clock();
  vio_makespan = init_DAG_mark(dag, mcts_time);
  end=clock();
  printf("vio_makespan = %lf\n", vio_makespan); 
  // printf("Viosearch time:%ld\n",(end-start));
  // fprintf(fp_am, "%lf %lf %lf %ld  ", vio_makespan, vio_makespan / SUM_CP, 
  //        min_sum_W / vio_makespan, (end - start));

  start = clock();
  heft_makespan = heft(dag);
  end = clock(); 
  printf("heft_makespan = %lf\n", heft_makespan);
  // fprintf(fp_am, "%lf %lf %lf %ld  ", heft_makespan, heft_makespan / SUM_CP, 
  //         min_sum_W / heft_makespan, (end - start));
  //printf("HEFT time:%ld\n",(end-start));
 
  start = clock();
  peft_makespan = peft(dag);
  end = clock();
  printf("peft_makespan = %lf\n", peft_makespan);
  // printf("PEFT time:%ld\n", (end - start));
  // fprintf(fp_am, "%lf %lf %lf %ld  ", peft_makespan, peft_makespan / SUM_CP, 
  //         min_sum_W / peft_makespan, (end - start));

  start = clock();
  cpop_makespan = cpop(dag);
  end = clock(); 
  printf("cpop_makespan = %lf\n", cpop_makespan);
  // printf("CPOP time:%ld\n", (end - start));
  // fprintf(fp_am, "%lf %lf %lf %ld  ", cpop_makespan, cpop_makespan / SUM_CP, 
  //         min_sum_W / cpop_makespan, (end - start));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

  fprintf(fp_am, "\n");
  fclose(fp_am);

  // FILE *fp_ms;
  // fp_ms = fopen("makespan.txt", "a+");
  // fprintf(fp_ms, "%lf\n", heft_makespan);
  // fclose(fp_ms);

  // d_value = heft_makespan - mcts_makespan;
  // reduce = 100 * d_value / heft_makespan;

  // printf("----------------\n");
  // printf("D_value = %.4lf\n", d_value);
  // printf("reduce = %.4lf%%\n", reduce);
  // printf("----------------\n");

  // FILE *info;
  // info = fopen("info.txt", "w");
  // fprintf(info, "%d %d %d %d %lf %lf %lf %lf\n", global.n, global.n_processors, global.n_simulations, global.cp, mcts_makespan, heft_makespan, d_value, reduce);
  // fclose(info);

  /* free all created data structures */
  freeDAG(dag);

  exit(0);
}

/********************/
/** DAG generation **/
/********************/

static DAG generateDAG(void)
{
  int i, j, k;
  DAG dag;

  dag = (DAG)calloc(1, sizeof(struct _DAG));

  /* Generating all the tasks */
  generateTasks(dag);

  /* Generating the Dependencies */
  generateDependencies(dag);

  /* Generating the transfer costs */
  generateTransfers(dag);

  computeDownRank(dag);

  computeUpRank(dag);

  computeSumRank(dag);

  computeCP(dag);

  generateCP(dag);

  k = 0;
  for (i = 0; i < dag->nb_levels; i++)
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      global.nb_parents[k] = dag->levels[i][j]->nb_parents;
      k++;
    }
 //print_dag_info(dag);

  return dag;
}

static void generateCP(DAG dag)
{
  int i, j, k;
  Task task;
  double min;
  //min = 999999;
  //temp=0;

  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      task = dag->levels[i][j];
      //printf("---------task_id = %lf\n",task->cost);
      if (task->flat == 1)
      {
        min = 999999;
        for (k = 0; k < global.n_processors; k++)
        {
          if (task->comp_costs[k] < min)
            min = task->comp_costs[k];
        }
        //temp += min;
        task->cp = min;
        //printf("task_id = %d , task_cp = %lf\n",task->tag,task->cp);
      }
      else
      {
        task->cp = 0;
        //printf("task_id = %d , task_cp = %lf\n",task->tag,task->cp);
      }
    }
  }
}

//find out all key node
static void computeCP(DAG dag)
{
  int i, j;
  double CP, temp, min;
  Task temp_task, task;
  CP = -999999;
  int node_count = 1;

  /* count and tag the nodes */
  for (i = 0; i < dag->nb_levels; i++)
  {
    //printf("level_id=%d\n",i);
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      dag->levels[i][j]->tag = node_count;
      //printf("%d ",dag->levels[i][j]->tag);
      node_count++;
    }
    //printf("\n");
  }
  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      task = dag->levels[i][j];
      if (task->sum_rank > CP)
      {
        CP = task->sum_rank;
      }
    }
  }
  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      task = dag->levels[i][j];
      task->flat = 0;
      if (task->sum_rank == CP)
        task->flat = 1;
    }
  }
}

static void computeSumRank(DAG dag)
{
  int i, j;
  Task task;
  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      task = dag->levels[i][j];
      task->sum_rank = task->urank + task->drank;
      printf("sum=%lf\n",task->sum_rank);
    }
  }
}

//computer up rank
static void computeUpRank(DAG dag)
{
  int i, j, k;
  Task task;
  i = dag->nb_levels - 1;
  for (k = 0; k < dag->nb_tasks_per_level[i]; k++)
  {
    task = dag->levels[i][k];
    task->urank = task->ave_comp_cost;
  }
  if (dag->nb_levels >= 2)
  {
    for (i = dag->nb_levels - 2; i >= 0; i--)
    {
      for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
      {
        task = dag->levels[i][j];
        task->urank = task->ave_comp_cost + MaxSucc(dag, task);
      }
    }
  }
}

static double MaxSucc(DAG dag, Task task)
{
  int i, j;
  double temp, max;
  Task child;
  max = 0;
  for (i = 0; i < task->nb_children; i++)
  {
    child = task->children[i];
    temp = (mcts_g.comm_costs[task->tag][child->tag] / computeAvebw(dag)) + child->urank;
    if (temp > max)
      max = temp;
  }
  return max;
}

//computer down rank
static void computeDownRank(DAG dag)
{
  int i, j, k;
  Task task;
  //double avg_bw;
  i = 0;
  for (k = 0; k < dag->nb_tasks_per_level[i]; k++)
  {
    dag->levels[i][k]->drank = 0;
  }
  if (dag->nb_levels >= 2)
  {
    for (i = 1; i < dag->nb_levels; i++)
    {
      for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
      {
        task = dag->levels[i][j];
        task->drank = MaxPred(dag, task);
      }
    }
  }
}

static double MaxPred(DAG dag, Task task)
{
  int i, j;
  double max, temp;
  Task parent, child;
  max = 0;
  temp = 0;
  for (i = 0; i < task->nb_parents; i++)
  {
    parent = task->parents[i];
    temp = (mcts_g.comm_costs[parent->tag][task->tag] / computeAvebw(dag)) + parent->drank + parent->ave_comp_cost;
    if (temp > max)
      max = temp;
  }
  return max;
}

static double computeAvebw(DAG dag)
{
  double avg_bw = 0;
  int i, j;
  for (i = 0; i < global.n_processors; i++)
    for (j = 0; j < global.n_processors; j++)
    {
      avg_bw += mcts_g.bandwidth[i][j];
    }
  avg_bw = avg_bw / (global.n_processors * global.n_processors);
  //printf("avg_bw=%lf\n",avg_bw);
  return avg_bw;
}
/*
 * generateTransfers()
 *
 * Enforces the CCR
 */
static void generateTransfers(DAG dag)
{
  int i, j, k;
  double temp;
  /* assign costs.
	 Get the data size handled by the parent
	 Compute its square (matrix #elements)
	 multiply by 8 (double -> bytes)
	 Costs are in bytes
   */
  for (i = 0; i < dag->nb_levels; i++)
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      temp = 0;
      for (k = 0; k < global.n_processors; k++)
      {
        dag->levels[i][j]->comp_costs[k] =
            dag->levels[i][j]->cost / mcts_g.processor_performance[k];
        temp += dag->levels[i][j]->comp_costs[k];
        // beta, heterogeneity factor 0.8 1.2
        // printf("%f ",dag->levels[i][j]->comp_costs[k]);
      }
      dag->levels[i][j]->ave_comp_cost = temp / global.n_processors;
      //printf("\n");
    }

  for (i = 0; i < dag->nb_levels - 1; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      for (k = 0; k < dag->levels[i][j]->nb_children; k++)
      {
        dag->levels[i][j]->comm_costs[k] = dag->levels[i][j]->data_size * getRandomNumberBetween(0.1, 1.2);
      }
    }
  }

  return;
}

/*
 * generateDependencies()
 */
static void generateDependencies(DAG dag)
{
  int i, j, k, l, m, parent_index, parent_level;
  int nb_parents;
  Task parent;
  int a[50] = {0}; // array of parent IDs
  int ij[2] = {0}; // ij[0] dag level, ij[1] level index

  /* for all levels but the last one */
  /* operate at the "parent" level   */
  for (i = 1; i < dag->nb_levels; i++)
  {

    int n_tasks_above;
    n_tasks_above = get_nb_tasks(dag, i);
    // printf("n_tasks_above level %d: %d\n", i+1, n_tasks_above);
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      /* compute how many parent the task should have,
       * at least one of course */
      nb_parents = MIN(1 + (int)getRandomNumberBetween(0.0,
                                                       global.density * (dag->nb_tasks_per_level[i - 1])),
                       dag->nb_tasks_per_level[i - 1]);
      dag->levels[i][j]->nb_parents = nb_parents;
      dag->levels[i][j]->nb_parents_r = nb_parents;
      dag->levels[i][j]->parents = (Task *)calloc(nb_parents, sizeof(Task));

      // printf("task %d parents: ", n_tasks_above + j+1);
      generate_randomID(0, n_tasks_above, nb_parents, a);

      for (k = 0; k < nb_parents; k++)
      {
        IDtoIJ(a[k], ij, dag);
        parent = dag->levels[ij[0]][ij[1]];
        dag->levels[i][j]->parents[k] = parent;
        parent->children = (Task *)realloc(parent->children,
                                           (parent->nb_children + 1) * sizeof(Task));
        parent->children[(parent->nb_children)] = dag->levels[i][j];
        (parent->nb_children)++;
      }
    }
  }

  /* Allocate memory for comm_costs and tags */
  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      dag->levels[i][j]->comm_costs = (double *)calloc(
          dag->levels[i][j]->nb_children, sizeof(double));
      dag->levels[i][j]->transfer_tags = (int *)calloc(
          dag->levels[i][j]->nb_children, sizeof(int));
      dag->levels[i][j]->comp_costs = (double *)calloc(
          global.n_processors, sizeof(double));
    }
  }
}

/*
 * generateTasks()
 */
static void generateTasks(DAG dag)
{
  int i, j;
  double integral_part;
  double op = 0;
  int nb_levels = 0;
  int *nb_tasks = NULL;
  int nb_tasks_per_level;
  int total_nb_tasks = 0;
  int tmp;

  /* compute the perfect number of tasks per levels */
  modf(exp(global.fat * log((double)global.n)), &integral_part);
  nb_tasks_per_level = (int)(integral_part);

  /* assign a number of tasks per level */
  while (1)
  {
    tmp = getIntRandomNumberAround(nb_tasks_per_level, 100.00 - 100.0 * global.regular);
    if (total_nb_tasks + tmp > global.n)
    {
      tmp = global.n - total_nb_tasks;
    }
    nb_tasks = (int *)realloc(nb_tasks, (nb_levels + 1) * sizeof(int));
    nb_tasks[nb_levels++] = tmp;
    total_nb_tasks += tmp;
    if (total_nb_tasks >= global.n)
      break;
  }

  /* Put info in the dag structure */
  dag->nb_levels = nb_levels;
  dag->levels = (Task **)calloc(dag->nb_levels, sizeof(Task *));
  dag->nb_tasks_per_level = nb_tasks;
  for (i = 0; i < dag->nb_levels; i++)
  {
    dag->levels[i] = (Task *)calloc(dag->nb_tasks_per_level[i],
                                    sizeof(Task));
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      dag->levels[i][j] = (Task)calloc(1, sizeof(struct _Task));
      /** Task cost computation                **/
      /** (1) pick a data size (in elements)   **/
      /** (2) pick a complexity                **/
      /** (3) add a factor for N_2 and N_LOG_N **/
      /** (4) multiply (1) by (2) and by (3)   **/
      /** Cost are in flops                    **/

      dag->levels[i][j]->data_size = (int)getRandomNumberBetween(
          (double)MINP * global.ccr, (double)MAXP * global.ccr);
      dag->levels[i][j]->cost = getRandomNumberBetween(MINP, MAXP);
      //printf("%f ", dag->levels[i][j]->cost);
      // op = getRandomNumberBetween(64.0, 512.0);

      // if (!global.ccr) {
      //   dag->levels[i][j]->complexity = ((int) getRandomNumberBetween(
      //       global.mindata, global.maxdata) % 3 + 1);
      // } else {
      //   dag->levels[i][j]->complexity = (int)global.ccr;
      // }

      // switch (dag->levels[i][j]->complexity) {
      // case N_2:
      //   dag->levels[i][j]->cost = (op * pow(
      //       dag->levels[i][j]->data_size, 2.0));
      //   break;
      // case N_LOG_N:
      //   dag->levels[i][j]->cost = (2 * op * pow(
      //       dag->levels[i][j]->data_size, 2.0)
      //   * (log(dag->levels[i][j]->data_size)/log(2.0)));
      //   break;
      // case N_3:
      //   dag->levels[i][j]->cost
      //   = pow(dag->levels[i][j]->data_size, 3.0);
      //   break;
      // case MIXED:
      //   fprintf(stderr, "Modulo error in complexity function\n");
      //   break;
      // }

      dag->levels[i][j]->alpha = getRandomNumberBetween(global.minalpha,
                                                        global.maxalpha);
    }
  }
}

void freeDAG(DAG dag)
{
  int i, j;
  for (i = 0; i < dag->nb_levels; i++)
  {
    for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
    {
      free(dag->levels[i][j]->transfer_tags);
      free(dag->levels[i][j]->children);
      free(dag->levels[i][j]->comm_costs);
      free(dag->levels[i][j]);
    }
    free(dag->levels[i]);
  }
  free(dag->levels);
  free(dag->nb_tasks_per_level);
  free(dag);
}