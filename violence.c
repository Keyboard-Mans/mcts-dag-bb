#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <signal.h>

#include "daggen_commons.h"
#include "mcts.h"

clock_t start, end;
int count = 0;
int flag = 0;

void handler()
{
    printf("hello\n");
    //exit(0);
    //return 1;
}

void vio_assign_processor_bandwidth()
{
    int i, j, k;
    mcts_g.edge_mark = (int **)calloc(global.n, sizeof(int *));
    for (i = 0; i < global.n; i++)
        mcts_g.edge_mark[i] = (int *)calloc(global.n, sizeof(int));

    mcts_g.visit_mark = (int *)calloc(global.n, sizeof(int));

    mcts_g.ans = (int *)calloc(global.n, sizeof(int));

    mcts_g.schedule_Task = (Task *)calloc(global.n + 1, sizeof(Task));

    mcts_g.makespan = 0;

    mcts_g.scheduled = (Scheduled *)calloc(global.n + 1, sizeof(Scheduled));
    mcts_g.availP = (double *)calloc(global.n_processors, sizeof(double));

    for (i = 0; i < global.n + 1; i++)
    {
        mcts_g.scheduled[i] = (Scheduled)calloc(1, sizeof(struct _scheduled));
    }
}

//求开始时间
double vio_get_start_time(DAG dag, Task task, int processor)
{
    int task_id, i, j, parent_id, parent_proc;
    double finish_time, comm_time, max_time = -1;
    //printf("get_start_time:1\n");
    task_id = task->tag;
    ////printf("get_start_time:2\n");
    for (i = 0; i < task->nb_parents; i++)
    {
        //printf("get_start_time:3\n");
        parent_id = task->parents[i]->tag;
        //printf("get_start_time:4\n");
        parent_proc = mcts_g.scheduled[parent_id]->processor;
        if (mcts_g.bandwidth[parent_proc][processor] != 0)
        {
            comm_time = mcts_g.comm_costs[parent_id][task_id] / mcts_g.bandwidth[parent_proc][processor];
        }
        else
        {
            comm_time = 0;
        }
        //printf("get_start_time:5\n");
        finish_time = mcts_g.scheduled[parent_id]->AFT + comm_time;
        //printf("get_start_time:6\n");
        if (finish_time > max_time)
        {
            max_time = finish_time;
        }
    }
    //printf("get_start_time:7\n");
    //没有父节点
    return max_time > mcts_g.availP[processor] ? max_time : mcts_g.availP[processor];
}

void vio_update_scheduled_task(DAG dag, Task task, int processor)
{
    double start_time;
    //FILE *fp_process;
    //fp_process = fopen("makes.txt","w");
    //mcts_g.count++;
    //printf("vio_update_scheduled_task:1\n");
    // printf("44444444\n");
    start_time = vio_get_start_time(dag, task, processor);
    //printf("vio_update_scheduled_task:1-1\n");
    //printf("vio_update_scheduled_task:2 %lf\n",mcts_g.scheduled[task->tag]->AFT);
    mcts_g.availP[processor] = start_time + task->comp_costs[processor];
    mcts_g.scheduled[task->tag]->processor = processor;
    mcts_g.scheduled[task->tag]->AFT = mcts_g.availP[processor];
}

//获得处理器
// int vio_get_best_processor(DAG dag, Task task)
// {
//     int bestp;
//     //srand((int)time(0));
//     //bestp = random(global.n_processors);
//     bestp = rand() % global.n_processors;
//     printf("bestp = %d\n",bestp);
//     return bestp;
// }

//获得最大的makespan
double vio_get_makespan()
{
    int i;
    double aft, max = -1;
    for (i = 0; i < global.n_processors; i++)
    {
        aft = mcts_g.availP[i];
        if (aft > max)
        {
            max = aft;
        }
    }
    return max;
    //可以记录最大的makespan值时的路径
}

int ok(int i, int cnt) //如果在ans[0,cnt-1]出现了一个本应在i后面才出现的字母,那么返回false
{
    for (int j = 0; j < cnt; j++)
        if (mcts_g.edge_mark[i][mcts_g.ans[j]])
            return 0;
    return 1;
}
void dfs(int cnt, DAG dag, double *min_makespan, int mcts_time)
//void dfs(int cnt, DAG dag)
{
    //signal(SIGALRM, handler);
    Task task1;
    double temp_makespan;
    //printf("22222222\n");
    //;
    if (flag == 1)
    {
        return;
    }
    if (cnt == global.n)
    {
        // printf("22222222\n");
        //初始化scheduled
        for (int i = 0; i < global.n_processors; i++)
        {
            mcts_g.availP[i] = 0;
            //printf("33333333\n");
        }
        for (int i = 0; i < global.n + 1; i++)
        {
            mcts_g.scheduled[i]->processor = 0;
            mcts_g.scheduled[i]->AFT = 0;
            //printf("33333333\n");
        }

        for (int i = 0; i < global.n; i++)
        {
            //printf("33333333\n");
            //printf("%d ", mcts_g.ans[i]);
            task1 = mcts_g.schedule_Task[mcts_g.ans[i]];
            //     //开始调度
            //printf("task_id:%d ", task1->tag);
            //     printf("33333333\n");
            //     //此处有问题
            vio_update_scheduled_task(dag, task1, rand() % global.n_processors);
        }
        //printf("\n");
        // count ++;
        // alarm(4);
        //for(int i=0;i<=1;i++){
        // if(count==global.n_simulations){
        end = clock();
        if ((end - start) > mcts_time)
        {
            //exit(0);
            flag = 1;
            return;
        }
        //
        //     //printf("end_time = %ld\n",end);
        //     if()
        //         return;

        // }
        //count++;
        //}

        temp_makespan = vio_get_makespan();
        //printf("makespan = %lf\n",temp_makespan);
        if (*min_makespan > temp_makespan)
        {
            *min_makespan = temp_makespan;
        }
    }

    else
    {
        //printf("#########\n");
        for (int i = 0; i < global.n; i++)
        {
            if (!mcts_g.visit_mark[i] && ok(i, cnt))
            {
                // printf("******\n");
                mcts_g.visit_mark[i] = 1;
                mcts_g.ans[cnt] = i;
                //printf("$$$$$$$$$$\n");
                dfs(cnt + 1, dag, min_makespan, mcts_time);
                //dfs(cnt + 1, dag);
                // printf("8888888888\n");
                mcts_g.visit_mark[i] = 0;
            }
        }
    }
}

double init_DAG_mark(DAG dag, int mcts_time)
{
    int i, j, k, p_id, c_id, t;
    double minmakespan;
    //signal(SIGALRM, handler);
    //Task *taskReady;
    //taskReady = (Task *)calloc(golbal.n, sizeof(Task));
    //printf("22222222\n");
    start = clock();
    t = 0;
    for (i = 0; i < dag->nb_levels; i++)
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            //taskReady
            mcts_g.schedule_Task[t] = dag->levels[i][j];
            t++;
            p_id = dag->levels[i][j]->tag - 1;
            //dag->levels[i][j]-> node_mark = true;
            for (k = 0; k < dag->levels[i][j]->nb_children; k++)
            {
                c_id = dag->levels[i][j]->children[k]->tag - 1;
                //printf("init_comm_costs:1 %d %d\n", p_id, c_id);
                //printf("init_comm_costs:2 %d\n", mcts_g.comm_costs[p_id][c_id]);
                mcts_g.edge_mark[p_id][c_id] = 1;
            }
        }
    //printf("22222222\n");
    minmakespan = 999999;

    //alarm(2);
    //  for(i = 1; i < 7; i++)
    //   {
    // printf("sleep %d ...\n", i);
    //sleep(1);
    //}

    dfs(0, dag, &minmakespan, mcts_time);
    //printf("count = %d\n",count);
    //dfs(0, dag);
    return minmakespan;
}