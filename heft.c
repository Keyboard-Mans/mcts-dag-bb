#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "daggen_commons.h"
#include "mcts.h"
#include "heft.h"

void init_htask(DAG dag, Htask *htask){
    int i,j,k,p,q;
    double temp_comp,avg_bw;

    heft_g.Ptime = (double *)calloc(global.n_processors,sizeof(double));
    heft_g.availP = (double *)calloc(global.n_processors,sizeof(double));
    heft_g.index = (int *)calloc(global.n,sizeof(int));
    avg_bw=0;
    for(i=0;i<global.n_processors;i++)
        for(j=0;j<global.n_processors;j++){
            avg_bw += mcts_g.bandwidth[i][j];
        }
    avg_bw = avg_bw / (global.n_processors * global.n_processors);
    k=0;
    for(i=0;i<dag->nb_levels;i++){
        for(j=0;j<dag->nb_tasks_per_level[i];j++){
            htask[k]->task = dag->levels[i][j];
            temp_comp = 0;
            for(p=0;p<global.n_processors;p++){
                temp_comp += dag->levels[i][j]->comp_costs[p]; 
            }
            htask[k]->comp_mean = temp_comp / global.n_processors;
            htask[k]->comm_means = (double *)calloc(
                                dag->levels[i][j]->nb_children, sizeof(double));
            htask[k]->efts = (double *)calloc(global.n_processors,sizeof(double));
            for(p=0;p<dag->levels[i][j]->nb_children;p++){
                htask[k]->comm_means[p] = dag->levels[i][j]->comm_costs[p] / avg_bw;
            }    
            htask[k]->ranku = 0;        
            k++;
        }
    }
}

double max_succ(DAG dag, Htask *htask, Task task){
    int i,j;
    double temp,max;
    Task child;
    max = 0;
    for(i=0;i<task->nb_children;i++){
        child = task->children[i];
        temp = htask[task->tag-1]->comm_means[i] + htask[child->tag-1]->ranku;
        if(temp > max)
            max = temp;
    }
    return max;
}

void compute_ranku(DAG dag, Htask *htask){
    int i,j,k;
    Task task;
    //printf("compute_ranku 1: %d\n",dag->nb_levels-1);
    i = dag->nb_levels-1;
    for(k=0;k<dag->nb_tasks_per_level[i];k++){
        task = dag->levels[i][k];
        htask[task->tag-1]->ranku = htask[task->tag-1]->comp_mean;
    }
    //printf("compute_ranku 2\n");
    if(dag->nb_levels >= 2){
        for(i=dag->nb_levels-2;i>=0;i--){
            for(j=0;j<dag->nb_tasks_per_level[i];j++){
                task = dag->levels[i][j];
                htask[task->tag-1]->ranku = htask[task->tag-1]->comp_mean 
                                            + max_succ(dag,htask,task);
            }
        }
    }
}
void sort_ranku(Htask *htask){
    int i,j;
    Htask ptask;
    for(i=0;i<global.n-1;i++)
        for(j=0;j<global.n-i-1;j++){
            if(htask[j]->ranku < htask[j+1]->ranku){
                ptask = htask[j];
                htask[j] = htask[j+1];
                htask[j+1] = ptask;
            } 
        }
}

double max_availP(){
    double max=0;
    for(int i=0;i<global.n_processors;i++){
        if(max < heft_g.availP[i])
            max = heft_g.availP[i];
    }
    return max;
}

double compute_est(Htask *htask, Task task, int processor){
    int p_id, t_id, pro_id;
    double max_pred=-1, temp_pred;
    Task parent;
    t_id = task->tag;   //mcts_g.comm_costs的下标从1开始
    max_pred = heft_g.availP[processor];
    for(int i=0; i<task->nb_parents; i++){
        parent = task->parents[i];
        p_id = parent->tag;
        pro_id = htask[heft_g.index[p_id-1]]->processor;  //htask的下标从0开始
        //printf("p%d -> p%d bandwidth:%d\n",pro_id,processor,mcts_g.bandwidth[pro_id][processor]);
        //printf("t%d -> t%d comm_costs: %.2lf\n",p_id-1, t_id-1, mcts_g.comm_costs[p_id][t_id]);
        //printf("pid:%d aft:%.2f pro_id:%d\n",p_id-1, htask[heft_g.index[p_id-1]]->aft, pro_id);
        //printf("availP:%.2f\n",heft_g.availP[processor]);
        if(mcts_g.bandwidth[pro_id][processor]!=0){
            temp_pred = htask[heft_g.index[p_id-1]]->aft + 
                        mcts_g.comm_costs[p_id][t_id]/mcts_g.bandwidth[pro_id][processor];
            //printf("temp_pred:%.2lf\n",temp_pred);
        }else{
            temp_pred = htask[heft_g.index[p_id-1]]->aft;
            //printf("temp_pred:%.2lf\n",temp_pred);
        }
        if(max_pred < temp_pred)
            max_pred = temp_pred;
        //printf("max_pred:%.2lf\n",max_pred);
    }
    return max_pred;
}
double compute_eft(DAG dag, Htask *htask, Task task, int processor){
    double est,eft;
    est = compute_est(htask,task,processor);
    eft = est + task->comp_costs[processor];
    //printf("Task %d est: %.2lf, eft: %.2lf\n",task->tag-1,est,eft);
    return eft;
}
double heft(DAG dag){
    int i,j,k,min_p;
    double min_eft=INFINITY, heft_makespan;
    Htask *htask;
    FILE *fp_heft;
    htask = (Htask *)calloc(global.n,sizeof(Htask));
    for(i=0;i<global.n;i++)
        htask[i] = (Htask)calloc(1,sizeof(struct _Htask));
    //printf("heft 1: %lf\n", htask[0]->comp_mean);
    init_htask(dag,htask);
    //printf("heft 2\n");
    compute_ranku(dag,htask);
    //按ranku降序调度每个任务
    sort_ranku(htask);
    
    //注意htask的下标与task->tag的转换关系，index存储该关系
    for(i=0;i<global.n;i++){
        heft_g.index[htask[i]->task->tag-1]=i;
    }
    for(i=0;i<global.n;i++){
        min_eft=INFINITY;
        //printf("$$$$$$task %d ranku:%.2lf\n",htask[i]->task->tag-1,htask[i]->ranku);
        for(j=0;j<global.n_processors;j++){
            htask[i]->efts[j] = compute_eft(dag,htask,htask[i]->task,j);
            if(min_eft > htask[i]->efts[j]){
                min_eft = htask[i]->efts[j];
                min_p = j;
            }
        }
        //printf("****task %d p:%d aft:%.2lf\n",htask[i]->task->tag-1,min_p,min_eft);
        htask[i]->processor = min_p;
        htask[i]->aft = min_eft;
        heft_g.availP[min_p] = min_eft;
        //printf("****task %d p:%d aft:%.2lf\n",htask[i]->task->tag-1,htask[i]->processor,htask[i]->aft);
    }
    heft_makespan = max_availP();
    
    fp_heft = fopen("Heft_SchedulingProcess.txt","w");
    fprintf(fp_heft,"---------------------------HEFT---------------------------\n");
    fprintf(fp_heft,"heft_makespan = %lf\n",heft_makespan);
    for(i=0;i<global.n;i++)
	fprintf(fp_heft,"Task ID:%d, processor ID:%d, EST:%lf, EFT:%lf\n",
	        htask[i]->task->tag-1,
			htask[i]->processor,
			htask[i]->aft - htask[i]->task->comp_costs[htask[i]->processor],
			htask[i]->aft);
	fprintf(fp_heft,"\n"); 
	fclose(fp_heft);
    //printf("HEFT makespan: %lf\n", heft_makespan);
    return heft_makespan;
    //for(i=0;i<global.n_processors;i++)
      //  printf("processor i:%lf\n",heft_g.availP[i]);
}
