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
#include "peft.h"


void init_ptask(DAG dag, Ptask *ptask){
    int i,j,k,p,q;
    double temp_comp,avg_bw;

    peft_g.availP = (double *)calloc(global.n_processors,sizeof(double));
    peft_g.index = (int *)calloc(global.n,sizeof(int));
    
    avg_bw=0;
    for(i=0;i<global.n_processors;i++)
        for(j=0;j<global.n_processors;j++){
            avg_bw += mcts_g.bandwidth[i][j];
        }
    avg_bw = avg_bw / (global.n_processors * global.n_processors);
    //����ƽ�����㿪�� 
    k=0;
    for(i=0;i<dag->nb_levels;i++){
        for(j=0;j<dag->nb_tasks_per_level[i];j++){
            ptask[k]->task = dag->levels[i][j];
            //printf("task id :%d\n",ptask[k]->task->tag ); //tag��1��ʼ�� 
            temp_comp = 0;
            for(p=0;p<global.n_processors;p++){
                temp_comp += dag->levels[i][j]->comp_costs[p]; 
            }
            ptask[k]->comp_mean = temp_comp / global.n_processors;
            //����ƽ��ͨ�ſ��� 
            ptask[k]->comm_means = (double *)calloc(
                                dag->levels[i][j]->nb_children, sizeof(double));
            ptask[k]->efts = (double *)calloc(global.n_processors,sizeof(double));
            ptask[k]->oct = (double *)calloc(global.n_processors,sizeof(double));
            for(p=0;p<dag->levels[i][j]->nb_children;p++){
                ptask[k]->comm_means[p] = dag->levels[i][j]->comm_costs[p] / avg_bw;
            }    
            ptask[k]->ave_oct = 0;        
            k++;
        }
    }
}

double compute_oct(DAG dag, Ptask *ptask, Task task, int processor){
	int i,j,k,t;
	double temp,min=INFINITY,max=0;
	Task child;
    //max = 0;
    //t=0;
    for (i = 0; i < task->nb_children; i++){
    child = task->children[i];
    min=INFINITY;
    for (k = 0; k < global.n_processors ; k++)
    {
    	if(mcts_g.bandwidth[processor][k]!=0){
    		temp = ptask[child->tag-1]->oct[k] + child->comp_costs[k] + ptask[task->tag-1]->comm_means[i];
    		//printf("--------%d %lf %lf\n",t,temp,child->comp_costs[k] );
        }
        else{
            temp = ptask[child->tag-1]->oct[k] + child->comp_costs[k];
            //printf("--------%d %lf %lf\n",t,temp,child->comp_costs[k] );
        }
        if(temp<min)
        	min = temp;
        //t++;
    }
    if(max<min)
    	max = min;
    	//printf("max=%lf\n",max);
    }
    return max;
}

void compute_Rankoct(DAG dag,Ptask *ptask){
	int i,j,k;
    Task task;
    double temp;
    i = dag->nb_levels-1;
    //printf("####%d\n",i);
    for(k=0;k<dag->nb_tasks_per_level[i];k++){
    	task = dag->levels[i][k];
    	//printf("$$$$$%d\n",task->tag-1);
    	for(j=0;j<global.n_processors;j++)
    	    ptask[task->tag-1]->oct[j] = 0;
        ptask[task->tag-1]->ave_oct = 0;
        //printf("Task ID:%d  ave_oct----%lf\n",task->tag-1,ptask[task->tag-1]->ave_oct);
    }
    if(dag->nb_levels >= 2){
        for(i=dag->nb_levels-2;i>=0;i--){
            for(j=0;j<dag->nb_tasks_per_level[i];j++){
            	task = dag->levels[i][j];
            	temp = 0;
            	for(k=0;k<global.n_processors;k++){
            		ptask[task->tag-1]->oct[k] = compute_oct(dag,ptask,task,k);
            		//printf("ptask_oct %d:%lf\n",k,ptask[task->tag-1]->oct[k]);
            		temp+=ptask[task->tag-1]->oct[k];
            	}
            	ptask[task->tag-1]->ave_oct = temp/global.n_processors; 
            	//printf("Task ID:%d  ave_oct----%lf\n",task->tag-1,ptask[task->tag-1]->ave_oct);
            }
        }      
    }
}

void sort_rankoct(Ptask *ptask){
    int i,j;
    Ptask temp_task;
    for(i=0;i<global.n-1;i++)
        for(j=0;j<global.n-i-1;j++){
            if(ptask[j]->ave_oct < ptask[j+1]->ave_oct){
                temp_task = ptask[j];
                ptask[j] = ptask[j+1];
                ptask[j+1] = temp_task;
            } 
        }
}
//准备工作已完成

/*double Max_availP(){
    double max=0;
    for(int i=0;i<global.n_processors;i++){
        if(max < peft_g.availP[i])
            max = peft_g.availP[i];
    }
    return max;
}*/

double compute_EST(Ptask *ptask, Task task, int processor){
    int p_id, t_id, pro_id;
    double max_pred=-1, temp_pred;
    Task parent;
    t_id = task->tag;   //mcts_g.comm_costs的下标从1开始
    max_pred = peft_g.availP[processor];
    for(int i=0; i<task->nb_parents; i++){
        parent = task->parents[i];
        p_id = parent->tag;
        pro_id = ptask[peft_g.index[p_id-1]]->processor;  //htask的下标从0开始
        //printf("p%d -> p%d bandwidth:%d\n",pro_id,processor,mcts_g.bandwidth[pro_id][processor]);
       // printf("t%d -> t%d comm_costs: %.2lf\n",p_id-1, t_id-1, mcts_g.comm_costs[p_id][t_id]);
       // printf("pid:%d aft:%.2f pro_id:%d\n",p_id-1, ptask[peft_g.index[p_id-1]]->aft, pro_id);
       // printf("availP:%.2f\n",peft_g.availP[processor]);
        if(mcts_g.bandwidth[pro_id][processor]!=0){
            temp_pred = ptask[peft_g.index[p_id-1]]->aft + 
                        mcts_g.comm_costs[p_id][t_id]/mcts_g.bandwidth[pro_id][processor];
            //printf("temp_pred:%.2lf\n",temp_pred);
        }else{
            temp_pred = ptask[peft_g.index[p_id-1]]->aft;
            //printf("temp_pred:%.2lf\n",temp_pred);
        }
        if(max_pred < temp_pred)
            max_pred = temp_pred;
        //printf("max_pred:%.2lf\n",max_pred);
    }
    return max_pred;
}
double compute_EFT(DAG dag, Ptask *ptask, Task task, int processor){
    double est,eft;
    est = compute_EST(ptask,task,processor);
    eft = est + task->comp_costs[processor];
    //printf("Task %d est: %.2lf, eft: %.2lf\n",task->tag-1,est,eft);
    return eft;
}

double peft(DAG dag){
    int i,j,k,min_p;
    double min_oeft=INFINITY,temp_oeft,min_eft;
    Ptask *ptask;
    FILE *fp_peft;
    ptask = (Ptask *)calloc(global.n,sizeof(Ptask));
    for(i=0;i<global.n;i++)
        ptask[i] = (Ptask)calloc(1,sizeof(struct _Ptask));
    //printf("heft 1: %lf\n", htask[0]->comp_mean);
    init_ptask(dag,ptask);
    //printf("peft 2\n");
    compute_Rankoct(dag,ptask);
    //按ranku降序调度每个任务
    sort_rankoct(ptask);
    
    //注意htask的下标与task->tag的转换关系，index存储该关系
    for(i=0;i<global.n;i++){
        peft_g.index[ptask[i]->task->tag-1]=i;
        //printf("ptask[i]->task->tag-1 %d:%d\n",i,ptask[i]->task->tag-1);
        //printf("ptask[i]->ave_oct :%lf\n",ptask[i]->ave_oct);
    }
    
    for(i=0;i<global.n;i++){
    	//printf("Task id : %d\n", ptask[i]->task->tag-1);
        min_oeft=INFINITY;
        min_eft = 0;
        min_p = 0;
        //printf("$$$$$$task %d ranku:%.2lf\n",htask[i]->task->tag-1,htask[i]-k.
         for(j=0;j<global.n_processors;j++){
            ptask[i]->efts[j] = compute_EFT(dag,ptask,ptask[i]->task,j);
            //printf("EST=%lf\n",ptask[i]->efts[j]);
            temp_oeft = ptask[i]->efts[j] + ptask[i]->oct[j];
            if(min_oeft > temp_oeft){
                min_oeft = temp_oeft;
                min_p = j;
                min_eft = ptask[i]->efts[j];
            }
        }
        //printf("****task %d p:%d aft:%.2lf\n",htask[i]->task->tag-1,min_p,min_eft);
        ptask[i]->processor = min_p;
        ptask[i]->aft = min_eft;
        peft_g.availP[min_p] = min_eft;
    }
    //printf("PEFT makespan: %lf\n", peft_g.availP[min_p]);
    return peft_g.availP[min_p];
    
    fp_peft = fopen("Peft_SchedulingProcess.txt","w");
    fprintf(fp_peft,"-------------------------PEFT-------------------------\n");
    fprintf(fp_peft,"peft_makespan = %lf\n",peft_g.availP[min_p]);
    for(i=0;i<global.n;i++)
	fprintf(fp_peft,"Task ID:%d, processor ID:%d, EST:%lf, EFT:%lf\n",
	        ptask[i]->task->tag-1,
			ptask[i]->processor,
			ptask[i]->aft - ptask[i]->task->comp_costs[ptask[i]->processor],
			ptask[i]->aft);
	fprintf(fp_peft,"\n"); 
	fclose(fp_peft);
	
    //for(i=0;i<global.n_processors;i++)
      //  printf("processor i:%lf\n",heft_g.availP[i]);
}
