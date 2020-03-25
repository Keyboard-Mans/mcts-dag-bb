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
#include "cpop.h"

//初始化
double Init_ctask(DAG dag, Ctask *ctask){
    int i,j,k,p,q;
    double temp_comp,avg_bw;

    //cpop_g.Ptime = (double *)calloc(global.n_processors,sizeof(double));
    cpop_g.availP = (double *)calloc(global.n_processors,sizeof(double));
    cpop_g.index = (int *)calloc(global.n,sizeof(int));

    avg_bw=0;
    for(i=0;i<global.n_processors;i++)
        for(j=0;j<global.n_processors;j++){
            avg_bw += mcts_g.bandwidth[i][j];
        }
    avg_bw = avg_bw / (global.n_processors * global.n_processors);
    
    k=0;
    for(i=0;i<dag->nb_levels;i++){
        for(j=0;j<dag->nb_tasks_per_level[i];j++){
            ctask[k]->task = dag->levels[i][j];
            ctask[k]->parent_num = dag->levels[i][j]->nb_parents; 
            // temp_comp = 0;
            // for(p=0;p<global.n_processors;p++){
            //     temp_comp += dag->levels[i][j]->comp_costs[p]; 
            // }
            // ctask[k]->comp_mean = temp_comp / global.n_processors;
            ctask[k]->comp_mean = dag->levels[i][j]->ave_comp_cost;
            ctask[k]->comm_means = (double *)calloc(dag->levels[i][j]->nb_children, sizeof(double));
            ctask[k]->efts = (double *)calloc(global.n_processors,sizeof(double));
            for(p=0;p<dag->levels[i][j]->nb_children;p++){
                ctask[k]->comm_means[p] = dag->levels[i][j]->comm_costs[p] / avg_bw;
            }    
            ctask[k]->ranku = 0;
			ctask[k]->rankd = 0;
			ctask[k]->rank_sum = 0; 
			ctask[k]->flag = 0;
			ctask[k]->mark_visit = 0;       
            k++;
        }
    }
    return avg_bw;
}

double Max_succ(DAG dag, Ctask *ctask, Task task){
    int i,j;
    double temp,max;
    Task child;
    max = 0;
    for(i=0;i<task->nb_children;i++){
        child = task->children[i];
        temp = ctask[task->tag-1]->comm_means[i] + ctask[child->tag-1]->ranku;
        if(temp > max)
            max = temp;
    }
    return max;
}

void Compute_ranku(DAG dag, Ctask *ctask){
    int i,j,k;
    Task task;
    //printf("compute_ranku 1: %d\n",dag->nb_levels-1);
    i = dag->nb_levels-1;
    for(k=0;k<dag->nb_tasks_per_level[i];k++){
        task = dag->levels[i][k];
        ctask[task->tag-1]->ranku = ctask[task->tag-1]->comp_mean;
    }
    //printf("compute_ranku 2\n");
    if(dag->nb_levels >= 2){
        for(i=dag->nb_levels-2;i>=0;i--){
            for(j=0;j<dag->nb_tasks_per_level[i];j++){
                task = dag->levels[i][j];
                ctask[task->tag-1]->ranku = ctask[task->tag-1]->comp_mean 
                                            + Max_succ(dag,ctask,task);
            }
        }
    }
}

double Max_pred(DAG dag, Ctask *ctask, Task task,double ave_b){
	int i, j ;
	double max,temp;
	Task parent,child;
	max =0;
	temp=0;
	for(i=0;i<task->nb_parents;i++){
		parent = task->parents[i];
		/*for(j=0;j<parent->nb_children;j++){
		child = parent->children[j];
		if((child->tag-1)==(task->tag-1))
		printf("%d-----%d\n",child->tag-1,task->tag-1);
			k=j;	
		}*/
		//�˴����޸� 
		temp = (mcts_g.comm_costs[parent->tag][task->tag]/ave_b) + 
               ctask[parent->tag-1]->rankd + 
		       ctask[parent->tag-1]->comp_mean;
		if(temp>max)
		max = temp;
	}
	return max;
}

void Compute_rankd(DAG dag, Ctask *ctask, double ave_comn){
    int i,j,k;
    Task task;
    //double t;
    //t=ave_comn;
    //printf("compute_ranku 1: %d\n",dag->nb_levels-1);
    i = 0;  //�˴����޸� 
    for(k=0;k<dag->nb_tasks_per_level[i];k++){
        task = dag->levels[i][k];
        ctask[task->tag-1]->rankd = 0;
    }
    //printf("compute_ranku 2\n");
    if(dag->nb_levels >= 2){
        for(i=1;i<dag->nb_levels;i++){
            for(j=0;j<dag->nb_tasks_per_level[i];j++){
                task = dag->levels[i][j];
                ctask[task->tag-1]->rankd = Max_pred(dag, ctask, task, ave_comn);
                printf("rankd = %lf\n",ctask[task->tag-1]->rankd);
            }
        }
    }
}

void Compute_sumRank(Ctask *ctask){
	int i;
	for(i=0;i<global.n;i++)
    {
        ctask[i]->rank_sum = ctask[i]->ranku + ctask[i]->rankd;
        printf("CPOP_rank_sum = %lf\n",ctask[i]->rank_sum);
    }
}

int Compute_cp(DAG dag,Ctask *ctask){
	int i,j,p,k,t;
	double cp,temp,min;
	Task temp_task;
    cp =ctask[0]->rank_sum;
    for(i=0;i<global.n;i++){
        if(ctask[i]->rank_sum > cp){
            cp = ctask[i]->rank_sum;
        }
    }
	// for(p=0;p<dag->nb_tasks_per_level[0];p++){
	// 	if(ctask[p]->rank_sum>cp)
	// 	   cp = ctask[p]->rank_sum;
	// } //�˴��д��޸� 
	//printf("cp = %lf\n",cp);
	for(i=0;i<global.n;i++){
		if(ctask[i]->rank_sum == cp)
        {
            ctask[i]->flag = 1;	
            //printf("Key_node=%d\n",i);
        }
			
	}
	min = 9999999;
    k = 0;
	for(t=0;t<global.n_processors;t++){
		temp = 0;
		for(i=0;i<dag->nb_levels;i++)
		{	
			for(j=0;j<dag->nb_tasks_per_level[i];j++){  
				temp_task = dag->levels[i][j]; 
				if(ctask[temp_task->tag-1]->flag == 1)
				temp += temp_task->comp_costs[t];
			} 	
		}	
		if(temp<min){
			min = temp;
			k = t;
		}
	}
	return k;	
}

void Sort_rank(Ctask *ctask){
    int i,j;
    Ctask ptask;
    for(i=0;i<global.n-1;i++)
        for(j=0;j<global.n-i-1;j++){
            if(ctask[j]->rank_sum < ctask[j+1]->rank_sum){
                ptask = ctask[j];
                ctask[j] = ctask[j+1];
                ctask[j+1] = ptask;
            } 
        }
}


void Update_parentNum(DAG dag, Ctask *ctask, Task task){
	Task child;
	int i;
	for(i=0;i<task->nb_children;i++){
		child = task->children[i];
		//printf("---------child id %d:%d\n",ctask[cpop_g.index[child->tag-1]]->task->tag-1 ,cpop_g.index[child->tag-1]);
		if(ctask[cpop_g.index[child->tag-1]]->parent_num!=0)
			ctask[cpop_g.index[child->tag-1]]->parent_num--;
			//printf("task %d----parent_num %d\n",ctask[cpop_g.index[child->tag-1]]->task->tag-1, ctask[cpop_g.index[child->tag-1]]->parent_num);
	}
} 

double Max_availP(){
    double max=0;
    for(int i=0;i<global.n_processors;i++){
        if(max < cpop_g.availP[i])
            max = cpop_g.availP[i];
    }
    return max;
}

double Compute_est(Ctask *ctask, Task task, int processor){
    int p_id, t_id, pro_id;
    double max_pred=-1, temp_pred;
    Task parent;
    t_id = task->tag;   //mcts_g.comm_costs的下标从1开始
    max_pred = cpop_g.availP[processor];
    for(int i=0; i<task->nb_parents; i++){
        parent = task->parents[i];
        p_id = parent->tag;
        pro_id = ctask[cpop_g.index[p_id-1]]->processor;  //htask的下标从0开始
        //printf("p%d -> p%d bandwidth:%d\n",pro_id,processor,mcts_g.bandwidth[pro_id][processor]);
        //printf("t%d -> t%d comm_costs: %.2lf\n",p_id-1, t_id-1, mcts_g.comm_costs[p_id][t_id]);
        //printf("pid:%d aft:%.2f pro_id:%d\n",p_id-1, htask[heft_g.index[p_id-1]]->aft, pro_id);
        //printf("availP:%.2f\n",heft_g.availP[processor]);
        if(mcts_g.bandwidth[pro_id][processor]!=0){
            temp_pred = ctask[cpop_g.index[p_id-1]]->aft + 
                        mcts_g.comm_costs[p_id][t_id]/mcts_g.bandwidth[pro_id][processor];
            //printf("temp_pred:%.2lf\n",temp_pred);
        }else{
            temp_pred = ctask[cpop_g.index[p_id-1]]->aft;
            //printf("temp_pred:%.2lf\n",temp_pred);
        }
        if(max_pred < temp_pred)
            max_pred = temp_pred;
        //printf("max_pred:%.2lf\n",max_pred);
    }
    return max_pred;
}

double Compute_eft(DAG dag, Ctask *ctask, Task task, int processor){
    double est,eft;
    est = Compute_est(ctask,task,processor);
    eft = est + task->comp_costs[processor];
    //printf("Task %d est: %.2lf, eft: %.2lf\n",task->tag-1,est,eft);
    return eft;
}

/*double Compute_slr(DAG dag, Ctask *ctask){
	int i,j,k;
	double temp=0,min_w;
	Task task;
	for(i=0;i<dag->nb_levels;i++){	
		for(j=0;j<dag->nb_tasks_per_level[i];j++){
			task = dag->levels[i][j];
			if(ctask[task->tag-1]->flag == 1){
			min_w =INFINITY;	
			for(k=0;k<global.n_processors;k++){
				if(min_w > task->comp_costs[k]){	
					min_w = task->comp_costs[k];	
				}	
			}
			temp += min_w;
			}
		}	
	}
	return temp;		
}
*/
double cpop(DAG dag){
    int i,j,k,min_p,temp_p,T,g;
    double min_eft=INFINITY,ave_comn,SLR;
    Ctask *ctask;
    FILE *fp_cpop;
    ctask = (Ctask *)calloc(global.n,sizeof(Ctask));
    for(i=0;i<global.n;i++)
        ctask[i] = (Ctask)calloc(1,sizeof(struct _Ctask));
    //printf("heft 1: %lf\n", ctask[0]->comp_mean);
    ave_comn=Init_ctask(dag,ctask);
    //printf("heft 2\n");
    Compute_ranku(dag,ctask);
    Compute_rankd(dag,ctask,ave_comn);
    //按ranku降序调度每个任务
    Compute_sumRank(ctask);
    
    temp_p = Compute_cp(dag,ctask);
    
    Sort_rank(ctask);
    
    // for(i=0;i<global.n;i++){
    //     printf("rank=%lf\n",ctask[i]->rank_sum );
    // }

    fp_cpop = fopen("Cpop_SchedulingProcess.txt","w");
    fprintf(fp_cpop,"---------------------------CPOP---------------------------\n");
    //注意htask的下标与task->tag的转换关系，index存储该关系
    for(i=0;i<global.n;i++){
        cpop_g.index[ctask[i]->task->tag-1]=i;
        //printf("----Task Id :%d  rank_sum: %lf\n",ctask[i]->task->tag-1,ctask[i]->rank_sum);
    }
    T = global.n;
    g=1;
    while(T>0){
	//g=1;
    	for(k=0;k<global.n;k++){
    		if(g==0)
    		{
    		 k=0;	
    		}
			i=k;
    		//printf("I = %d\n", i);
    		//printf("++++++++++task id:%d \n",ctask[i]->task->tag-1);
    		if((ctask[i]->parent_num==0)&&(ctask[i]->mark_visit==0)){
    			//printf("task id:%d\n",ctask[i]->task->tag-1);
    			T--;
    			g=0;
    			ctask[i]->mark_visit = 1; 
    			Update_parentNum(dag, ctask, ctask[i]->task);
    			if(ctask[i]->flag==1){ //���Ϊ�ؼ�·����� 
    				ctask[i]->efts[temp_p] = Compute_eft(dag,ctask,ctask[i]->task,temp_p);
				    ctask[i]->processor = temp_p;
        			ctask[i]->aft = ctask[i]->efts[temp_p];
        			cpop_g.availP[temp_p] = ctask[i]->efts[temp_p];	

                // printf("Task ID:%d, processor ID:%d, EST:%lf, EFT:%lf\n",
	            //             ctask[i]->task->tag-1,
			    //             ctask[i]->processor,
			    //             ctask[i]->aft - ctask[i]->task->comp_costs[ctask[i]->processor],
			    //             ctask[i]->aft);


        			fprintf(fp_cpop,"Task ID:%d, processor ID:%d, EST:%lf, EFT:%lf\n",
	                        ctask[i]->task->tag-1,
			                ctask[i]->processor,
			                ctask[i]->aft - ctask[i]->task->comp_costs[ctask[i]->processor],
			                ctask[i]->aft);
			                
    			}
    			else{
    				min_eft=INFINITY;
        			//printf("$$$$$$task %d ranku:%.2lf\n",htask[i]->task->tag-1,htask[i]->ranku);
        			for(j=0;j<global.n_processors;j++){
            			ctask[i]->efts[j] = Compute_eft(dag,ctask,ctask[i]->task,j);
            			if(min_eft > ctask[i]->efts[j]){
                			min_eft = ctask[i]->efts[j];
                			min_p = j;
            			}
        			}
        			//printf("****task %d p:%d aft:%.2lf\n",htask[i]->task->tag-1,min_p,min_eft);
        			ctask[i]->processor = min_p;
        			ctask[i]->aft = min_eft;
        			cpop_g.availP[min_p] = min_eft;
        			//printf("****task %d p:%d aft:%.2lf\n",htask[i]->task->tag-1,htask[i]->processor,htask[i]->aft);
    			    
                    
                    
                    // printf("Task ID:%d, processor ID:%d, EST:%lf, EFT:%lf\n",
	                //         ctask[i]->task->tag-1,
			        //         ctask[i]->processor,
			        //         ctask[i]->aft - ctask[i]->task->comp_costs[ctask[i]->processor],
			        //         ctask[i]->aft);
                    
                    
                    fprintf(fp_cpop,"Task ID:%d, processor ID:%d, EST:%lf, EFT:%lf\n",
	                        ctask[i]->task->tag-1,
			                ctask[i]->processor,
			                ctask[i]->aft - ctask[i]->task->comp_costs[ctask[i]->processor],
			                ctask[i]->aft);
				}
    		//i=0;	
    		} 
    		else{
    			g=1;
    			//continue;
    		}
        }
    }
    //printf("CPOP makespan: %lf\n", Max_availP());
    return  Max_availP();
    fprintf(fp_cpop,"cpop_makespan = %lf\n",Max_availP());
    //for(i=0;i<global.n;i++)
	fprintf(fp_cpop,"\n"); 
	fclose(fp_cpop);
    //SLR = Compute_slr(dag,ctask);
    //for(i=0;i<global.n_processors;i++)
    //  printf("processor i:%lf\n",heft_g.availP[i]);
}
