#ifndef CPOP_H_
#define CPOP_H_

typedef struct _Ctask *Ctask;

typedef struct{
    double *Ptime;         //处理器时间（暂未使用）
    double *availP;        //处理器就绪时间，计算eft使用
    int *index; 
}CPOP_g;
CPOP_g cpop_g;

struct _Ctask{
    Task task;
    double comp_mean;      //任务平均计算开销
    double *comm_means;    //边平均通信开销
    double ranku;          //upward rank
    double rankd;          //down rank 
    double rank_sum;       //the sum of  upward rank and down rank
    int flag;              //��ǵ���cpֵ�������� 
    int parent_num;
	int  mark_visit;       //��¼�����ʹ� 
    double aft;            //actual finish time
    double *efts;          //任务在每个处理器上的最早开始时间
    int processor;         //任务在哪个处理器上运行
};

double cpop(DAG dag);

#endif /*CPOP_H_*/
