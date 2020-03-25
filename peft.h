#ifndef PEFT_H_
#define PEFT_H_

typedef struct _Ptask *Ptask;

typedef struct{
    double *Ptime;         //处理器时间（暂未使用）
    double *availP;        //处理器就绪时间，计算eft使用
    int *index; 
}PEFT_g;
PEFT_g peft_g;

struct _Ptask{
    Task task;
    double comp_mean;      //任务平均计算开销
    double *comm_means;    //边平均通信开销
    double ranku;          //upward rank
    double aft;            //actual finish time
    double *efts;          //任务在每个处理器上的最早开始时间
    int processor;         //任务在哪个处理器上运行
    double *oct;
    double ave_oct;
    double oeft; 
};

double peft(DAG dag);

#endif /*PEFT_H_*/
