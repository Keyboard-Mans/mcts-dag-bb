#ifndef MCTS_H_
#define MCTS_H_

#define MINP 50
#define MAXP 300
#define MAXtask 500
//#define fp_process (mcts_g.output_process)
typedef struct _Node *Node;
typedef struct _Tree *Tree;
typedef struct _Queue *Queue;
typedef struct _Element *Element;
typedef struct _scheduled *Scheduled;
//typedef struct _TempScheduled *TempScheduled;

typedef struct{
    int *processor_performance;
    int **bandwidth;
    double **comm_costs;
    double makespan;
    int node_id;
    int nb_node;
    Node *nodetrace;
    Scheduled *scheduled;
    double *availP;
    int count;
    int *ans;
    int *visit_mark;
    Task *schedule_Task;
    int **edge_mark;
    //int count; 
    //FILE *output_process;
//    int **comp_costs;
//    int **comm_costs;
} MCTS_g;
MCTS_g mcts_g;

struct _scheduled{
    int processor;
    double AFT;
};

struct _Element{
    Task task;
    Element next;
};

struct _Queue{
    int n;
    Element head;
};

struct _Node{
    int node_id;
    int task_id;
    Task task;
    int processor;
    double Q;
    int visits;
    int nb_children;
    Node *children;  
	int marker;  
};

struct _Tree{
    int depth;
    Node root;
};

void create_root_node(DAG dag, Tree tree, Queue queue);
void assign_processor_bandwidth();
//void Assign_processor_bandwidth();
void BB_assign_processor_bandwidth();
void vio_assign_processor_bandwidth(); 
double UCTsearch(DAG dag, Tree tree, Queue queue);
double BB_UCTsearch(DAG dag, Tree tree, Queue queue);
double init_DAG_mark(DAG dag, int mcts_time); 
Tree init_tree(DAG dag, Queue queue);
Tree BB_init_tree(DAG dag, Queue queue);
Queue init_queue();
Queue BB_init_queue();
void init_processor();
//void BB_init_processor();
void init_bandwidth();
//void BB_init_processor();

#endif /*MCTS_H_*/
