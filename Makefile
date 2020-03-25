all: daggen

daggen_OBJECTS = daggen.o daggen_commons.o mcts.o mcts_bb.o violence.o heft.o peft.o cpop.o 
# Tool invocations
daggen: $(daggen_OBJECTS) 
	gcc $(CFLAGS) -o daggen $(daggen_OBJECTS) -lm

daggen.o:daggen.c
	gcc -c daggen.c
daggen_commons.o:daggen_commons.c
	gcc -c daggen_commons.c
mcts.o:mcts.c
	gcc -c mcts.c daggen_commons.c
mcts_bb.o:mcts_bb.c
	gcc -c mcts_bb.c mcts.c daggen_commons.c
violence.o:violence.c
	gcc -c violence.c mcts.c daggen_commons.c
heft.o:heft.c
	gcc -c heft.c mcts.c daggen_commons.c
peft.o:peft.c
	gcc -c peft.c mcts.c daggen_commons.c
cpop.o:cpop.c
	gcc -c cpop.c mcts.c daggen_commons.c
.PHONY:clean
	rm *.o
