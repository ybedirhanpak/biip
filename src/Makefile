INC      =  $(GUROBI_HOME)/include/
CC       =  gcc 
CPP      =  g++ 
CARGS    = -m64 -g
CLIB     = -L$(GUROBI_HOME)/lib/ -lgurobi65
CPPLIB   = -L$(GUROBI_HOME)/lib/ -lgurobi_c++ -lgurobi65

all:	biip


biip:	biip.cpp
		$(CPP) $(CARGS) -o biip biip.cpp -I$(INC) $(CPPLIB) -lpthread -lm
