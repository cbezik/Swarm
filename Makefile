CMPLR = g++
CFLGS = -O3 -w -std=c++11 -g

EXECU = ./swarm

SRCES = swarm.cpp \
		main.cpp \

OBJCS := $(addsuffix .o, $(basename $(SRCES)))

all: ${EXECU}

%.o : %.cpp
	${CMPLR} ${CFLGS} -c $< -o $@

%.o : %.c
	${CMPLR} ${CFLGS} -c $< -o $@

$(EXECU): ${OBJCS}
	${CMPLR} ${CFLGS} ${OBJCS} -o $@

clean:
	rm -f ${OBJCS} swarm
